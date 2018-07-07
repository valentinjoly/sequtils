import multiprocessing
import os

import Bio.SeqIO

from vjoly import open_file, partition_list, merge_dictionaries, STDIO_PATH


def make_tmp_output_paths(output_path, nb_files=1):
    if output_path is None:
        return [None] * nb_files
    if output_path == STDIO_PATH:
        d, name, ext = '.', 'stdout', ''
    else:
        d, f = os.path.split(output_path)
        if '.' in f:
            name, ext = f[:f.rfind('.')], f[f.rfind('.'):]
        else:
            name, ext = f, ''
    return [os.path.join(d, ''.join([name, '.', str(i), ext]))
            for i in range(nb_files)]


def merge_outputs(tmp_output_paths, output_path):
    with open_file(output_path, 'w') as output_file:
        for tmp_output_path in tmp_output_paths:
            with open(tmp_output_path) as tmp_output_file:
                for line in tmp_output_file:
                    output_file.write(line)
    for tmp_output_path in tmp_output_paths:
        os.remove(tmp_output_path)


def dump_queue(queue):
    results_parts = []
    while not queue.empty():
        results_parts.append(queue.get())
    return results_parts


def parallelize_unindexed(
        target_function, nb_cpus, input_paths, output_paths,
        input_format, output_format, merge_func=None, **kwargs):
    if output_paths is None:
        output_paths = [None]
    if len(output_paths) == 1:
        output_path = output_paths[0]
        results = _parallelize_unindexed_merge(
            target_function, nb_cpus, input_paths, output_path,
            input_format, output_format, merge_func, **kwargs)
    else:
        results = _parallelize_unindexed_sep(
            target_function, nb_cpus, input_paths, output_paths,
            input_format, output_format, merge_func, **kwargs)
    return results


def _parallelize_unindexed_sep(
        target_function, nb_cpus, input_paths, output_paths,
        input_format, output_format, merge_func=None, **kwargs):
    input_paths_parts = partition_list(input_paths, nb_cpus)
    output_paths_parts = partition_list(output_paths, nb_cpus)
    parts = zip(input_paths_parts, output_paths_parts)
    queue = multiprocessing.Queue()
    processes = []
    for input_paths_part, output_paths_part in parts:
        processes.append(multiprocessing.Process(
            target=_parallelize_unindexed_sep_part,
            args=(target_function, input_paths_part, output_paths_part,
                  input_format, output_format, queue),
            kwargs=kwargs))
    for p in processes:
        p.start()
    for p in processes:
        p.join()

    results_parts = dump_queue(queue)
    return merge_dictionaries(results_parts, merge_func)


def _parallelize_unindexed_sep_part(
        target_function, input_paths_part, output_paths_part,
        input_format, output_format, queue, **kwargs):
    for input_path, output_path in zip(input_paths_part, output_paths_part):
        with open_file(input_path) as input_file, \
             open_file(output_path, 'w') as output_file:
            records = Bio.SeqIO.parse(input_file, input_format)
            result = target_function(
                records, output_file, output_format, **kwargs)
            queue.put(result)


def _parallelize_unindexed_merge(
        target_function, nb_cpus, input_paths, output_path,
        input_format, output_format, merge_func=None, **kwargs):
    input_paths_parts = partition_list(input_paths, nb_cpus)
    nb_parts = len(input_paths_parts)
    if nb_parts > 1:
        tmp_output_paths = make_tmp_output_paths(output_path, nb_parts)
    else:
        tmp_output_paths = [output_path]
    parts = zip(input_paths_parts, tmp_output_paths)
    queue = multiprocessing.Queue()
    processes = []
    for input_paths_part, tmp_output_path in parts:
        processes.append(multiprocessing.Process(
            target=_parallelize_unindexed_merge_part,
            args=(target_function, input_paths_part, tmp_output_path,
                  input_format, output_format, queue),
            kwargs=kwargs))
    for p in processes:
        p.start()
    for p in processes:
        p.join()

    if nb_parts > 1:
        merge_outputs(tmp_output_paths, output_path)
    results_parts = dump_queue(queue)
    return merge_dictionaries(results_parts, merge_func)


def _parallelize_unindexed_merge_part(
        target_function, input_paths_part, tmp_output_path,
        input_format, output_format, queue, **kwargs):
    with open_file(tmp_output_path, 'w') as output_file:
        for input_path in input_paths_part:
            with open_file(input_path) as input_file:
                records = Bio.SeqIO.parse(input_file, input_format)
                result = target_function(
                    records, output_file, output_format, **kwargs)
                queue.put(result)


def parallelize_indexed(
        target_function, nb_cpus, input_paths, index_paths, output_paths,
        input_format, output_format, merge_func=None, **kwargs):
    if output_paths is None:
        output_paths = [None]
    if len(output_paths) == 1:
        output_path = output_paths[0]
        results = _parallelize_indexed_merge(
            target_function, nb_cpus, input_paths, index_paths, output_path,
            input_format, output_format, merge_func, **kwargs)
    else:
        results = _parallelize_indexed_sep(
            target_function, nb_cpus, input_paths, index_paths, output_paths,
            input_format, output_format, merge_func, **kwargs)
    return results


def _parallelize_indexed_merge(
        target_function, nb_cpus, input_paths, index_paths, output_path,
        input_format, output_format, merge_func=None, **kwargs):
    if index_paths is None:
        index_paths = [None] * len(input_paths)
    index_paths_parts = partition_list(index_paths, nb_cpus)
    nb_parts = len(index_paths_parts)

    if nb_parts > 1:
        tmp_output_paths = make_tmp_output_paths(output_path, nb_parts)
    else:
        tmp_output_paths = [output_path]

    if len(index_paths) > 1:
        input_paths_parts = partition_list(input_paths, nb_cpus)
    elif len(input_paths) > 1:
        input_paths_parts = [[input_paths]]
    else:
        input_paths_parts = [input_paths]

    parts = zip(input_paths_parts, index_paths_parts, tmp_output_paths)
    queue = multiprocessing.Queue()
    processes = []
    for input_paths_part, index_paths_part, tmp_output_path in parts:
        processes.append(multiprocessing.Process(
            target=_parallelize_indexed_merge_part,
            args=(target_function, input_paths_part, index_paths_part,
                  tmp_output_path, input_format, output_format, queue),
            kwargs=kwargs))
    for p in processes:
        p.start()
    for p in processes:
        p.join()

    if nb_parts > 1:
        merge_outputs(tmp_output_paths, output_path)
    results_parts = dump_queue(queue)
    return merge_dictionaries(results_parts, merge_func)


def _parallelize_indexed_merge_part(
        target_function, input_paths_part, index_paths_part, tmp_output_path,
        input_format, output_format, queue, **kwargs):
    paths = zip(input_paths_part, index_paths_part)
    with open_file(tmp_output_path, 'w') as output_file:
        for input_path, index_path in paths:
            if index_path is None:
                records = Bio.SeqIO.index(input_path, input_format)
            else:
                records = Bio.SeqIO.index_db(
                    index_path, input_path, input_format)
            results = target_function(
                records, output_file, output_format, **kwargs)
            queue.put(results)


def _parallelize_indexed_sep(
        target_function, nb_cpus, input_paths, index_paths, output_paths,
        input_format, output_format, merge_func=None, **kwargs):
    if index_paths is not None and len(index_paths) == 1:
        input_paths_parts = [[input_paths]]
        index_paths_parts = [index_paths]
        output_paths_parts = [output_paths]
    else:
        input_paths_parts = partition_list(input_paths, nb_cpus)
        output_paths_parts = partition_list(output_paths, nb_cpus)
        if index_paths is not None:
            index_paths_parts = partition_list(index_paths, nb_cpus)
        else:
            index_paths_parts = partition_list(
                [None]*len(input_paths_parts), nb_cpus)

    parts = zip(input_paths_parts, index_paths_parts, output_paths_parts)
    queue = multiprocessing.Queue()
    processes = []
    for input_paths_part, index_paths_part, output_paths_part in parts:
        processes.append(multiprocessing.Process(
            target=_parallelize_indexed_sep_part,
            args=(target_function, input_paths_part, index_paths_part,
                  output_paths_part, input_format, output_format, queue),
            kwargs=kwargs))
    for p in processes:
        p.start()
    for p in processes:
        p.join()

    results_parts = dump_queue(queue)
    return merge_dictionaries(results_parts, merge_func)


def _parallelize_indexed_sep_part(
        target_function, input_paths_part, index_paths_part,
        output_paths_part, input_format, output_format, queue, **kwargs):
    paths = zip(input_paths_part, index_paths_part, output_paths_part)
    for input_path, index_path, output_path in paths:
        with open_file(output_path, 'w') as output_file:
            if index_path is None:
                records = Bio.SeqIO.index(input_path, input_format)
            else:
                records = Bio.SeqIO.index_db(
                    index_path, input_path, input_format)
            results = target_function(
                records, output_file, output_format, **kwargs)
            queue.put(results)
