import multiprocessing

import Bio.SeqIO

from vjoly import (check_file_arg, check_num_arg, check_seq_file_format,
                   SEQFILE_EXTS, ZIP_EXTS, SEQFILE_FORMAT, STDIO_PATH,
                   FASTA_FORMAT, GENBANK_FORMAT)

INDEX_PATHS = None
OUTPUT_PATH = None
OUTPUT_PATHS = [None]
INPUT_FORMAT = 'input_format'
OUTPUT_FORMAT = 'output_format'
CPUS = 1
ALPHABET = None


def add_seqio_arguments(parser, add_input=True, add_index=True,
                        add_output=True, multiple_outputs=True):
    if add_input:
        parser.add_argument(
            'input_paths', metavar='FILE', nargs='+',
            help='Input sequence file(s).')
    if add_index:
        parser.add_argument(
            '-i', '--index_path', metavar='FILE',
            action='append', default=INDEX_PATHS,
            help='Index file(s) for input file(s). [{}]'.format(INDEX_PATHS))
    if add_output:
        if multiple_outputs:
            parser.add_argument(
                '-o', '--output_path', metavar='FILE',
                action='append', default=OUTPUT_PATHS,
                help='Output sequence file(s). [{}]'.format(STDIO_PATH))
        else:
            parser.add_argument(
                '-o', '--output_path', metavar='FILE', default=OUTPUT_PATH,
                help='Output sequence file. [{}]'.format(STDIO_PATH))
    if add_input:
        input_format = parser.add_mutually_exclusive_group()
        input_format.add_argument(
            '-f', '--fasta_input',
            action='store_const', const=FASTA_FORMAT, dest=INPUT_FORMAT,
            help='Input is in FASTA format.')
        input_format.add_argument(
            '-g', '--genbank_input',
            action='store_const', const=GENBANK_FORMAT, dest=INPUT_FORMAT,
            help='Input is in Genbank format.')
    if add_output:
        output_format = parser.add_mutually_exclusive_group()
        output_format.add_argument(
            '-F', '--fasta_output',
            action='store_const', const=FASTA_FORMAT, dest=OUTPUT_FORMAT,
            help='Write output in FASTA format.')
        output_format.add_argument(
            '-G', '--genbank_output',
            action='store_const', const=GENBANK_FORMAT, dest=OUTPUT_FORMAT,
            help='Write output in Genbank format.')


def add_other_arguments(parser, add_cpus=True, add_help=True):
    if add_cpus:
        parser.add_argument(
            '-c', '--cpus', metavar='INT', type=int, default=CPUS,
            help='Number of CPUs allocated for the execution.')
    if add_help:
        parser.add_argument(
            '-h', '--help', action='help',
            help='Print this help page and exit.')


def check_inputs(args, ext=SEQFILE_EXTS, zip_ext=ZIP_EXTS,
                 stdin_allowed=False):
    errors = []
    if args.input_paths is None:
        args.input_paths = [STDIO_PATH]
    try:
        input_files_indexed = args.index_paths is not None
    except AttributeError:
        input_files_indexed = False
    input_formats = set()
    for input_path in args.input_paths:
        error = check_file_arg(input_path, ext=ext, zip_ext=zip_ext,
                               stdio_allowed=stdin_allowed)
        if error:
            errors += error
        else:
            try:
                input_format = check_seq_file_format(
                    input_path, args.input_format,
                    default_format=SEQFILE_FORMAT,
                    file_indexed=input_files_indexed)
                input_formats |= {input_format}
            except ValueError as exc:
                error = ['Input file ', input_path, ': ', str(exc)]
                errors.append(''.join(error))
    if (stdin_allowed and len(args.input_paths) > 1 and
            STDIO_PATH in args.input_paths):
        errors.append('Standard input ({}) cannot be part of a list of '
                      'multiple input files.'.format(STDIO_PATH))
    if len(input_formats) > 1:
        errors.append('Input files must all have the same format.')
    elif len(input_formats) == 1:
        args.input_format = list(input_formats)[0]

    if not input_files_indexed:
        return errors
    if len(args.index_paths) not in [1, len(args.input_paths)]:
        errors.append('There must be one index file for each input file '
                      'or one for all.')
        return errors
    for index_path in args.index_paths:
        errors += check_file_arg(index_path, prefix='-i/--index_paths')
    return errors


def check_outputs(args, ext=SEQFILE_EXTS, zip_ext=ZIP_EXTS,
                  stdout_allowed=False):
    errors = []
    if isinstance(args.output_path, list):
        if args.output_path is None:
            if not hasattr(args, 'no_output') or not args.no_output:
                args.output_path = STDIO_PATH
        output_paths = args.output_path
    else:
        if args.output_path is None:
            if not hasattr(args, 'no_output') or not args.no_output:
                args.output_path = [STDIO_PATH]
        output_paths = [args.output_path]

    try:
        default_format = args.input_format
    except AttributeError:
        default_format = SEQFILE_FORMAT

    output_formats = set()
    for output_path in output_paths:
        error = check_file_arg(
            output_path, mode='w', ext=ext, zip_ext=zip_ext,
            stdio_allowed=stdout_allowed, prefix='-o/--output_path')
        if error:
            errors += error
        else:
            try:
                output_format = check_seq_file_format(
                    output_path, args.output_format,
                    default_format=default_format)
                output_formats |= {output_format}
            except ValueError as exc:
                error = ['-o/--output_path: ', output_path, ': ', str(exc)]
                errors.append(''.join(error))

    if len(output_formats) > 1:
        errors.append(prefix + ': Output files must all have the same format.')
    elif len(output_formats) == 1:
        args.output_format = list(output_formats)[0]

    if stdout_allowed and len(output_paths) > 1 and STDIO_PATH in output_paths:
        errors.append('Standard output ({}) cannot be part of a list of '
                      'multiple output files.'.format(STDIO_PATH))

    try:
        nb_inputs = len(args.input_paths)
    except AttributeError:
        return errors

    try:
        nb_indices = len(args.index_paths)
    except (AttributeError, TypeError):
        nb_indices = 0

    if nb_inputs > 1 and nb_indices == 1 and len(output_paths) != 1:
        errors.append('If a common index file is specified for all input '
                      'files, only one output file must be specified.')
    elif len(output_paths) not in [1, nb_inputs]:
        errors.append('There must be one output file for each input file '
                      'or one for all.')
    return errors


def check_cpu_count(args):
    return check_num_arg(args.cpus, mini=1, maxi=multiprocessing.cpu_count(),
                         prefix='-c/--cpus')
