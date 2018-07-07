#!/usr/bin/env python3

import argparse
import re
import sys

from sequtils.common import (
    add_seqio_arguments, add_other_arguments, check_inputs, check_outputs,
    check_cpu_count)

from sequtils.parallel import (parallelize_indexed, parallelize_unindexed)

from vjoly import (
    check_file_arg, open_file, print_stderr, write_records,)

SEQIDS_PATH = None
DESCRIPTIONS_PATH = None
NAMES_PATH = None
MOL_TYPES_PATH = None
TAXIDS_PATH = None
EQUAL = 'equal'
REGEX = 'regex'
NOT_EQUAL = 'not_equal'
CONTAIN = 'contain'
BEGIN = 'begin'
END = 'end'
GEQ = 'geq'
LEQ = 'leq'
COMPARISON = EQUAL
INVERSE = False
SEARCH_IN_DESC = False
VALUES_IN_FILE = False
START_WITH = False
END_WITH = False

INPUT_PATHS = None
OUTPUT_PATH = None


def add_arguments(parser):
    parser.set_defaults(prog=parser.prog, check=check_arguments, process=main)

    io_group = parser.add_argument_group(title='Input/output options')
    add_seqio_arguments(io_group)

    criteria_group = parser.add_argument_group(
        title='Values to look for (repeat the flag to specify '
              'multiple values)')
    criteria_group.add_argument(
        '-s', '--seqid', metavar='STR', action='append', dest='seqids',
        help='Sequence identifier.')
    criteria_group.add_argument(
        '-d', '--desc', metavar='STR', action='append', dest='descriptions',
        help='Sequence description.')
    criteria_group.add_argument(
        '-n', '--name', metavar='STR', action='append', dest='names',
        help='Sequence name (for GenBank files only).')
    criteria_group.add_argument(
        '-m', '--mol_type', metavar='STR', action='append', dest='mol_types',
        help='Molecule type (for GenBank files only).')
    criteria_group.add_argument(
        '-t', '--taxid', metavar='INT', action='append',
        type=int, dest='taxids',
        help='Taxonomic ID (for GenBank files only).')

    criteria_path_group = parser.add_argument_group(
        'Path to files containing values to look for (one value per line)')
    criteria_path_group.add_argument(
        '-S', '--seqids_path', metavar='FILE', default=SEQIDS_PATH,
        help='List of sequence identifiers '
             '[{}].'.format(SEQIDS_PATH))
    criteria_path_group.add_argument(
        '-D', '--descriptions_path', metavar='FILE', default=DESCRIPTIONS_PATH,
        help='List of sequence descriptions '
             '[{}].'.format(DESCRIPTIONS_PATH))
    criteria_path_group.add_argument(
        '-N', '--names_path', metavar='FILE', default=NAMES_PATH,
        help='List of sequence names (for GenBank files only) '
             '[{}].'.format(NAMES_PATH))
    criteria_path_group.add_argument(
        '-M', '--mol_types_path', metavar='FILE', default=MOL_TYPES_PATH,
        help='List of molecule types (for GenBank files only) '
             '[{}].'.format(MOL_TYPES_PATH))
    criteria_path_group.add_argument(
        '-T', '--taxids_path', metavar='FILE', default=TAXIDS_PATH,
        help='List of taxonomic IDs (for GenBank files only) '
             '[{}].'.format(TAXIDS_PATH))

    comparison_group = parser.add_argument_group(
        'Comparison to be performed [default: {}]'.format(COMPARISON))
    comparisons = comparison_group.add_mutually_exclusive_group()
    comparisons.add_argument(
        '--equal', action='store_const', const=EQUAL,
        dest='comparison', help='Equals...')
    comparisons.add_argument(
        '--regex', action='store_const', const=REGEX,
        dest='comparison', help='Matches the regex...')
    comparisons.add_argument(
        '--contain', action='store_const', const=CONTAIN,
        dest='comparison', help='Contains...')
    comparisons.add_argument(
        '--begin', action='store_const', const=BEGIN,
        dest='comparison', help='Begins with...')
    comparisons.add_argument(
        '--end', action='store_const', const=END,
        dest='comparison', help='Ends with...')
    comparisons.add_argument(
        '--leq', action='store_const', const=LEQ,
        dest='comparison', help='Is lesser than or equal to...')
    comparisons.add_argument(
        '--geq', action='store_const', const=GEQ,
        dest='comparison', help='Is greater than or equal to...')

    result_group = parser.add_argument_group(
        'Expected result to retain a sequence')
    result_group.add_argument(
        '--inverse', action='store_true', default=INVERSE,
        help='Add this flag to retain only sequences NOT matching '
             'search criteria. [{}]'.format(INVERSE))

    other_group = parser.add_argument_group('Other options')
    add_other_arguments(other_group)


def check_arguments(args):
    errors = []
    errors += check_inputs(args, stdin_allowed=True)
    errors += check_outputs(args, stdout_allowed=True)
    errors += check_cpu_count(args)
    errors += check_file_arg(args.seqids_path, none_allowed=True,
                             prefix='-S/--seqids_path')
    errors += check_file_arg(args.descriptions_path, none_allowed=True,
                             prefix='-D/--descriptions_path')
    errors += check_file_arg(args.names_path, none_allowed=True,
                             prefix='-N/--names_path')
    errors += check_file_arg(args.mol_types_path, none_allowed=True,
                             prefix='-M/--mol_types_path')
    errors += check_file_arg(args.taxids_path, none_allowed=True,
                             prefix='-T/--taxids_path')
    if args.input_format != 'genbank':
        values = [args.names, args.mol_types, args.taxids]
        paths = [args.names_path, args.mol_types_path, args.taxids_path]
        prefixes = ['-n/-N', '-m/-M', '-t/-T']
        for value, path, prefix in zip(values, paths, prefixes):
            if value or path is not None:
                errors += ('{}: These options require input file(s) '
                           'in GenBank format.'.format(prefix))
    if args.comparison is None:
        args.comparison = COMPARISON
    if errors:
        for error in errors:
            print_stderr(error, prefix='ERROR')
        sys.exit(1)


def import_values(values_path, integers=False):
    values = []
    with open_file(values_path) as input_file:
        for line in input_file:
            value = line.strip()
            if not value:
                continue
            if integers:
                value = int(value)
            values.append(value)
    return values


def get_values(values_list, values_path, integers=False):
    values = []
    if values_list is not None:
        values += values_list
    if values_path is not None:
        values += import_values(values_path, integers)
    values.sort()
    return values


def is_equal(value, reference):
    return value == reference


def matches_regex(value, reference):
    return re.search(reference, value)


def contains(value, reference):
    return str(reference) in str(value)


def begins_with(value, reference):
    return str(value).startswith(str(reference))


def ends_with(value, reference):
    return str(value).endswith(str(reference))


def is_geq(value, reference):
    return value >= reference


def is_leq(value, reference):
    return value <= reference


def set_comparison_func(comparison=COMPARISON):
    if comparison == 'equal':
        return is_equal
    if comparison == 'regex':
        return matches_regex
    if comparison == 'contain':
        return contains
    if comparison == 'begin':
        return begins_with
    if comparison == 'end':
        return ends_with
    if comparison == 'geq':
        return is_geq
    if comparison == 'leq':
        return is_leq


def get_record_desc(record):
    return record.description[len(record.id):].strip()


def get_record_source(record):
    if not record.features:
        return None
    for i in range(len(record.features)):
        if record.features[i].type == 'source':
            return record.features[i]
    return None


def get_record_mol_type(record_source):
    try:
        mol_type = record_source.qualifiers['mol_type'][0]
    except (KeyError, IndexError):
        return None
    else:
        return mol_type


def get_record_taxid(record_srouce):
    try:
        db_xref = record_srouce.qualifiers['db_xref']
    except KeyError:
        return None
    for taxdb_ref in db_xref:
        if taxdb_ref.startswith('taxon:'):
            taxid = int(taxdb_ref[taxdb_ref.index(':')+1:])
            return taxid
    return None


def good_value(record_value, reference_values, comparison_func):
    for reference_value in reference_values:
        if comparison_func(record_value, reference_value):
            return True
    return False


def good_record(record, seqids, descriptions, names,
                mol_types, taxids, comparison_func):
    if seqids:
        rec_seqid = record.id
        if not good_value(rec_seqid, seqids, comparison_func):
            return False
    if descriptions:
        rec_description = get_record_desc(record)
        if not good_value(rec_description, descriptions, comparison_func):
            return False
    if names:
        rec_name = record.name
        if not good_value(rec_name, names, comparison_func):
            return False
    if mol_types or taxids:
        rec_source = get_record_source(record)
        if rec_source is None:
            return False
    if mol_types:
        rec_mol_type = get_record_mol_type(rec_source)
        if rec_mol_type is None:
            return False
        if not good_value(rec_mol_type, mol_types, comparison_func):
            return False
    if taxids:
        rec_taxid = get_record_taxid(rec_source)
        if rec_taxid is None:
            return False
        if not good_value(rec_taxid, taxids, comparison_func):
            return False
    return True


def extract_records_complex(records, output_file, output_format, seqids=None,
                            descriptions=None, names=None, mol_types=None,
                            taxids=None, comparison_func=None, inverse=None):
    for record in records:
        keep_record = good_record(
            record, seqids, descriptions, names,
            mol_types, taxids, comparison_func)
        if inverse:
            keep_record = not keep_record
        if not keep_record:
            continue
        write_records(record, output_file, output_format)


def extract_records_by_seqid(records, output_file, output_format,
                             seqids=None):
    seqids_found = {seqid: False for seqid in seqids}
    for seqid in seqids:
        try:
            record = records[seqid]
        except KeyError:
            pass
        else:
            write_records(record, output_file, output_format)
            seqids_found[seqid] = True
    return seqids_found


def check_if_simple_case(seqids, descriptions, names, mol_types, taxids,
                         inverse, comparison):
    if comparison != EQUAL:
        return False
    if inverse:
        return False
    if not seqids:
        return False
    if descriptions or names or mol_types or taxids:
        return False
    return True


def print_missing_seqids(seqids_found):
    not_found = []
    for seqid in seqids_found:
        if not seqids_found[seqid]:
            not_found.append(seqid)
    if not_found:
        print_stderr('{:d} identifier(s) not found: {}'.format(
                        len(not_found), ', '.join(not_found)))


def get_search_settings(args):
    seqids = get_values(args.seqids, args.seqids_path)
    descriptions = get_values(args.descriptions, args.descriptions_path)
    names = get_values(args.names, args.names_path)
    mol_types = get_values(args.mol_types, args.mol_types_path)
    taxids = get_values(args.taxids, args.taxids_path)
    if not (seqids or descriptions or names or mol_types or taxids):
        print_stderr('ERROR: No search criterion specified.')
        sys.exit(0)

    simple_case = check_if_simple_case(
        seqids, descriptions, names, mol_types, taxids,
        args.inverse, args.comparison)

    comparison_func = set_comparison_func(args.comparison)
    return (seqids, descriptions, names, mol_types, taxids,
            comparison_func, simple_case)


def main(args):
    (seqids, descriptions, names, mol_types,
        taxids,  comparison_func, simple_case) = get_search_settings(args)

    if simple_case:
        seqids_found = parallelize_indexed(
            extract_records_by_seqid,
            args.cpus,
            args.input_paths,
            args.index_paths,
            args.output_paths,
            args.input_format,
            args.output_format,
            merge_results_function=any,
            seqids=seqids)
        print_missing_seqids(seqids_found)
    else:
        parallelize_unindexed(
            extract_records_complex,
            args.cpus,
            args.input_paths,
            args.output_paths,
            args.input_format,
            args.output_format,
            seqids=seqids,
            descriptions=descriptions,
            names=names,
            mol_types=mol_types,
            taxids=taxids,
            comparison_func=comparison_func,
            inverse=args.inverse)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Extract sequences by their ID or information',
        add_help=False)
    add_arguments(parser)
    args = parser.parse_args()
    args.check(args)
    args.process(args)
