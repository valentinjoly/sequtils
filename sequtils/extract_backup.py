#!/usr/bin/env python3

import argparse
import re
import sys

import Bio.SeqIO

from sequtils.common import (check_sequtils_inputs, check_sequtils_outputs,
                             make_indexed_records_list, get_indexed_record,
                             INPUT_FORMAT, OUTPUT_FORMAT, STDIO_PATH)
from vjoly import (check_file_arg, check_num_arg, open_file, print_stderr,
                   write_records, FASTA_FORMAT, GENBANK_FORMAT, SEQFILE_FORMAT)

INPUT_PATH = STDIO_PATH
INPUT_PATHS = [INPUT_PATH]
INDEX_PATH = None
OUTPUT_PATH = STDIO_PATH
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


def add_arguments(parser):
    parser.set_defaults(prog=parser.prog, check=check_arguments, process=main)

    seqio_group = parser.add_argument_group(title='Input/output options')
    seqio_group.add_argument(
        'input_paths', metavar='FILE', nargs='+', default=INPUT_PATHS,
        help='Input sequence file(s). [{}]'.format(INPUT_PATH))
    seqio_group.add_argument(
        '-i', '--index_path', metavar='FILE', default=INDEX_PATH,
        help='Index file for input file(s). [{}]'.format(INDEX_PATH))
    seqio_group.add_argument(
        '-o', '--output_path', metavar='FILE', default=OUTPUT_PATH,
        help='Output sequence file. [{}]'.format(OUTPUT_PATH))
    input_format = seqio_group.add_mutually_exclusive_group()
    input_format.add_argument(
        '-f', '--fasta_input', action='store_const', const=FASTA_FORMAT,
        dest=INPUT_FORMAT, help='Input is in FASTA format.')
    input_format.add_argument(
        '-g', '--genbank_input', action='store_const', const=GENBANK_FORMAT,
        dest=INPUT_FORMAT, help='Input is in Genbank format.')
    output_format = seqio_group.add_mutually_exclusive_group()
    output_format.add_argument(
        '-F', '--fasta_output', action='store_const', const=FASTA_FORMAT,
        dest=OUTPUT_FORMAT, help='Write output in FASTA format.')
    output_format.add_argument(
        '-G', '--genbank_output', action='store_const', const=GENBANK_FORMAT,
        dest=OUTPUT_FORMAT, help='Write output in Genbank format.')

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
        '-x', '--equal', action='store_const', const=EQUAL,
        dest='comparison', help='Equals...')
    comparisons.add_argument(
        '-r', '--regex', action='store_const', const=REGEX,
        dest='comparison', help='Matches the regex...')
    comparisons.add_argument(
        '-c', '--contain', action='store_const', const=CONTAIN,
        dest='comparison', help='Contains...')
    comparisons.add_argument(
        '-b', '--begin', action='store_const', const=BEGIN,
        dest='comparison', help='Begins with...')
    comparisons.add_argument(
        '-e', '--end', action='store_const', const=END,
        dest='comparison', help='Ends with...')
    comparisons.add_argument(
        '-l', '--leq', action='store_const', const=LEQ,
        dest='comparison', help='Is lesser than or equal to...')
    comparisons.add_argument(
        '-u', '--geq', action='store_const', const=GEQ,
        dest='comparison', help='Is greater than or equal to...')

    result_group = parser.add_argument_group(
        'Expected result to retain a sequence')
    result_group.add_argument(
        '-z', '--inverse', action='store_true', default=INVERSE,
        help='Add this flag to retain only sequences NOT matching '
             'search criteria. [{}]'.format(INVERSE))

    other_group = parser.add_argument_group(title="Other options")
    other_group.add_argument(
        '-h', '--help', action='help', help='Print this help page and exit.')


def check_arguments(args):
    errors = []
    errors += check_sequtils_inputs(args, stdin_allowed=True)
    errors += check_sequtils_outputs(args, stdout_allowed=True)
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
    if errors:
        for error in errors:
            print_stderr(error, prefix='ERROR')
        sys.exit(1)


def import_values(values_path, integers=False):
    values = []
    with open_file(values_path) as input_file:
        for line in file:
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
        db_xref = source.qualifiers['db_xref']
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


def extract_records_complex(seqids, descriptions, names, mol_types, taxids,
                            comparison_func, inverse=INVERSE,
                            input_paths=INPUT_PATHS,
                            input_format=SEQFILE_FORMAT,
                            output_path=OUTPUT_PATH,
                            output_format=SEQFILE_FORMAT):
    with open_file(output_path, 'w') as output_file:
        output_alphabet = None
        for input_path in input_paths:
            with open_file(input_path) as input_file:
                records = Bio.SeqIO.parse(input_file, input_format)
                for record in records:
                    keep_record = good_record(
                        record, seqids, descriptions, names,
                        mol_types, taxids, comparison_func)
                    if inverse:
                        keep_record = not keep_record
                    if not keep_record:
                        continue
                    write_records(record, output_file, output_format)


def extract_records_by_seqid(seqids, input_paths=INPUT_PATHS,
                             input_format=SEQFILE_FORMAT,
                             index_path=INDEX_PATH, output_path=OUTPUT_PATH,
                             output_format=SEQFILE_FORMAT):
    records_list = make_indexed_records_list(
        input_paths, index_path, input_format)
    with open_file(output_path, 'w') as output_file:
        not_found = []
        for seqid in seqids:
            try:
                record = get_indexed_record(seqid, records_list)
            except KeyError:
                not_found.append(seqid)
            else:
                write_records(record, output_file, output_format)
    if not_found:
        print_stderr('{:d} identifier(s) not found: {}'.format(
                        len(not_found), ', '.join(not_found)))


def simple_case(args, seqids, descriptions, names, mol_types, taxids):
    if args.comparison != EQUAL:
        return False
    if args.inverse:
        return False
    if not seqids:
        return False
    if descriptions or names or mol_types or taxids:
        return False
    return True


def main(args):
    seqids = get_values(args.seqids, args.seqids_path)
    descriptions = get_values(args.descriptions, args.descriptions_path)
    names = get_values(args.names, args.names_path)
    mol_types = get_values(args.mol_types, args.mol_types_path)
    taxids = get_values(args.taxids, args.taxids_path)
    if args.comparison is None:
        args.comparison = COMPARISON
    if simple_case(args, seqids, descriptions, names, mol_types, taxids):
        extract_records_by_seqid(
            seqids,
            input_paths=args.input_paths,
            input_format=args.input_format,
            index_path=args.index_path,
            output_path=args.output_path,
            output_format=args.output_format)
    else:
        comparison_func = set_comparison_func(args.comparison)
        extract_records_complex(
            seqids,
            descriptions,
            names,
            mol_types,
            taxids,
            comparison_func,
            inverse=args.inverse,
            input_paths=args.input_paths,
            input_format=args.input_format,
            output_path=args.output_path,
            output_format=args.output_format)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Extract sequences by their ID or information',
        add_help=False)
    add_arguments(parser)
    args = parser.parse_args()
    args.check(args)
    args.process(args)
