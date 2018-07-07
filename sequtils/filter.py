#!/usr/bin/env python3

import argparse
import sys

import Bio.SeqIO

from sequtils.common import (check_sequtils_inputs, check_sequtils_outputs,
                             INPUT_FORMAT, OUTPUT_FORMAT, STDIO_PATH)
from vjoly import (check_file_arg, check_num_arg, open_file, print_stderr,
                   write_records, FASTA_FORMAT, GENBANK_FORMAT, SEQFILE_FORMAT)

INPUT_PATH = STDIO_PATH
INPUT_PATHS = [INPUT_PATH]
OUTPUT_PATH = STDIO_PATH
ALPHABET = None
MIN_LENGTH = 1
MAX_LENGTH = None
COUNTED_CHAR = None
MIN_COUNT = 0
MAX_COUNT = None
MIN_PROP = 0.0
MAX_PROP = 100.0


def add_arguments(parser):
    parser.set_defaults(prog=parser.prog, check=check_arguments, process=main)

    seqio_group = parser.add_argument_group(title='Input/output options')
    seqio_group.add_argument(
        'input_paths', metavar='FILE', nargs='+', default=INPUT_PATHS,
        help='Input sequence file(s). [{}]'.format(INPUT_PATH))
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

    settings_group = parser.add_argument_group(title='Filtering options')
    settings_group.add_argument(
        '-l', '--min_length', metavar='INT', type=int, default=MIN_LENGTH,
        help='Minimum sequence length [{}]'.format(MIN_LENGTH))
    settings_group.add_argument(
        '-L', '--max_length', metavar='INT', type=int, default=MAX_LENGTH,
        help='Maximum sequence length [{}]'.format(MAX_LENGTH))
    settings_group.add_argument(
        '-x', '--counted_char', metavar='STR', default=COUNTED_CHAR,
        help='Character to be counted in sequences [{}]'.format(COUNTED_CHAR))
    settings_group.add_argument(
        '-n', '--min_count', metavar='INT', type=int, default=MIN_COUNT,
        help='Minimum char count [{}]'.format(MIN_COUNT))
    settings_group.add_argument(
        '-N', '--max_count', metavar='INT', type=int, default=MAX_COUNT,
        help='Maximum char count [{}]'.format(MAX_COUNT))
    settings_group.add_argument(
        '-p', '--min_prop', metavar='FLOAT', type=float, default=MIN_PROP,
        help='Minimum char prop [{}]'.format(MIN_PROP))
    settings_group.add_argument(
        '-P', '--max_prop', metavar='FLOAT', type=float, default=MAX_PROP,
        help='Maximum char prop [{}]'.format(MAX_PROP))

    other_group = parser.add_argument_group(title="Other options")
    other_group.add_argument(
        '-h', '--help', action='help', help='Print this help page and exit.')


def check_arguments(args):
    errors = []
    errors += check_sequtils_inputs(args, stdin_allowed=True)
    errors += check_sequtils_outputs(args, stdout_allowed=True)
    errors += check_num_arg(
        args.min_length, number_type=int, mini=1,
        prefix='-l/--min_length')
    errors += check_num_arg(
        args.max_length, number_type=int, mini=args.min_length,
        none_allowed=True, prefix='-L/--max_length')
    errors += check_num_arg(
        args.min_count, number_type=int, mini=0,
        prefix='-n/--min_count')
    errors += check_num_arg(
        args.max_count, number_type=int, mini=args.min_count,
        none_allowed=True, prefix='-N/--max_count')
    errors += check_num_arg(
        args.min_prop, number_type=float, mini=0.0, maxi=100.0,
        prefix='-p/--min_prop')
    errors += check_num_arg(
        args.max_prop, number_type=float, mini=args.min_prop, maxi=100.0,
        prefix='-P/--max_prop')
    if args.counted_char is not None and len(args.counted_char) != 1:
        errors.append('-c/--counted_char: Must be one character.')

    if errors:
        for error in errors:
            print_stderr(error, prefix='ERROR')
        sys.exit(1)


def filter_seqs(input_paths=INPUT_PATHS, input_format=SEQFILE_FORMAT,
                output_path=OUTPUT_PATH, output_format=SEQFILE_FORMAT,
                min_length=MIN_LENGTH, max_length=MAX_LENGTH,
                counted_char=COUNTED_CHAR, min_count=MIN_COUNT,
                max_count=MAX_COUNT, min_prop=MIN_PROP, max_prop=MAX_PROP):
    with open_file(output_path, 'w') as output_file:
        for input_path in input_paths:           
            with open_file(input_path) as input_file:
                for record in Bio.SeqIO.parse(input_file, input_format):
                    length = len(record)
                    if length < min_length:
                        continue
                    if max_length is not None and length > max_length:
                        continue
                    if counted_char is not None:
                        count = record.seq.count(counted_char)
                        prop = 100.0 * count/length
                        if count < min_count:
                            continue
                        if max_count is not None and count > max_count:
                            continue
                        if prop < min_prop:
                            continue
                        if prop > max_prop:
                            continue
                    write_records(record, output_file, output_format)


def main(args):
    filter_seqs(
        input_paths=args.input_paths,
        input_format=args.input_format,
        output_path=args.output_path,
        output_format=args.output_format,
        min_length=args.min_length,
        max_length=args.max_length,
        counted_char=args.counted_char,
        min_count=args.min_count,
        max_count=args.max_count,
        min_prop=args.min_prop,
        max_prop=args.max_prop)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Filter sequences', add_help=False)
    add_arguments(parser)
    args = parser.parse_args()
    args.check(args)
    args.process(args)
