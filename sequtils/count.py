#!/usr/bin/env python3

import argparse
import os
import sys

import Bio.SeqIO

from sequtils.common import check_sequtils_inputs, INPUT_FORMAT
from vjoly import (check_file_arg, export_table, open_file, print_table,
                   print_stderr, SEQFILE_FORMAT, FASTA_FORMAT, GENBANK_FORMAT,
                   STDIO_PATH)

COUNT_SEQS = True
COUNT_BASES = False
IGNORE_CASE = False
SEQ_COUNTS_PATH = None
BASE_COUNTS_PATH = None


def add_arguments(parser):
    parser.set_defaults(prog=parser.prog, check=check_arguments, process=main)

    seqio_group = parser.add_argument_group(title='Input options')
    seqio_group.add_argument(
        'input_paths', metavar='FILE', nargs='+',
        help='Sequence FASTA file(s)')
    input_format = seqio_group.add_mutually_exclusive_group()
    input_format.add_argument(
        '-f', '--fasta_input', action='store_const', const=FASTA_FORMAT,
        dest=INPUT_FORMAT, help='Input is in FASTA format.')
    input_format.add_argument(
        '-g', '--genbank_input', action='store_const', const=GENBANK_FORMAT,
        dest=INPUT_FORMAT, help='Input is in Genbank format.')

    settings_group = parser.add_argument_group(title='Counting options')
    settings_group.add_argument(
        '-s', '--count_seqs', action='store_true', default=COUNT_SEQS,
        help='Count sequences [{}]'.format(COUNT_SEQS))
    settings_group.add_argument(
        '-b', '--count_bases', action='store_true', default=COUNT_BASES,
        help='Count bases/aminoacids [{}]'.format(COUNT_BASES))
    settings_group.add_argument(
        '-i', '--ignore_case', action='store_true',
        help='Ignore case when counting bases [{}]'.format(IGNORE_CASE))

    output_group = parser.add_argument_group(title='Output options')
    output_group.add_argument(
        '-S', '--seq_counts_path', metavar='FILE', default=None,
        help='Write seq counts to this file [{}].'.format(SEQ_COUNTS_PATH))
    output_group.add_argument(
        '-B', '--base_counts_path', metavar='FILE', default=None,
        help='Write base counts to this file [{}].'.format(BASE_COUNTS_PATH))

    other_group = parser.add_argument_group(title="Other options")
    other_group.add_argument(
        '-h', '--help', action='help', help='Print this help page and exit.')


def check_arguments(args):
    errors = []
    errors += check_sequtils_inputs(args, stdin_allowed=True)
    errors += check_file_arg(
        args.seq_counts_path, mode='w', none_allowed=True,
        prefix='-S/--seq_counts_path')
    errors += check_file_arg(
        args.base_counts_path, mode='w', none_allowed=True,
        prefix='-C/--base_counts_path')
    if errors:
        for error in errors:
            print_stderr(error, prefix='ERROR')
        sys.exit(1)


def count_seqs(input_paths, input_format=SEQFILE_FORMAT, count_seqs=COUNT_SEQS,
               count_bases=COUNT_BASES, ignore_case=IGNORE_CASE):
    seq_counts, base_counts = {}, {}
    for input_path in input_paths:
        seq_count = 0
        with open_file(input_path) as input_file:
            records = Bio.SeqIO.parse(input_file, input_format)
            for record in records:
                if count_seqs:
                    seq_count += 1
                if count_bases:
                    seq = record.seq
                    if ignore_case:
                        seq = seq.upper()
                    for base in seq:
                        try:
                            base_counts[base][input_path] += 1
                        except KeyError:
                            try:
                                base_counts[base][input_path] = 1
                            except KeyError:
                                base_counts[base] = {input_path: 1}
        seq_counts[input_path] = seq_count
    return seq_counts, base_counts


def make_seq_counts_table(seq_counts, input_paths):
    header = ['File', 'Seq count']
    table = []
    for input_path in input_paths:
        if input_path != STDIO_PATH:
            file_name = os.path.basename(input_path)
        else:
            file_name = '<stdin>'
        row = [file_name, seq_counts[input_path]]
        table.append(row)
    return table, header


def make_base_counts_table(base_counts, input_paths):
    bases = sorted(list(base_counts.keys()))
    header = ['File'] + bases + ['Total']
    table = []
    for input_path in input_paths:
        counts = []
        for base in bases:
            try:
                count = base_counts[base][input_path]
            except:
                count = 0
            counts.append(count)
        total = sum(counts)
        counts.append(total)
        props = ['{: >5.1f}%'.format(100.0*count/total) for count in counts]
        if input_path != STDIO_PATH:
            file_name = os.path.basename(input_path)
        else:
            file_name = '<stdin>'
        table.append([file_name] + counts)
        table.append([''] + props)
    return table, header


def main(args):
    seq_counts, base_counts = count_seqs(
        args.input_paths,
        input_format=args.input_format,
        count_seqs=args.count_seqs,
        count_bases=args.count_bases,
        ignore_case=args.ignore_case)

    if args.count_seqs:
        table, header = make_seq_counts_table(seq_counts, args.input_paths)
        print_table(table, header=header)
        if args.seq_counts_path is not None:
            export_table(table, args.seq_counts_path, header=header)

    if args.count_bases:
        table, header = make_base_counts_table(base_counts, args.input_paths)
        print_table(table, header=header)
        if args.base_counts_path is not None:
            export_table(table, args.base_counts_path, header=header)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Count sequences', add_help=False)
    add_arguments(parser)
    args = parser.parse_args()
    args.check(args)
    args.process(args)
