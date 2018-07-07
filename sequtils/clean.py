#!/usr/bin/env python3

import argparse
import re
import sys

import Bio.SeqIO

from sequtils.common import (
    add_seqio_arguments, add_other_arguments,
    check_inputs, check_outputs, check_cpu_count)

from sequtils.parallel import parallelize_indexed

from vjoly import (
    check_file_arg, check_num_arg, check_seq_file_format,
    open_file, print_stderr, write_records, FASTA_FORMAT,
    GENBANK_FORMAT, SEQFILE_FORMAT, STDIO_PATH)

NO_OUTPUT = False
REMOVED_PATH = None


def add_arguments(parser):
    parser.set_defaults(prog=parser.prog, check=check_arguments, process=main)

    io_group = parser.add_argument_group(title='Input/output options')
    add_seqio_arguments(io_group)
    io_group.add_argument(
        '-r', '--removed_path', metavar='FILE', default=REMOVED_PATH,
        help='Output list of removed ids. [{}]'.format(REMOVED_PATH))
    io_group.add_argument(
        '-n', '--no_output', action='store_true', default=NO_OUTPUT,
        help='Do not output clean sequences. [{}]'.format(NO_OUTPUT))

    other_group = parser.add_argument_group(title="Other options")
    add_other_arguments(other_group)


def check_arguments(args):
    errors = []
    errors += check_inputs(args)
    errors += check_outputs(args, stdout_allowed=True)
    errors += check_cpu_count(args)
    errors += check_file_arg(
        args.removed_path, mode='w', none_allowed=True, prefix='-r/--removed')
    if args.no_output:
        args.output_paths = None
    if errors:
        for error in errors:
            print_stderr(error, prefix='ERROR')
        sys.exit(1)


def clean_sequences(records, output_file, output_format):
    bad_records = []
    for seqid in records:
        try:
            record = records[seqid]
        except UnicodeDecodeError:
            bad_records.append(seqid)
        else:
            write_records(record, output_file, output_format)
    return bad_records


def export_removed_list(bad_records, removed_path=REMOVED_PATH):
    if not bad_records:
        print_stderr('No bad record found.', time=False)
    else:
        nb_bad = len(bad_records)
        bad = '\n'.join(sorted(bad_records))
        if nb_bad == 1:
            msg = '1 bad record found: {}'.format(bad)
        else:
            msg = '{:d} bad records found:\n{}'.format(nb_bad, bad)
        print_stderr(msg, time=False)
        if removed_path is not None:
            with open(removed_path, 'w') as removed_file:
                removed_file.write(bad)


def main(args):
    bad_records = parallelize_indexed(
        clean_sequences,
        args.cpus,
        args.input_paths,
        args.index_paths,
        args.output_paths,
        args.input_format,
        args.output_format)
    export_removed_list(
        bad_records,
        removed_path=args.removed_path)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Remove sequences with encoding errors',
        add_help=False)
    add_arguments(parser)
    args = parser.parse_args()
    args.check(args)
    args.process(args)
