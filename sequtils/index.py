#!/usr/bin/env python3

import argparse
import os
import sys

import Bio.SeqIO

from sequtils.common import (
    add_seqio_arguments, add_other_arguments, check_inputs)

from vjoly import (
    check_file_arg, print_stderr, BGZIP_EXT, SEQFILE_FORMAT)


def add_arguments(parser):
    parser.set_defaults(prog=parser.prog, check=check_arguments, process=main)
    io_group = parser.add_argument_group(title='Input/output options')
    add_seqio_arguments(io_group, add_index=False, add_output=False)
    io_group.add_argument(
        '-o', '--output_path', metavar='FILE', default=None,
        help='Output index file.')

    other_group = parser.add_argument_group('Other options')
    add_other_arguments(other_group, add_cpus=False)


def check_arguments(args):
    errors = []
    errors += check_inputs(args, zip_ext=BGZIP_EXT)
    if args.output_path is None:
        if len(args.input_paths) > 1:
            args.output_path = 'index.idx'
        else:
            input_basename = os.path.basename(args.input_paths[0])
            args.output_path = input_basename + '.idx'
    errors += check_file_arg(
        args.output_path, mode='w', prefix='-o/--output_path')

    if errors:
        for error in errors:
            print_stderr(error, prefix='ERROR')
        sys.exit(1)


def build_index(input_paths, output_path, input_format=SEQFILE_FORMAT):
    if os.path.exists(output_path):
        os.remove(output_path)
    Bio.SeqIO.index_db(output_path, input_paths, input_format)


def main(args):
    build_index(
        args.input_paths,
        output_path=args.output_path,
        input_format=args.input_format)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='', add_help=False)
    add_arguments(parser)
    args = parser.parse_args()
    args.check(args)
    args.process(args)
