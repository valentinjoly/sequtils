#!/usr/bin/env python3

import argparse
import sys

from vjoly import print_stderr


def add_arguments(parser):
    parser.set_defaults(prog=parser.prog, check=check_arguments, process=main)
    # Insert arguments here
    parser.add_argument(
        '-h', '--help', action='help', help='Print this help page and exit.')


def check_arguments(args):
    errors = []
    # Insert argument verifications here
    if errors:
        for error in errors:
            print_stderr(error, prefix='ERROR')
        sys.exit(1)


def main(args):
    pass


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='', add_help=False)
    add_arguments(parser)
    args = parser.parse_args()
    args.check(args)
    args.process(args)


