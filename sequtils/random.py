#!/usr/bin/env python3

import argparse
import multiprocessing
import random
import sys

import Bio.Seq
import Bio.SeqIO
import Bio.SeqRecord

from sequtils.common import (
    add_seqio_arguments, add_other_arguments, check_outputs, check_cpu_count,
    INPUT_FORMAT, OUTPUT_FORMAT)

from sequtils.parallel import make_tmp_output_paths, merge_outputs

from vjoly import (
    check_file_arg, check_num_arg, open_file, print_stderr, write_records,
    FASTA_FORMAT, GENBANK_FORMAT, SEQFILE_FORMAT, STDIO_PATH,
    DNA_ALPHA, PROT_ALPHA, RNA_ALPHA)

OUTPUT_PATH = STDIO_PATH
MIN_LENGTH = 100
MAX_LENGTH = None
NB_SEQUENCES = 1
PREFIX = 'seq'
DESCRIPTION = 'random sequence'

DNA = 'dna'
RNA = 'rna'
CDS = 'cds'
AA = 'aa'
PROT = 'prot'
PROT_STOP = 'prot_stop'
SEQ_TYPE = DNA

DNA_NTS = list(DNA_ALPHA.letters)
RNA_NTS = list(RNA_ALPHA.letters)
AMINOACIDS = list(PROT_ALPHA.letters)

CODONS = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA',
          'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATT', 'CAA', 'CAC', 'CAG',
          'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT',
          'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA',
          'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC',
          'GTG', 'GTT', 'TAC', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGC',
          'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']
START_CODONS = ['ATG']
STOP_CODONS = ['TAA', 'TAG', 'TGA']
START_AA = ['M']
STOP_AA = ['*']


def add_arguments(parser):
    parser.set_defaults(prog=parser.prog, check=check_arguments, process=main)

    io_group = parser.add_argument_group(title='Input/output options')
    add_seqio_arguments(io_group, add_input=False, add_index=False,
                        multiple_outputs=False)

    settings_group = parser.add_argument_group(
        title='Sequence generation options')
    settings_group.add_argument(
        '-n', '--nb_sequences', metavar='INT', type=int, default=NB_SEQUENCES,
        help='Number of sequences to generate [{}]'.format(NB_SEQUENCES))
    settings_group.add_argument(
        '-l', '--min_length', metavar='INT', type=int, default=MIN_LENGTH,
        help='Minimum sequence length [{}]'.format(MIN_LENGTH))
    settings_group.add_argument(
        '-L', '--max_length', metavar='INT', type=int, default=MAX_LENGTH,
        help='Maximum sequence length [{}]'.format(MIN_LENGTH))
    settings_group.add_argument(
        '-p', '--prefix', metavar='STR', default=PREFIX,
        help='Prefix for sequence IDs [{}]'.format(PREFIX))

    seq_type_group = parser.add_argument_group(
        title='Sequence type [default: {}]'.format(SEQ_TYPE))
    seq_types = seq_type_group.add_mutually_exclusive_group()
    seq_types.add_argument(
        '--dna', action='store_const', const=DNA,
        dest='seq_type', help='Generate random DNA sequences')
    seq_types.add_argument(
        '--rna', action='store_const', const=RNA,
        dest='seq_type', help='Generate random RNA sequences')
    seq_types.add_argument(
        '--cds', action='store_const', const=CDS,
        dest='seq_type', help='Generate random protein-coding DNA sequences')
    seq_types.add_argument(
        '--aa', action='store_const', const=AA,
        dest='seq_type', help='Generate random aminoacid sequences')
    seq_types.add_argument(
        '--prot', action='store_const', const=PROT,
        dest='seq_type', help='Generate random protein sequences (M...)')
    seq_types.add_argument(
        '--prot_stop', action='store_const', const=PROT_STOP,
        dest='seq_type', help='Generate random protein sequences (M...*)')

    other_group = parser.add_argument_group(title="Other options")
    add_other_arguments(other_group)


def check_arguments(args):
    errors = []
    errors += check_outputs(args, stdout_allowed=True)
    errors += check_cpu_count(args)
    errors += check_num_arg(
        args.min_length, number_type=int, mini=1,
        prefix='-l/--min_length')
    if args.max_length is None:
        args.max_length = args.min_length
    else:
        errors += check_num_arg(
            args.max_length, number_type=int, mini=args.min_length,
            prefix='-L/--max_length')
    errors += check_num_arg(
        args.nb_sequences, number_type=int, mini=1,
        prefix='-n/--nb_sequences')
    if args.seq_type is None:
        args.seq_type = SEQ_TYPE
    if errors:
        for error in errors:
            print_stderr(error, prefix='ERROR')
        sys.exit(1)


def define_length(min_length=MIN_LENGTH, max_length=MAX_LENGTH):
    return random.choice(range(min_length, max_length + 1))


def make_seqid(prefix, index, index_width):
    return prefix + str(index).zfill(index_width)


def generate_random_seq(length=MIN_LENGTH, characters=DNA_NTS,
                        start=None, end=None):
    prefix, suffix = '', ''
    if start is not None:
        prefix = random.choice(start)
    if end is not None:
        suffix = random.choice(end)
    return prefix + ''.join(random.choices(characters, k=length)) + suffix


def partition_counts(total, nb_parts=1):
    if total <= nb_parts:
        return [1] * total
    parts = [total // nb_parts] * nb_parts
    for i in range(total % nb_parts):
        parts[i] += 1
    return parts


def generate_random_sequences_part(
        output_path, output_format, seq_type, min_length, max_length,
        nb_sequences, prefix, start_index, index_width):
    with open_file(output_path, 'w') as output_file:
        for i in range(nb_sequences):
            seqid = make_seqid(prefix, i + start_index, index_width)
            length = define_length(min_length, max_length)
            if seq_type == 'dna':
                seq = generate_random_seq(length, DNA_NTS)
                alphabet = DNA_ALPHA()
            elif seq_type == 'rna':
                seq = generate_random_seq(length, RNA_NTS)
                alphabet = RNA_ALPHA()
            elif seq_type == 'aa':
                seq = generate_random_seq(length, AMINOACIDS)
                alphabet = PROT_ALPHA()
            elif seq_type == 'prot':
                seq = generate_random_seq(length-1, AMINOACIDS, start=START_AA)
                alphabet = PROT_ALPHA()
            elif seq_type == 'prot_stop':
                seq = generate_random_seq(length-1, AMINOACIDS,
                                          start=START_AA, end=STOP_AA)
                alphabet = PROT_ALPHA()
            elif seq_type == 'cds':
                seq = generate_random_seq(length//3-2, CODONS,
                                          start=START_CODONS, end=STOP_CODONS)
                alphabet = DNA_ALPHA()

            record = Bio.SeqRecord.SeqRecord(
                id=seqid,
                seq=Bio.Seq.Seq(seq, alphabet),
                description=DESCRIPTION)
            write_records(record, output_file, output_format)


def generate_random_sequences(
        output_path, output_format, seq_type, min_length,
        max_length, nb_sequences, prefix, nb_cpus):
    nb_sequences_parts = partition_counts(nb_sequences, nb_cpus)
    start_index_parts = [sum(nb_sequences_parts[:i]) + 1
                         for i in range(len(nb_sequences_parts) + 1)]
    nb_parts = len(nb_sequences_parts)
    if nb_parts > 1:
        tmp_output_paths = make_tmp_output_paths(output_path, nb_parts)
    else:
        tmp_output_paths = [output_path]
    parts = zip(tmp_output_paths, nb_sequences_parts, start_index_parts)
    index_width = len(str(nb_sequences))
    processes = []
    for tmp_output_path, nb_sequences_part, start_index_part in parts:
        processes.append(multiprocessing.Process(
            target=generate_random_sequences_part,
            args=(tmp_output_path, output_format, seq_type, min_length,
                  max_length, nb_sequences_part, prefix, start_index_part,
                  index_width)))
    for p in processes:
        p.start()
    for p in processes:
        p.join()

    if nb_parts > 1:
        merge_outputs(tmp_output_paths, output_path)


def main(args):
    generate_random_sequences(
        args.output_path,
        args.output_format,
        args.seq_type,
        args.min_length,
        args.max_length,
        args.nb_sequences,
        args.prefix,
        args.cpus)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Generate random sequences', add_help=False)
    add_arguments(parser)
    args = parser.parse_args()
    args.check(args)
    args.process(args)
