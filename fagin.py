#!/usr/bin/env python3

import argparse
import sys
import os
import itertools
import collections
import math

import lib.util           as util
import lib.intervals      as intervals
import lib.genome         as genome
import lib.synteny        as synteny
import lib.exonerate      as exonerate
import lib.hit_merger     as hit_merger
import lib.syn_merger     as syn_merger
import lib.result_manager as result_manager

__version__ = '0.0.1'

def parse(argv=None):
    parser = argparse.ArgumentParser(
        description='Discover and categorize orphan genes',
        usage='fagin [options]'
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {}'.format(__version__)
    )

    parser.add_argument(
        '-o', '--output_dir',
        help='output directory (must be empty or nonexistant)',
        default='output'
    )

    parser.add_argument(
        '--gff',
        type=argparse.FileType('r'),
        help='gff formated gene models for the query species'
    )

    parser.add_argument(
        '--synteny',
        type=argparse.FileType('r'),
        help='output tabular output from SatsumaSynteny (query versus target)'
    )

    parser.add_argument(
        '--exonerate',
        type=argparse.FileType('r'),
        help='the parsed output of exonerate'
    )

    parser.add_argument(
        '--nstring',
        type=argparse.FileType('r'),
        help='tab-delimited file representing chr, start, and length of N repeats in target genome'
    )

    parser.add_argument(
        '-w', '--context_width',
        help='the number of upstream and downstream synteny blocks to include in the analysis',
        metavar='N',
        type=int,
        default=10
    )

    args = parser.parse_args(argv)
    return(args)


if __name__ == '__main__':
    args = parse()
    try:
        os.mkdir(args.output_dir)
    except FileExistsError:
        if(os.listdir(args.output_dir)):
            err('Output directory must be empty')
    except PermissionError:
        err("You don't have permission to make directory '%s'" % args.output_dir)

    res = result_manager.ResultManager(
        gen        = genome.Genome(args.gff),
        syn        = synteny.Synteny(args.synteny),
        exo        = exonerate.Exonerate(args.exonerate),
        hit_merger = hit_merger.HitMerger(),
        syn_merger = syn_merger.SynMerger()
    )
    res.write()
