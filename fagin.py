#!/usr/bin/env python3

import argparse
import os

import lib.util           as util
import lib.genome         as genome
import lib.synteny        as synteny
import lib.exonerate      as exonerate
import lib.hit_merger     as hit_merger
import lib.syn_merger     as syn_merger
import lib.hit_analyzer   as hit_analyzer
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

    # === OUTPUTS ===

    parser.add_argument(
        '-o', '--output_dir',
        help='output directory (must be empty or nonexistant)',
        default='output'
    )

    parser.add_argument(
        '-q', '--quiet',
        help='suppress "contig with no syntenic block" warnings',
        action="store_true",
        default=False
    )

    # === INPUTS ===

    parser.add_argument(
        '-g', '--gen-file',
        help='gff formated gene models for the query species',
        type=argparse.FileType('r')
    )

    parser.add_argument(
        '-s', '--syn-file',
        help='output tabular output from SatsumaSynteny (query versus target)',
        type=argparse.FileType('r')
    )

    parser.add_argument(
        '-t', '--hit-file',
        help='the parsed output of Exonerate',
        type=argparse.FileType('r')
    )

    parser.add_argument(
        '-N', '--nstr-file',
        help='tab-delimited file representing chr, start, and length of N repeats in target genome',
        type=argparse.FileType('r')
    )

    # === PARAMETERS ===

    parser.add_argument(
        '-w', '--hit-flank-width',
        help='width of the upstream and downstream flanks surrounding a gene in which to search for syntenic blocks',
        type=int,
        default=25000
    )

    parser.add_argument(
        '-b', '--hit-min-neighbors',
        help='the minimum number of syntenic blocks near a gene which must be near a hit in order to keep the hit',
        type=int,
        default=3
    )

    parser.add_argument(
        '-r', '--hit-target-flank-ratio',
        help='the ratio between the query and target context widths',
        type=float,
        default=2
    )

    parser.add_argument(
        '-c', '--syn-context-width',
        help='the number of upstream and downstream synteny blocks to include in the analysis',
        metavar='N',
        type=int,
        default=10
    )

    args = parser.parse_args(argv)
    return(args)

def prepare_output_directory(args):
    try:
        os.mkdir(args.output_dir)
    except FileExistsError:
        if(os.listdir(args.output_dir)):
            err('Output directory must be empty')
    except PermissionError:
        err("You don't have permission to make directory '%s'" % args.output_dir)


if __name__ == '__main__':
    args = parse()

    syn_merger = syn_merger.SynMerger(
        width = args.syn_context_width
    )

    hit_merger = hit_merger.HitMerger(
        flank_width        = args.hit_flank_width,
        min_neighbors      = args.hit_min_neighbors,
        target_flank_ratio = args.hit_target_flank_ratio,
        quiet              = args.quiet
    )

    hit_analyzer = hit_analyzer.HitAnalyzer()

    res = result_manager.ResultManager(
        gen          = genome.Genome(args.gen_file),
        syn          = synteny.Synteny(args.syn_file),
        exo          = exonerate.Exonerate(args.hit_file),
        hit_merger   = hit_merger,
        syn_merger   = syn_merger,
        hit_analyzer = hit_analyzer
    )
    res.write()
