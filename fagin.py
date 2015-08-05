#!/usr/bin/env python3

import argparse
import sys
import os
import itertools
import collections
import math

import src.util      as util
import src.intervals as intervals
import src.genome    as genome
import src.synteny   as synteny

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

class Result:
    def __init__(self, gene, syn, width=1):
        self.anchor = syn.anchor_query(gene)
        self.gene = gene
        self.width = width
        if self.anchor:
            self.links = self._get_links()

            self.context = self._get_context()

            # does the gene overlap a syntenic block?
            self.is_present = bool(self.links)

            # are the upstream and downstream blocks in order on the same
            # chromosome?
            self.is_simple = self._get_is_simple()
        else:
            self.context = None
            self.is_present = False
            self.is_simple = False

    def _get_links(self):
        m = self.anchor
        links = []
        while m and intervals.overlaps(m, self.gene):
            links.append(m)
            m = m.last
        m = self.anchor.next
        while m and intervals.overlaps(m, self.gene):
            links.append(m)
            m = m.next
        links = sorted(links, key=lambda x: (x.start, x.stop))
        return(links)

    def _get_context(self):
        # Get the syntenic blocks flanking (but not overlapping) the gene
        if not self.links:
            if self.gene.stop < self.anchor.start:
                lower, upper = (self.anchor.last, self.anchor)
            else:
                lower, upper = (self.anchor, self.anchor.next)
        else:
            lower, upper = (self.links[0].last, self.links[-1].next)

        # get the context of gene
        lower_context = intervals.get_preceding(lower, self.width)
        upper_context = intervals.get_following(upper, self.width)
        everything = [x for x in itertools.chain(lower_context, self.links, upper_context) if x]

        return(everything)

    def _get_is_simple(self):
        # are all the intervals on the same contig?
        all_on_same_contig = intervals.allequal((x.contig for x in self.context))

        # make an interval describing the start and stop of the query context
        minstart = min(x.start for x in self.context)
        maxstop  = max(x.stop  for x in self.context)
        query_bound = intervals.Interval(contig=self.anchor.contig, start=minstart, stop=maxstop)

        # do all the syntenic blocks within the target range map to regions within the query range?
        t = sorted(self.context, key=lambda x: (x.over.start))
        q = t[0]
        has_outer = False
        while True:
            if not q:
                break
            elif q.start > t[-1].stop:
                break
            elif not intervals.overlaps(q.over, query_bound):
                has_outer = True
                break
            else:
                q = q.next

        is_simple = all_on_same_contig and not has_outer
        return(is_simple)

    def __str__(self):
        out = '\t'.join((self.gene.name, str(self.is_present), str(self.is_simple)))
        return(out)


if __name__ == '__main__':
    args = parse()
    try:
        os.mkdir(args.output_dir)
    except FileExistsError:
        if(os.listdir(args.output_dir)):
            err('Output directory must be empty')
    except PermissionError:
        err("You don't have permission to make directory '%s'" % args.output_dir)

    gen = genome.Genome(args.gff)
    syn = synteny.Synteny(args.synteny)
    for g in gen.intervals():
        print(Result(gene=g, syn=syn))

