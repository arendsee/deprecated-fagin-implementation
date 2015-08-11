#!/usr/bin/env python3

import argparse
import sys
import os
import itertools
import collections
import math

import lib.util      as util
import lib.intervals as intervals
import lib.genome    as genome
import lib.synteny   as synteny
import lib.exonerate as exonerate

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


class ResultSet:
    def __init__(self, gen, syn, exo, width=1):
        self.results = {g.name: Result(g, syn, width) for g in gen.intervals()}
        for hit in exo.generator():
            self.results[hit.name].add_exonerate_hit(hit, syn)

    def get(self, name):
        return self.results[name]

class Result:
    '''
    Synthesizes all data sets into final report
    '''
    def __init__(self, gene, syn, width=1):
        self.name = gene.name
        self.gene = gene

        # variables set in _syntenic_analysis
        self.context = None
        self.is_present = False
        self.is_simple = False
        self._syntenic_analysis(syn=syn, width=width)

        # variables modified when adding exonerate hits
        self.hits = []

    def add_exonerate_hit(self, hit, syn):
        assert(hit.name == self.name)

        anchor = syn.anchor_query(hit.target)

    def _syntenic_analysis(self, syn, width):
        anchor = syn.anchor_query(self.gene)
        if anchor:
            links = self._get_links(anchor=anchor)

            context = self._get_context(anchor=anchor, links=links, width=width)

            # does the gene overlap a syntenic block?
            is_present = bool(links)

            # are the upstream and downstream blocks in order on the same
            # chromosome?
            is_simple = self._get_is_simple(anchor=anchor, context=context)

        self.context    = context
        self.is_present = is_present
        self.is_simple  = is_simple

    def _get_links(self, anchor):
        '''
        find contiguous synteny blocks in the target that all overlap query gene
        '''
        m = anchor
        links = []
        while m and intervals.overlaps(m, self.gene):
            links.append(m)
            m = m.last
        m = anchor.next
        while m and intervals.overlaps(m, self.gene):
            links.append(m)
            m = m.next
        links = sorted(links, key=lambda x: (x.start, x.stop))
        return(links)

    def _get_context(self, anchor, links, width):
        # Get the syntenic blocks flanking (but not overlapping) the gene
        if not links:
            if self.gene.stop < anchor.start:
                lower, upper = (anchor.last, anchor)
            else:
                lower, upper = (anchor, anchor.next)
        else:
            lower, upper = (links[0].last, links[-1].next)

        # get the context of gene
        lower_context = intervals.get_preceding(lower, width)
        upper_context = intervals.get_following(upper, width)
        everything = [x for x in itertools.chain(lower_context, links, upper_context) if x]

        return(everything)

    def _get_is_simple(self, anchor, context):
        # are all the intervals on the same contig?
        all_on_same_contig = intervals.allequal((x.contig for x in context))

        # make an interval describing the start and stop of the query context
        minstart = min(x.start for x in context)
        maxstop  = max(x.stop  for x in context)
        query_bound = intervals.Interval(contig=anchor.contig, start=minstart, stop=maxstop)

        # do all the syntenic blocks within the target range map to regions within the query range?
        t = sorted(context, key=lambda x: (x.over.start))
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
    exo = exonerate.Exonerate(args.exonerate, genome=gen)
    res = ResultSet(gen=gen, syn=syn, exo=exo, width=10)
