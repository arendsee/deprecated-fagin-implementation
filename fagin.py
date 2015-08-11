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

        # --- variables set in _syntenic_analysis ---
        # -------------------------------------------
        self.is_present = False
        self.is_simple  = False
        self.lower = None  # the block downstream of the last homologous block
        self.upper = None  # the block ustream of the last homologous block
        self._syntenic_analysis(syn=syn, width=width)

        # --- variables modified when adding exonerate hits ---
        # -----------------------------------------------------
        # declaration of an interval containing the flanks around the gene
        # in the query
        self.query_flanks = None
        self.hits = []

    def add_exonerate_hit(self, hit, syn, flank_width=25000):
        '''
        Input one hit from an exonerate output.
        '''
        assert(hit.name == self.name)

        if not self.query_flanks:
            self.query_flanks = intervals.Interval(
                contig=self.gene.contig,
                start=max(0, self.gene.start - flank_width),
                stop=(self.gene.stop + flank_width))

        anchor = syn.anchor_query(hit.target)
        if anchor:
            target_flanks = intervals.Interval(
                contig=anchor.contig,
                start=max(0, anchor.start - 2 * flank_width),
                stop=(anchor.stop + 2 * flank_width))
        else:
            print("%s has no match" % self.gene.name)

    def _syntenic_analysis(self, syn, width):
        '''
        Analyzes the synteny data, setting the following variables
        1. is_present - does at least one syntenic block overlap the query gene?
        2. is_simple - all the up and downstream syntenic blocks on the same
              target contig and do they all map to the query region?
        3. lower - an internal variable storing the first non-overlapping block before the gene
        4. upper - an internal variable storing the first non-overlapping block after the gene
        '''
        anchor = syn.anchor_query(self.gene)
        if anchor:
            links = self._get_links(anchor=anchor)
            self.lower, self.upper = self._get_flanks(links=links, anchor=anchor)

            context = self._get_context(anchor=anchor, links=links, width=width)

            # does the gene overlap a syntenic block?
            is_present = bool(links)

            # are the upstream and downstream blocks in order on the same
            # chromosome?
            is_simple = self._get_is_simple(anchor=anchor, context=context)

        self.is_present = is_present
        self.is_simple  = is_simple

    def _get_links(self, anchor):
        '''
        find contiguous synteny blocks in the query that all overlap gene
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

    def _get_flanks(self, links, anchor):
        '''
        Get the syntenic blocks flanking (but not overlapping) the gene
        '''
        if not links:
            if self.gene.stop < anchor.start:
                lower, upper = (anchor.last, anchor)
            else:
                lower, upper = (anchor, anchor.next)
        else:
            lower, upper = (links[0].last, links[-1].next)
        return((lower, upper))

    def _get_context(self, anchor, links, width):
        '''
        get all blocks on the query between, but not including, the flanking genes
        '''
        lower_context = intervals.get_preceding(self.lower, width)
        upper_context = intervals.get_following(self.upper, width)
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
