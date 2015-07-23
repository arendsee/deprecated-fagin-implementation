#!/usr/bin/env python3

import argparse
import sys
import os
import itertools
import collections
import math
# from memory_profiler import profile

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

def err(msg):
    sys.exit(msg)

def toint(x):
    return tuple(int(i) for i in x)


class Tabular:
    def __init__(self, tab_data, validate=True):
        rows = tuple(s.split() for s in tab_data if s[0] != '#')
        self._load_rows(rows)
        if validate:
            self._validate_row_lengths(rows)
            self._validate_data()

    def _validate_row_lengths(self, rows):
        lengths = [len(s) for s in rows]
        all_equal_length = all([i == lengths[0] for i in lengths])
        if not all_equal_length:
            msg = 'Unequal rows of unequal length in synteny file'

    def _validate_data(self):
        pass

    def _load_rows(self):
        raise NotImplemented


class Genome(Tabular):
    def _load_rows(self, rows):
        try:
            rows = ((a, b, c, int(d), int(e), f, g, h, i) for a,b,c,d,e,f,g,h,i in rows)
            rows = sorted(rows, key=lambda x: (x[0], x[3]))
            self.genes = [Gene(*[r[i] for i in (0, 3, 4, 8)]) for r in rows]
        except IndexError:
            err('GFF formated files must have 9 columns')

    def add_context(self, syn, width=10):
        for gene in self.genes:
            gene.get_context(syn, width)

    def print_query_context(self, filename=sys.stdout):
        try:
            f = open(filename, 'w')
        except TypeError:
            f = filename
        for gene in self.genes:
            for line in gene.query_context_string():
                print(line, file=f)

    def __str__(self):
        return '\n'.join([str(g) for g in self.genes])

class Gene:
    def __init__(self, chrm, start, stop, name):
        self.chrm = chrm
        self.start = start
        self.stop = stop
        self.name = name
        self.context = None

    def get_context(self, syn, width=10):
        self.context = Context(gene=self, syn=syn, width=width)

    def query_context_string(self):
        for s in self.context.to_string():
            yield '%s\t%s' % (self.name, s)

    def __str__(self):
        return('%s' % (self.name))


class Synteny(Tabular):
    def _load_rows(self, rows):
        '''
        Load the synteny file. This file must have the following columns:
            1. target scaffold
            2. target start
            3. target stop
            4. query scaffold
            5. query start
            6. query stop
            7. proportion identity (0 <= x <= 1)
            8. orientation (+/-)
        All start and stop locations are indexed from 0.
        '''
        try:
            # convert to appropriate types
            rows = ((a, int(b), int(c), d, int(e), int(f), float(g), h) for a,b,c,d,e,f,g,h in rows)
            # load rows into Block objects
            self.blocks = [Block(*[r[i] for i in (3, 4, 5, 0, 1, 2, 6)]) for r in rows]
        except IndexError:
            err('Synteny block file must have 8 columns')
        except ValueError:
            err('Columns 1,2,4,5 of the synteny file must be integers, column 6 must be numeric')

        # create an linked list within the blocks ordered by target start sites
        self._index_target()
        # create an linked list within the blocks ordered by query start sites
        self._index_query()
        # order the list by query start sites
        self.blocks.sort(key = lambda x: (x.qchr, x.qstart, -x.qstop))
        # find query chromosome intervals
        self.qchr_intervals = self._calculate_chromosome_intervals()
        # number of blocks
        self.length = len(self.blocks)

    def _index_target(self):
        # sort by target chromosome and then by target start
        self.blocks.sort(key = lambda x: (x.tchr, x.tstart, -x.tstop))
        prior = self.blocks[0]
        for b in self.blocks[1:]:
            if prior.tchr == b.tchr:
                prior.tnext = b
                b.tprevious = prior
            prior = b

    def _index_query(self):
        # sort by query chromosome and then by target start
        self.blocks.sort(key = lambda x: (x.qchr, x.qstart, -x.qstop))
        prior = self.blocks[0]
        for b in self.blocks[1:]:
            if prior.qchr == b.qchr:
                prior.qnext = b
                b.qprevious = prior
            prior = b

    def _calculate_chromosome_intervals(self):
        qchr_intervals = dict()
        qchr = None
        for i, b in enumerate(self.blocks):
            if qchr == b.qchr:
                continue
            elif i == 0:
                qstart = 0
                qchr = b.qchr
            else:
                qchr_intervals[qchr] = (qstart, i-1)
                qstart = i
                qchr = b.qchr
        else:
            qchr_intervals[qchr] = (qstart, i)
        return(qchr_intervals)

    def _validate_data(self):
        if not all(0 <= b.pident <= 1 for b in self.blocks):
            err('In synteny file, proportion identity column (column 7) must be between 0 and 1')

        # Assert the mapping from chromosome region to block is consistent
        for k,ids in self.qchr_intervals.items():
            if not all(k == self.blocks[i].qchr for i in range(ids[0], ids[1]+1)):
                err('Internal error in Synteny class: chromosome mismatch')

    def _anchor(self, gene):
        '''
        Returns the index of a block overlapping the gene or, if the gene
        overlaps no syntenic block, an adjacent block in log(n) time.
        '''
        # the first and last indices bounding the gene's chromosome
        # (the Synteny constructor sorts the blocks tuple by chromosome)
        low, high = self.qchr_intervals[gene.chrm]

        # initial block index (blocks are sorted by start point)
        i = low + (high - low) // 2

        for _ in range(math.ceil(math.log2(high - low)) + 1):
            block = self.blocks[i]
            # If query endpoint is before the target startpoint
            if block.before(gene):
                high = i
                i = high - math.ceil((high - low) / 2)
            # If query startpoint is after the target endpoint
            elif block.after(gene):
                low = i
                i = low + math.ceil((high - low) / 2)
            # If neither of the above is true, the intervals must overlap
            else:
                return block
        else:
            return block

    def overlaps(self, gene, block_id):
        return(self.blocks[block_id].overlaps(gene=gene))

    def after(self, gene, block_id):
        return(self.blocks[block_id].after(gene=gene))

    def before(self, gene, block_id):
        return(self.blocks[block_id].before(gene=gene))

class Block:
    def __init__(self, qchr, qstart, qstop, tchr, tstart, tstop, pident):
        self.qchr = qchr
        self.qstart = qstart
        self.qstop = qstop
        self.tchr = tchr
        self.tstart = tstart
        self.tstop = tstop
        self.pident = pident
        self.tnext     = None # next block in target order
        self.tprevious = None # previous block in target order
        self.qnext     = None # next block in query order
        self.qprevious = None # previous block in query order

    def overlaps(self, gene):
        '''
        Return true if the gene and block overlap
        '''
        return(gene.start <= self.qstop and self.qstart <= gene.stop)

    def after(self, gene):
        '''
        Return true if the gene follows the block
        '''
        return(gene.start > self.qstop)

    def before(self, gene):
        '''
        Return true if the gene precedes the block
        '''
        return(gene.stop < self.qstart)

    def __str__(self):
        return '\t'.join((self.qchr,
                          str(self.qstart),
                          str(self.qstop),
                          self.tchr,
                          str(self.tstart),
                          str(self.tstop),
                          str(self.pident)))

class Context:
    def __init__(self, gene, syn, width=10):
        self.width = width
        anchor = syn._anchor(gene=gene)
        self.match = anchor.overlaps(gene)
        self.query_context = self._get_query_context(gene=gene, syn=syn, anchor=anchor)

    def _get_query_context(self, gene, syn, anchor):
        '''
        Given the index of one overlapping or adjacent gene, find the k
        upstream, all overlapping, and k downstream syntenic blocks. Where k is
        the context size set by the user (default=10)
        '''
        context = list()
        nup = 0
        block = anchor
        while nup != self.width and block:
            if block.overlaps(gene):
                context.append(('overlap', block))
            elif block.after(gene):
                context.append(('upstream', block))
                nup += 1
            block = block.qprevious

        ndown = 0
        block = anchor.qnext
        while ndown != self.width and block:
            if block.overlaps(gene):
                context.append(('overlap', block))
            elif block.before(gene):
                context.append(('downstream', block))
                ndown += 1
            block = block.qnext

        return(context)

    def _get_target_context(self, syn):
        '''
        Find the context of the homologous synteny region in the target
        '''
        raise NotImplemented

    def _classify_gene(self):
        '''
        Classify the gene into the following categories:
            1. simple matching
            2. simple missing
        '''
        raise NotImplemented

    def to_string(self):
        for b in self.query_context:
            yield '%s\t%s' % (b[0], str(b[1]))


class NStrings(Tabular):
    def _assign_colnames(self, columns):
        try:
            self.chr = columns[0]
            self.start = toint(columns[1])
            self.length = toint(columns[2])
        except IndexError:
            err('NStrings file must have 3 columns: chr_name, start, and length')
        except ValueError:
            err('in Nstrings file, columns 2 and 3 must be counting numbers')


if __name__ == '__main__':
    args = parse()
    try:
        os.mkdir(args.output_dir)
    except FileExistsError:
        if(os.listdir(args.output_dir)):
            err('Output directory must be empty')
    except PermissionError:
        err("You don't have permission to make directory '%s'" % args.output_dir)

    gen = Genome(args.gff)
    syn = Synteny(args.synteny)
    gen.add_context(syn, width=10)
    gen.print_query_context(os.path.join(args.output_dir, 'query_context.tab'))
