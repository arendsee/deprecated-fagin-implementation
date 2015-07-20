#!/usr/bin/env python3

import argparse
import sys
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

    def add_synteny_annotations(self, syn):
        for gene in self.genes:
            gene.add_synteny(syn)

    def __str__(self):
        return '\n'.join([str(g) for g in self.genes])

class Gene:
    def __init__(self, chrm, start, stop, name):
        self.chrm = chrm
        self.start = start
        self.stop = stop
        self.name = name
        self.synteny = None

    def add_synteny(self, syn):
        self.synteny = GeneSynteny(self, syn)

    def __str__(self):
        return('%s %s' % (self.name, self.synteny.matched))

class GeneSynteny:
    def __init__(self, gene, syn):
        self.gene = gene
        self.matched = syn.overlaps(gene)

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
            # sort by query chromosome and then by query start
            rows = sorted(rows, key=lambda x: (x[3], x[4]))
            # load rows into Block objects
            self.blocks = tuple(Block(*[r[i] for i in (3, 4, 5, 0, 1, 2, 6)]) for r in rows)
        except IndexError:
            err('Synteny block file must have 8 columns')
        except ValueError:
            err('Columns 1,2,4,5 of the synteny file must be integers, column 6 must be numeric')

        self.qchr_intervals = dict()
        self.tchr_map = dict()
        self._calculate_chromosome_intervals()
        self.length = len(self.blocks)

    def _calculate_chromosome_intervals(self):
        qchr = None
        for i, b in enumerate(self.blocks):
            self.tchr_map[i] = b.tchr
            if qchr == b.qchr:
                continue
            elif i == 0:
                qstart = 0
                qchr = b.qchr
            else:
                self.qchr_intervals[qchr] = (qstart, i-1)
                qstart = i
                qchr = b.qchr
        else:
            self.qchr_intervals[qchr] = (qstart, i)

    def _validate_data(self):
        if not all(0 <= b.pident <= 1 for b in self.blocks):
            err('In synteny file, proportion identity column (column 7) must be between 0 and 1')

        # Assert the mapping from chromosome region to block is consistent
        for k,ids in self.qchr_intervals.items():
            if not all(k == self.blocks[i].qchr for i in range(ids[0], ids[1]+1)):
                err('Internal error in Synteny class: chromosome mismatch')

    def overlaps(self, gene):
        '''
        Return true if input gene overlaps a synteny block
        '''
        low, high = self.qchr_intervals[gene.chrm]
        block_id = low + (high - low) // 2
        max_iterations = math.ceil(math.log2(high - low)) + 1
        for i in range(max_iterations):
            block = self.blocks[block_id]
            # If query endpoint is before the target startpoint
            if block.before(gene):
                high = block_id
                block_id = high - math.ceil((high - low) / 2)
            # If query startpoint is after the target endpoint
            elif block.after(gene):
                low = block_id
                block_id = low + math.ceil((high - low) / 2)
            # If neither of the above is true, the intervals must overlap
            else:
                return True
        return False

class Block:
    def __init__(self, qchr, qstart, qstop, tchr, tstart, tstop, pident):
        self.qchr = qchr
        self.qstart = qstart
        self.qstop = qstop
        self.tchr = tchr
        self.tstart = tstart
        self.tstop = tstop
        self.pident = pident

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
    gen = Genome(args.gff)
    syn = Synteny(args.synteny)
    gen.add_synteny_annotations(syn)
    print(gen)
    # nstrings = NStrings(args.nstring)
