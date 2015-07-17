#!/usr/bin/env python3

import argparse
import sys
import itertools
import collections

__version__ = '0.0.0'

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


class IntervalSet:
    def __init__(self, c, s, e):
        self.intervals = dict()
        d = sorted(zip(c,s,e), key=lambda x: x[0])
        for k,g in itertools.groupby(d, lambda x: x[0]):
            start, end = tuple(zip(*g))[1:3]
            self.intervals[k] = Interval(start, end)

class Interval:
    def __init__(self, start, stop):
        try:
            self.start = toint(start)
            self.stop = toint(stop)
        except ValueError:
            err('start and stop indices must be integers')
        self.N = len(self.stop)
        self._validate()

    def _validate(self):
        if any((self.stop[i] < self.start[i] for i in range(self.N))):
            err('stop positions must be greater than start positions')

class Tabular:
    def __init__(self, fileobj, validate=True):
        rows = tuple(s.split() for s in fileobj if s[0] != '#')
        columns = tuple(zip(*rows))
        self.filename = fileobj.name
        self._assign_colnames(columns)
        if validate:
            self._validate_row_lengths(rows)
            self._validate_data()

    def _validate_row_lengths(self, rows):
        lengths = [len(s) for s in rows]
        all_equal_length = all([i == lengths[0] for i in lengths])
        if not all_equal_length:
            msg = 'Unequal rows of unequal length in "%s"'
            err(msg % self.filename)

    def _validate_data(self):
        pass

    def _assign_colnames(self):
        raise NotImplemented

class GFF(Tabular):
    def _assign_colnames(self, columns):
        try:
            self.type = columns[2]
            self.bounds = IntervalSet(c=columns[0], s=columns[3], e=columns[4])
            self.strand = columns[6]
            self.seqid = columns[8]
        except IndexError:
            err('GFF formated files must have 9 columns')

class Synblocks(Tabular):
    def _assign_colnames(self, columns):
        try:
            self.tbounds = IntervalSet(columns[0], columns[1], columns[2])
            self.qbounds = IntervalSet(columns[3], columns[4], columns[5])
            self.strand = columns[7]
        except IndexError:
            err('Synteny block file must have 8 columns')

        try:
            self.pident = tuple(float(x) for x in columns[6])
        except ValueError:
            err('In synteny block file, column 7 must be numeric')

    def _validate_data(self):
        if not all(0 <= x <= 1 for x in self.pident):
            err('In synteny file, proportion identity column (column 7) must be between 0 and 1')

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
    gff = GFF(args.gff)
    syn = Synblocks(args.synteny)
    nstrings = NStrings(args.nstring)
