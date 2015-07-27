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


class Interval:
    def __init__(self, contig, start, stop):
        self.start = start
        self.stop = stop
        self.contig = contig

    def overlaps(self, other):
        return((other.stop > self.start) and (other.start <= self.stop))

    def __str__(self):
        return('\t'.join((self.contig, str(self.start), str(self.stop))))

class IntervalSet:
    def __init__(self, intervals):
        intervals = sorted(intervals, key = lambda x: (x.contig, x.start, x.stop))
        self.contigs = collections.defaultdict(list)
        for interval in intervals:
            self.contigs[interval.contig].append(interval)
        for group in self.contigs.values():
            prior = group[0]
            for this in group[1:]:
                prior.next = this
                this.last = prior
                prior = this

    def intervals(self):
        for interval in itertools.chain(*self.contigs.values()):
            yield interval

    def anchor(self, other):
        '''
        Returns the index of a block overlapping the gene or, if the gene
        overlaps no syntenic block, an adjacent block in log(n) time.
        '''
        con = self.contigs[other.contig]
        low, high = 0, len(con) - 1
        i = high // 2
        for _ in range(math.ceil(math.log2(high - low)) + 1):
            this = con[i]
            if other.stop < this.start:
                high = i
                i = high - math.ceil((high - low) / 2)
            elif other.start > this.stop:
                low = i
                i = low + math.ceil((high - low) / 2)
            else:
                return this
        else:
            return this

class OrderedInterval(Interval):
    def __init__(self, last=None, next=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.last = last
        self.next = next

    @classmethod
    def remove(cls, obj):
        if obj.last:
            obj.last.next = obj.next
        obj.next.last = obj.last
        del obj

class MappedInterval(OrderedInterval):
    def __init__(self, over=None, score=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.over = over
        self.score = score

    @classmethod
    def remove(cls, obj):
        if obj.last:
            obj.last.next = obj.next
            obj.over.last.next = obj.over.next
        obj.next.last = obj.last
        obj.over.next.last = obj.over.last
        del obj.over
        del obj

class Gene(OrderedInterval):
    def __init__(self, name, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = name
        self.links = []

    def add_links(self, syn):
        anchor = syn.anchor_query(interval=self)
        m = anchor
        while m and m.overlaps(self):
            self.links.append(m)
            m = m.last
        m = anchor.next
        while m and m.overlaps(self):
            self.links.append(m)
            m = m.next

    def print_links(self):
        if self.links:
            rows = ['\t'.join(self.name, str(link)) for link in self.links]
            print('\n'.join(rows))


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

class Genome(Tabular, IntervalSet):
    def __init__(self, filename):
        Tabular.__init__(self, filename)
        IntervalSet.__init__(self, self.genes)

    def _load_rows(self, rows):
        try:
            self.genes = (Gene(name=i, contig=a, start=int(d), stop=int(e)) for a,b,c,d,e,f,g,h,i in rows)
        except IndexError:
            err('GFF formated files must have 9 columns')
        except ValueError:
            err("Start and stop positions must be integers")

    def add_links(self, syn):
        for g in self.intervals():
            g.add_links(syn)

    def __str__(self):
        rows = ((g.name, g.contig, str(g.start), str(g.stop)) for g in self.intervals())
        out = '\n'.join(['\t'.join(x) for x in rows])
        return(out)

class Synteny(Tabular):
    def _load_rows(self, rows):
        '''
        Load the synteny file. This file must have the following columns:
            0. target scaffold
            1. target start
            2. target stop
            3. query scaffold
            4. query start
            5. query stop
            6. proportion identity (0 <= x <= 1)
            7. orientation (+/-)
        All start and stop locations are indexed from 0.
        '''
        try:
            # convert to appropriate types
            qint = [MappedInterval(contig=d, start=int(e), stop=int(f)) for a,b,c,d,e,f,g,h in rows]
            tint = [MappedInterval(contig=a, start=int(b), stop=int(c)) for a,b,c,d,e,f,g,h in rows]
            scores = [float(g) for a,b,c,d,e,f,g,h in rows]
        except IndexError:
            err('Synteny block file must have 8 columns')
        except ValueError:
            err('Columns 1,2,4,5 of the synteny file must be integers, column 6 must be numeric')

        for q,t,s in zip(qint, tint, scores):
            q.over, t.over = t, q
            t.score, q.score = s, s
        self.query, self.target = IntervalSet(qint), IntervalSet(tint)

    def _validate_data(self):
        if not all(0 <= b.score <= 1 for b in itertools.chain(self.query.intervals(), self.target.intervals())):
            err('In synteny file, proportion identity column (column 7) must be between 0 and 1')

    def anchor_query(self, interval):
        return(self.query.anchor(interval))

    def anchor_target(self, interval):
        return(self.target.anchor(interval))

    def cut_low_score_pairs(self, minscore):
        for interval in self.query:
            if interval.score < minscore:
                interval.remove()


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
    gen.add_links(syn)
    for gene in gen.intervals():
        print(gene.name, len(gene.links))
