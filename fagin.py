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

def allequal(x):
    return(len(set(x)) == 1)

def overlaps(a, b):
    return((a.stop >= b.start) and (a.start <= b.stop))

def get_preceding(ordint, n):
    '''
    Retrieve the preceding n intervals
    '''
    for _ in range(n):
        try:
            ordint = ordint.last
            yield ordint
        except AttributeError:
            break

def get_following(ordint, n):
    '''
    Retrieve the following n intervals
    '''
    for _ in range(n):
        try:
            ordint = ordint.next
            yield ordint
        except AttributeError:
            break

class Interval:
    def __init__(self, contig, start, stop):
        self.start = start
        self.stop = stop
        self.contig = contig

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

        # return None if the contig of the input is not in the interval set
        if other.contig in self.contigs:
            con = self.contigs[other.contig]
        else:
            return None

        # original bounds of the search space (these will be tightened until a
        # solution is found
        low, high = 0, len(con) - 1

        # starting index (start in the middle of the vector)
        i = high // 2

        # maximum number of steps to find solution
        steps = math.ceil(math.log2(high - low)) + 1

        # binary search
        for _ in range(steps):
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

class OrderedInterval:
    def __init__(self, contig, start, stop, last, next):
        self.start = start
        self.stop = stop
        self.contig = contig
        self.last = last
        self.next = next

    @classmethod
    def remove(cls, obj):
        if obj.last:
            obj.last.next = obj.next
        obj.next.last = obj.last
        del obj

class MappedInterval:
    def __init__(self, contig, start, stop, last=None, next=None, over=None, score=None):
        self.start = start
        self.stop = stop
        self.contig = contig
        self.last = last
        self.next = next
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

class Gene:
    def __init__(self, contig, start, stop, last=None, next=None, name=None):
        self.start = start
        self.stop = stop
        self.contig = contig
        self.last = last
        self.next = next
        self.name = name

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

    def analyze(self, syn):
        for g in self.intervals():
            yield Result(gene=g, syn=syn)

    def __str__(self):
        rows = ((g.name, g.contig, str(g.start), str(g.stop)) for g in self.intervals())
        out = '\n'.join(['\t'.join(x) for x in rows])
        return(out)

class Synteny(Tabular):
    def _load_rows(self, rows):
        '''
        Load the synteny file. This file must have the following columns:
            0. query scaffold
            1. query start
            2. query stop
            3. target scaffold
            4. target start
            5. target stop
            6. proportion identity (0 <= x <= 1)
            7. orientation (+/-)
        All start and stop locations are indexed from 0.
        '''
        try:
            # convert to appropriate types
            qint = [MappedInterval(contig=a, start=int(b), stop=int(c)) for a,b,c,d,e,f,g,h in rows]
            tint = [MappedInterval(contig=d, start=int(e), stop=int(f)) for a,b,c,d,e,f,g,h in rows]
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
        while m and overlaps(m, self.gene):
            links.append(m)
            m = m.last
        m = self.anchor.next
        while m and overlaps(m, self.gene):
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
        lower_context = get_preceding(lower, self.width)
        upper_context = get_following(upper, self.width)
        everything = [x for x in itertools.chain(lower_context, self.links, upper_context) if x]

        return(everything)

    def _get_is_simple(self):
        # are all the intervals on the same contig?
        all_on_same_contig = allequal((x.contig for x in self.context))

        # make an interval describing the start and stop of the query context
        minstart = min(x.start for x in self.context)
        maxstop  = max(x.stop  for x in self.context)
        query_bound = Interval(contig=self.anchor.contig, start=minstart, stop=maxstop)

        # do all the syntenic blocks within the target range map to regions within the query range?
        t = sorted(self.context, key=lambda x: (x.over.start))
        q = t[0]
        has_outer = False
        while True:
            if not q:
                break
            elif q.start > t[-1].stop:
                break
            elif not overlaps(q.over, query_bound):
                has_outer = True
                break
            else:
                q = q.next

        is_simple = all_on_same_contig and not has_outer
        return(is_simple)

    def __str__(self):
        out = '\t'.join((self.gene.name, str(self.is_present), str(self.is_simple)))
        return(out)

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
    g = list(gen.analyze(syn))
    for result in gen.analyze(syn):
        print(result)
