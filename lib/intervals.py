import collections
import itertools
import math

def allequal(x):
    return(len(set(x)) == 1)

def overlaps(a, b):
    try:
        return (a.stop >= b.start) and (a.start <= b.stop)
    except TypeError:
        return False
    except AttributeError:
        return False

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

    def __eq__(self, other):
        return all((self.start == other.start,
                   self.stop  == other.stop,
                   self.contig == other.contig))

    def __ne__(self, other):
        return not self.__eq__(other)

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

