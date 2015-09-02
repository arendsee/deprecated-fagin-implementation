#!/usr/bin/env python3

import fagin
import lib.intervals as intervals
import lib.genome as genome
import lib.syn_merger as syn_merger
import lib.synteny as synteny
import lib.result_manager as rMan
import unittest

class TestIntervals(unittest.TestCase):
    def setUp(self):
        self.a = intervals.Interval(contig='a', start=1, stop=9)
        self.b = intervals.Interval(contig='a', start=1, stop=10)
        self.c = intervals.Interval(contig='a', start=10, stop=11)
        self.d = intervals.Interval(contig='a', start=11, stop=19)
        self.e = intervals.Interval(contig='a', start=20, stop=21)
        self.f = intervals.Interval(contig='a', start=1, stop=30)
        self.g = intervals.Interval(contig='a', start=21, stop=30)
        self.h = intervals.Interval(contig='b', start=7, stop=14)

    def test_overlaps(self):
        a = intervals.Interval(contig='a', start=10, stop=20)
        self.assertFalse(intervals.overlaps(a, self.a))
        self.assertTrue(intervals.overlaps( a, self.b))
        self.assertTrue(intervals.overlaps( a, self.c))
        self.assertTrue(intervals.overlaps( a, self.d))
        self.assertTrue(intervals.overlaps( a, self.e))
        self.assertTrue(intervals.overlaps( a, self.f))
        self.assertFalse(intervals.overlaps(a, self.g))
        self.assertFalse(intervals.overlaps(a, self.h))
    def test_allequal(self):
        self.assertTrue(intervals.allequal([1,1,1,1]))
        self.assertTrue(intervals.allequal(['a', 'a', 'a', 'a']))

class TestIntervalSet(unittest.TestCase):
    def setUp(self):
        self.a = genome.Gene(contig='c1', start=1,  stop=10, name='a')
        self.b = genome.Gene(contig='c1', start=12, stop=15, name='b')
        self.c = genome.Gene(contig='c1', start=12, stop=17, name='c')
        self.d = genome.Gene(contig='c1', start=13, stop=16, name='d')
        self.e = genome.Gene(contig='c1', start=13, stop=16, name='e')
        self.f = genome.Gene(contig='c2', start=5,  stop=10, name='f')
        self.intset = intervals.IntervalSet([self.a, self.b, self.c, self.d, self.e, self.f])

        self.anchor_set = intervals.IntervalSet((
            genome.Gene(name='a', contig='a', start=10, stop=20),
            genome.Gene(name='b', contig='a', start=30, stop=40),
            genome.Gene(name='c', contig='a', start=50, stop=60),
            genome.Gene(name='d', contig='a', start=70, stop=80),
            genome.Gene(name='e', contig='b', start=0,  stop=5 ),
            genome.Gene(name='f', contig='c', start=5,  stop=10),
            genome.Gene(name='g', contig='c', start=15, stop=20)
        ))

    def test_sorting(self):
        i1 = intervals.IntervalSet([self.b, self.a, self.c, self.d, self.e, self.f])
        self.assertTrue(list(i1.intervals()) == list(self.intset.intervals()))

        i2 = intervals.IntervalSet([self.f, self.a, self.b, self.c, self.d, self.e])
        self.assertTrue(list(i2.intervals()) == list(self.intset.intervals()))

        i3 = intervals.IntervalSet([self.a, self.c, self.b, self.d, self.e, self.f])
        self.assertTrue(list(i3.intervals()) == list(self.intset.intervals()))

        i4 = intervals.IntervalSet([self.a, self.b, self.c, self.e, self.d, self.f])
        self.assertTrue(list(i4.intervals()) == list(self.intset.intervals()))

    def test_order(self):
        self.assertTrue(self.a.last == None)
        self.assertTrue(self.b.last.name == 'a')
        self.assertTrue(self.c.last.name == 'b')
        self.assertTrue(self.d.last.name == 'c')
        self.assertTrue(self.e.last.name == 'd')

        self.assertTrue(self.a.next.name == 'b')
        self.assertTrue(self.b.next.name == 'c')
        self.assertTrue(self.c.next.name == 'd')
        self.assertTrue(self.d.next.name == 'e')
        self.assertTrue(self.e.next == None)

        self.assertTrue(self.f.last == None)
        self.assertTrue(self.f.next == None)

    def test_anchor_for_overlap_cases(self):

        a  = self.anchor_set.anchor(intervals.Interval('a',  0, 5  ))
        b  = self.anchor_set.anchor(intervals.Interval('a', 21, 31 ))
        d  = self.anchor_set.anchor(intervals.Interval('a', 75, 95 ))
        cd = self.anchor_set.anchor(intervals.Interval('a', 61, 69 ))
        self.assertTrue(a.name == 'a')
        self.assertTrue(b.name == 'b')
        self.assertTrue(d.name == 'd')
        self.assertTrue(cd.name == 'c' or cd.name == 'd')

    def test_anchor_for_single_interval_contigs(self):
        e  = self.anchor_set.anchor(intervals.Interval('b', 90, 100))
        self.assertTrue(e.name == 'e')

    def test_anchor_for_two_interval_contigs(self):
        f  = self.anchor_set.anchor(intervals.Interval('c', 0, 2))
        fg = self.anchor_set.anchor(intervals.Interval('c', 12, 13))
        g  = self.anchor_set.anchor(intervals.Interval('c', 21, 22))
        self.assertTrue(f.name == 'f')
        self.assertTrue(fg.name == 'f' or fg.name == 'g')
        self.assertTrue(g.name == 'g')
    def test_get_preceding(self):
        self.assertTrue([x.name for x in intervals.get_preceding(self.c, 2)] == ['b', 'a'])
    def test_get_following(self):
        self.assertTrue([x.name for x in intervals.get_following(self.c, 2)] == ['d', 'e'])
    def test_get_overlapping_single(self):
        bound = intervals.Interval('c1', start=17, stop=20)
        self.assertTrue([x.name for x in self.intset.get_overlapping(bound)] == ['c'])
    def test_get_overlapping_none(self):
        bound = intervals.Interval('c1', start=20, stop=30)
        self.assertTrue([x.name for x in self.intset.get_overlapping(bound)] == [])
    def test_get_overlapping_multiple(self):
        bound = intervals.Interval('c1', start=11, stop=13)
        self.assertTrue(set([x.name for x in self.intset.get_overlapping(bound)]) == {'b', 'c', 'd', 'e'})

class TestContext(unittest.TestCase):
    def setUp(self):
        self.syn_simple = synteny.Synteny(rows = (
            ('q1', 10, 20, 'c1', 10, 20, 1, '+'),
            ('q1', 50, 60, 'c1', 50, 60, 1, '+')
        ))
        self.syn_not_simple = synteny.Synteny(rows = (
            ('q1', 10, 20, 'c2', 10, 20, 1, '+'),
            ('q1', 50, 60, 'c1', 50, 60, 1, '+')
        ))
        self.syn_not_simple_insertion = synteny.Synteny(rows = (
            ('q1',  10,  20, 'c1', 10, 20, 1, '+'),
            ('q1', 100, 180, 'c1', 30, 40, 1, '+'),
            ('q1',  50,  60, 'c1', 50, 60, 1, '+')
        ))
        self.missing_gene = genome.Gene(contig='q1', start=30, stop=40, name='missing')
        self.present_gene = genome.Gene(contig='q1', start=15, stop=25, name='present')

    def _get_result(self, gene, syn, width=1):
        synmer = syn_merger.SynMerger(width)
        result = rMan.Result(gene=gene)
        synmer.merge(result=result, syn=syn)
        return(result)

    def test_is_present(self):
        result = self._get_result(self.missing_gene, self.syn_simple)
        self.assertFalse(result.is_present)

        result = self._get_result(self.present_gene, self.syn_simple)
        self.assertTrue(result.is_present)

    def test_upper_and_lower(self):
        result = self._get_result(self.missing_gene, self.syn_simple)
        self.assertTrue(result.upper)
        self.assertTrue(result.lower)

        # Since the gene overlaps the lowest syntenic block, no non-overlapped
        # region is lower. Having this test prevents much mischief.
        result = self._get_result(self.present_gene, self.syn_simple)
        self.assertTrue(result.upper)
        self.assertFalse(result.lower)

    def test_is_simple(self):
        result = self._get_result(self.present_gene, self.syn_simple)
        self.assertTrue(result.is_simple)

        result = self._get_result(self.missing_gene, self.syn_simple)
        self.assertTrue(result.is_simple)

        result = self._get_result(self.missing_gene, self.syn_not_simple)
        self.assertFalse(result.is_simple)

        result = self._get_result(self.missing_gene, self.syn_not_simple_insertion)
        self.assertFalse(result.is_simple)


if __name__ == '__main__':
    unittest.main()
