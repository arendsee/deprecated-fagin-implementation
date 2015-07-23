#!/usr/bin/env python3

import fagin
import unittest

class TestSynteny(unittest.TestCase):
    def setUp(self):
        synshort = (
            't1 10  20  q1 10  20  .5 +',
            't1 30  40  q1 30  40  .5 +',
            't1 50  60  q1 50  60  .5 +',
        )
        self.syn = fagin.Synteny(synshort)
    def test_overlaps(self):
        self.assertTrue(self.syn.overlaps(fagin.Gene('q1', 25, 35, '.'), 1))
    def test_after(self):
        self.assertTrue(self.syn.before(fagin.Gene('q1', 10, 20, '.'), 1))
    def test_before(self):
        self.assertTrue(self.syn.after(fagin.Gene('q1', 50, 60, '.'), 1))

class TestBlock(unittest.TestCase):
    def setUp(self):
        self.block = fagin.Block('qqq', 10, 20, 'ttt', 100, 120, 0.5)
        self.before_gene = fagin.Gene('qqq', 1, 5, 'a')
        self.low_overlap_gene = fagin.Gene('qqq', 5, 12, 'a')
        self.inside_gene = fagin.Gene('qqq', 12, 18, 'a')
        self.high_overlap_gene = fagin.Gene('qqq', 15, 25, 'a')
        self.over_gene = fagin.Gene('qqq', 5, 25, 'a')
        self.after_gene = fagin.Gene('qqq', 30, 40, 'a')
    def test_overlaps(self):
        self.assertTrue(not self.block.overlaps(self.before_gene))
        self.assertTrue(self.block.overlaps(self.low_overlap_gene))
        self.assertTrue(self.block.overlaps(self.inside_gene))
        self.assertTrue(self.block.overlaps(self.high_overlap_gene))
        self.assertTrue(self.block.overlaps(self.over_gene))
        self.assertTrue(not self.block.overlaps(self.after_gene))
    def test_after(self):
        self.assertTrue(not self.block.after(self.before_gene))
        self.assertTrue(not self.block.after(self.low_overlap_gene))
        self.assertTrue(not self.block.after(self.inside_gene))
        self.assertTrue(not self.block.after(self.high_overlap_gene))
        self.assertTrue(not self.block.after(self.over_gene))
        self.assertTrue(self.block.after(self.after_gene))
    def test_before(self):
        self.assertTrue(self.block.before(self.before_gene))
        self.assertTrue(not self.block.before(self.low_overlap_gene))
        self.assertTrue(not self.block.before(self.inside_gene))
        self.assertTrue(not self.block.before(self.high_overlap_gene))
        self.assertTrue(not self.block.before(self.over_gene))
        self.assertTrue(not self.block.before(self.after_gene))

class TestContext(unittest.TestCase):
    def setUp(self):
        syndata = (
            't1 10  20  q1 10  20  .5 +',
            't1 30  40  q1 30  40  .5 +',
            't1 50  60  q1 50  60  .5 +',
            't1 50  60  q1 50  60  .5 +',
            't1 80  90  q1 80  90  .5 +',
            't1 100 120 q2 100 120 .5 +',
            't1 150 200 q2 150 200 .5 +',
            't1 250 270 q2 250 270 .5 +',
            't1 10  20  q3 10  20  .5 +',
            't1 30  40  q3 30  40  .5 +',
            't1 50  60  q3 50  60  .5 +'
        )
        synshort = (
            't1 10  20  q1 10  20  .5 +',
            't1 30  40  q1 30  40  .5 +',
            't1 50  60  q1 50  60  .5 +',
        )
        self.syn = fagin.Synteny(syndata)
        self.synone = fagin.Synteny(synshort)
    def test_synteny_overlap_TRUE(self):
        self.assertTrue(fagin.Context(syn=self.syn, gene=fagin.Gene('q1', 25,  35, '.')).match)
        self.assertTrue(fagin.Context(syn=self.syn, gene=fagin.Gene('q1', 32,  38, '.')).match)
        self.assertTrue(fagin.Context(syn=self.syn, gene=fagin.Gene('q1', 35,  45, '.')).match)
        self.assertTrue(fagin.Context(syn=self.syn, gene=fagin.Gene('q1', 25,  45, '.')).match)
        self.assertTrue(fagin.Context(syn=self.syn, gene=fagin.Gene('q1',  0, 200, '.')).match)
    def test_synteny_overlap_FALSE(self):
        self.assertTrue(not fagin.Context(syn=self.syn, gene=fagin.Gene('q1',  0, 5,   '.')).match)
        self.assertTrue(not fagin.Context(syn=self.syn, gene=fagin.Gene('q1', 22, 28,  '.')).match)
        self.assertTrue(not fagin.Context(syn=self.syn, gene=fagin.Gene('q1', 95, 100, '.')).match)
    def test_synteny_overlap_different_chromosomes(self):
        self.assertTrue(not fagin.Context(syn=self.syn, gene=fagin.Gene('q1', 175, 225, '.')).match)
        self.assertTrue(not fagin.Context(syn=self.syn, gene=fagin.Gene('q2',  22,  28, '.')).match)
    def test_synteny_edges(self):
        self.assertTrue(fagin.Context(syn=self.synone, gene=fagin.Gene('q1',  5, 15, '.'  )).match)
        self.assertTrue(fagin.Context(syn=self.syn,    gene=fagin.Gene('q1', 85, 95, '.'  )).match)

    def _context_lengths(self, context_obj, up, over, down):
        b = context_obj.query_context
        return(sum([x[0] == 'upstream'   for x in b]) == up   and
               sum([x[0] == 'overlap'    for x in b]) == over and
               sum([x[0] == 'downstream' for x in b]) == down)
    def test_query_context_overlapping(self):
        c = fagin.Context(syn=self.syn, gene=fagin.Gene('q1', 25, 35, '.'))
        self.assertTrue(self._context_lengths(context_obj=c, up=1, over=1, down=3))
    def test_query_context_non_overlapping(self):
        c = fagin.Context(syn=self.syn, gene=fagin.Gene('q1', 22, 28, '.'))
        self.assertTrue(self._context_lengths(context_obj=c, up=1, over=0, down=4))
    def test_query_context_width_down(self):
        c = fagin.Context(syn=self.syn, gene=fagin.Gene('q1', 0, 5, '.'), width=3)
        self.assertTrue(self._context_lengths(context_obj=c, up=0, over=0, down=3))
    def test_query_context_width_up(self):
        c = fagin.Context(syn=self.syn, gene=fagin.Gene('q1', 1000, 1005, '.'), width=3)
        self.assertTrue(self._context_lengths(context_obj=c, up=3, over=0, down=0))


if __name__ == '__main__':
    unittest.main()
