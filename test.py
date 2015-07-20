#!/usr/bin/env python3

import fagin
import unittest

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

class TestSynteny(unittest.TestCase):
    def setUp(self):
        syndata = (
            't1 10  20  q1 10  20  .5 +',
            't1 30  40  q1 30  40  .5 +',
            't1 50  60  q1 50  60  .5 +',
            't1 50  60  q1 50  60  .5 +',
            't1 50  60  q1 50  60  .5 +',
            't1 50  60  q1 50  60  .5 +',
            't1 80  90  q1 80  90  .5 +',
            't1 100 120 q2 100 120 .5 +',
            't1 150 200 q2 150 200 .5 +',
            't1 250 270 q2 250 270 .5 +'
        )
        synshort = (
            't1 10  20  q1 10  20  .5 +',
            't1 30  40  q1 30  40  .5 +',
            't1 50  60  q1 50  60  .5 +',
        )
        self.syn = fagin.Synteny(syndata)
        self.synone = fagin.Synteny(synshort)
    def test_synteny_overlap_TRUE(self):
        self.assertTrue(self.syn.overlaps(fagin.Gene('q1', 25, 35, '.')))
        self.assertTrue(self.syn.overlaps(fagin.Gene('q1', 32, 38, '.')))
        self.assertTrue(self.syn.overlaps(fagin.Gene('q1', 35, 45, '.')))
        self.assertTrue(self.syn.overlaps(fagin.Gene('q1', 25, 45, '.')))
        self.assertTrue(self.syn.overlaps(fagin.Gene('q1', 0, 200, '.')))
    def test_synteny_overlap_FALSE(self):
        self.assertTrue(not self.syn.overlaps(fagin.Gene('q1', 0, 5, '.')))
        self.assertTrue(not self.syn.overlaps(fagin.Gene('q1', 22, 28, '.')))
        self.assertTrue(not self.syn.overlaps(fagin.Gene('q1', 95, 100, '.')))
    def test_synteny_overlap_different_chromosomes(self):
        self.assertTrue(not self.syn.overlaps(fagin.Gene('q1', 175, 225, '.')))
        self.assertTrue(not self.syn.overlaps(fagin.Gene('q2', 22, 28, '.')))
    def test_synteny_edges(self):
        self.assertTrue(self.synone.overlaps(fagin.Gene('q1', 5, 15, 'q')))
        self.assertTrue(self.syn.overlaps(fagin.Gene('q1', 85, 95, 'q')))


if __name__ == '__main__':
    unittest.main()
