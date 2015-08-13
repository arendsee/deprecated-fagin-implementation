import sys
import lib.intervals as intervals

class HitMerger:
    def __init__(self, flank_width=25000, min_neighbors=3, winnow=False, target_flank_ratio=2):
        self.flank_width = flank_width
        self.min_neighbors = min_neighbors
        self.winnow = winnow
        self.target_flank_ratio = target_flank_ratio

    def merge(self, result, hit, syn):
        '''
        Input one hit from an exonerate output.
        '''
        assert(hit.name == result.name)

        result.total_hits += 1

        # If this hit is inferior to a prior hit, stop
        if self.winnow:
            return False

        anchor = syn.anchor_target(hit.target)

        # If no blocks map to the specified target contig, stop
        if not anchor:
            print("%s is on a contig with no syntenic blocks: %s" % (result.gene.name, str(hit.target)), file=sys.stderr)
            return False

        # Define the interval in which to search for neighbors in the target
        target_flanks = intervals.Interval(
            contig=anchor.contig,
            start=max(0, anchor.start - self.target_flank_ratio * self.flank_width),
            stop=(anchor.stop + self.target_flank_ratio * self.flank_width))

        # If the query flanks have not yet been measured, do so
        if not result.query_flanks:
            result.query_flanks = intervals.Interval(
                contig=result.gene.contig,
                start=max(0, result.gene.start - self.flank_width),
                stop=(result.gene.stop + self.flank_width))

        # === my hacky first order solution ===

        # I simply count the number of nearby (by an arbitrary cutoff) blocks
        # that map to a region near the target interval. If there are more than
        # a certain number, I keep the exonerate hit.

        matching, total = (0, 0)
        a, b = (result.lower, result.upper)

        while intervals.overlaps(result.query_flanks, a):
            total += 1
            matching += intervals.overlaps(target_flanks, a.over)
            a = a.last

        while intervals.overlaps(result.query_flanks, b):
            total += 1
            matching += intervals.overlaps(target_flanks, b.over)
            b = b.next

        if matching >= self.min_neighbors:
            result.hits.append(hit)

    def _winnow(result, hit):
        '''
        skip if this hit overlaps an existing hit and the existing hit
        1) has a higher score
        2) does not have a stop or frameshift unless this hit does as well
        3) the longest intron is under 100000 bases unless this hit is as well
        '''
        for i,h in enumerate(result.hits):
            if intervals.overlaps(h.target, hit.target) and   \
               hit.score <= h.score and                       \
               (not (h.has_stop or h.has_frameshift) or       \
                    (hit.has_stop or hit.has_frameshift)) and \
               (h.max_intron < 1e5 or hit.max_intron >= 1e5):
                return True

