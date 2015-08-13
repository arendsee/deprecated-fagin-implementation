import sys
import lib.intervals as intervals

class HitMerger:
    def __init__(self, flank_width, min_neighbors, target_flank_ratio, quiet=False):
        self.flank_width = flank_width
        self.min_neighbors = min_neighbors
        self.target_flank_ratio = target_flank_ratio
        self.quiet = quiet

    def merge(self, result, hit, syn):
        '''
        Input one hit from an exonerate output.
        '''
        assert(hit.name == result.name)

        result.total_hits += 1

        anchor = syn.anchor_target(hit.target)

        # If no blocks map to the specified target contig, stop
        if not anchor:
            if not self.quiet:
                msg = "%s is on a contig with no syntenic blocks: %s"
                print(msg % (result.gene.name, str(hit.target)), file=sys.stderr)
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


