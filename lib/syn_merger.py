import itertools
import lib.intervals as intervals

class SynMerger:
    def __init__(self, width):
        self.width=width

    def merge(self, result, syn):
        self._syntenic_analysis(result=result, syn=syn)

    def _syntenic_analysis(self, result, syn):
        '''
        Analyzes the synteny data, setting the following variables
        1. is_present - does at least one syntenic block overlap the query gene?
        2. is_simple - all the up and downstream syntenic blocks on the same
              target contig and do they all map to the query region?
        3. lower - an internal variable storing the first non-overlapping block before the gene
        4. upper - an internal variable storing the first non-overlapping block after the gene
        '''
        anchor = syn.anchor_query(result.gene)
        if anchor:
            links = self._get_links(result=result, anchor=anchor)
            result.lower, result.upper = self._get_flanks(result=result, links=links, anchor=anchor)

            context = self._get_context(result=result, links=links)

            # does the gene overlap a syntenic block?
            result.is_present = bool(links)

            # are the upstream and downstream blocks in order on the same
            # chromosome?
            result.is_simple = self._get_is_simple(anchor=anchor, context=context, syn=syn)

    def _get_links(self, result, anchor):
        '''
        find contiguous synteny blocks in the query that all overlap gene
        '''
        m = anchor
        links = []
        while m and intervals.overlaps(m, result.gene):
            links.append(m)
            m = m.last
        m = anchor.next
        while m and intervals.overlaps(m, result.gene):
            links.append(m)
            m = m.next
        links = sorted(links, key=lambda x: (x.start, x.stop))
        return(links)

    def _get_flanks(self, result, links, anchor):
        '''
        Get the syntenic blocks flanking (but not overlapping) the gene
        '''
        if not links:
            if result.gene.stop < anchor.start:
                lower, upper = (anchor.last, anchor)
            else:
                lower, upper = (anchor, anchor.next)
        else:
            lower, upper = (links[0].last, links[-1].next)
        return((lower, upper))

    def _get_context(self, result, links):
        '''
        get all syntenic blocks that overlap the gene along with the WIDTH blocks up and down stream
        '''
        lower_context = intervals.get_preceding(result.lower, self.width - 1)
        upper_context = intervals.get_following(result.upper, self.width - 1)
        everything = itertools.chain(lower_context, [result.lower], links, [result.upper], upper_context)
        return([x for x in everything if x])

    def _get_is_simple(self, anchor, context, syn):
        # are all the intervals on the same contig?
        all_on_same_contig = intervals.allequal((x.over.contig for x in context))

        # make an interval describing the start and stop of the query context
        qminstart = min(x.start for x in context)
        qmaxstop  = max(x.stop  for x in context)
        tminstart = min(x.over.start for x in context)
        tmaxstop  = max(x.over.stop  for x in context)
        query_bound  = intervals.Interval(contig=anchor.contig, start=qminstart, stop=qmaxstop)
        target_bound = intervals.Interval(contig=anchor.over.contig, start=tminstart, stop=tmaxstop)
        has_outer = False
        for q in syn.target.get_overlapping(target_bound):
            if not intervals.overlaps(query_bound, q.over):
                has_outer = True
                break

        print([all_on_same_contig, has_outer])

        is_simple = all_on_same_contig and not has_outer
        return(is_simple)

