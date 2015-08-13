import itertools
import lib.intervals as intervals

class SynMerger:
    def __init__(self, width=1):
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

            context = self._get_context(result=result, anchor=anchor, links=links)

            # does the gene overlap a syntenic block?
            result.is_present = bool(links)

            # are the upstream and downstream blocks in order on the same
            # chromosome?
            result.is_simple = self._get_is_simple(anchor=anchor, context=context)

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

    def _get_context(self, result, anchor, links):
        '''
        get all blocks on the query between, but not including, the flanking genes
        '''
        lower_context = intervals.get_preceding(result.lower, self.width)
        upper_context = intervals.get_following(result.upper, self.width)
        everything = [x for x in itertools.chain(lower_context, links, upper_context) if x]
        return(everything)

    def _get_is_simple(self, anchor, context):
        # are all the intervals on the same contig?
        all_on_same_contig = intervals.allequal((x.contig for x in context))

        # make an interval describing the start and stop of the query context
        minstart = min(x.start for x in context)
        maxstop  = max(x.stop  for x in context)
        query_bound = intervals.Interval(contig=anchor.contig, start=minstart, stop=maxstop)

        # do all the syntenic blocks within the target range map to regions within the query range?
        t = sorted(context, key=lambda x: (x.over.start))
        q = t[0]
        has_outer = False
        while True:
            if not q:
                break
            elif q.start > t[-1].stop:
                break
            elif not intervals.overlaps(q.over, query_bound):
                has_outer = True
                break
            else:
                q = q.next

        is_simple = all_on_same_contig and not has_outer
        return(is_simple)

