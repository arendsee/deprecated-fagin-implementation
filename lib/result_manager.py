class ResultManager:
    def __init__(self, gen, syn, exo, syn_merger, hit_merger, hit_analyzer):
        self.results = {g.name: Result(g) for g in gen.intervals()}

        # merge in the synteny data
        for result in self.results.values():
            syn_merger.merge(result=result, syn=syn)

        # merge in the exonerate hit data
        for hit in exo.generator():
            try:
                hit_merger.merge(result=self.results[hit.name], hit=hit, syn=syn)
            except KeyError:
                err('The gene {} in the hit file is missing from the gff file' % hit.name)

        # hit_analyzer.filter(self.results)

    def get(self, name):
        return self.results[name]

    def write(self):
        for r in self.results.values():
            s = str(r)
            if s:
                print(s)

class Result:
    '''
    Stores all results and handles record output
    '''
    def __init__(self, gene):
        self.name = gene.name
        self.gene = gene

        # --- variables set by SynMerger ---
        # -------------------------------------------
        self.is_present = False
        self.is_simple  = False
        self.lower = None # the block downstream of the last homologous block
        self.upper = None # the block upstream of the last homologous block

        # --- variables set by HitMerger ---
        # -------------------------------------------
        # declaration of an interval containing the flanks around the gene
        # in the query
        self.query_flanks = None
        self.total_hits = 0
        self.hits = []

    def __str__(self):
        # out = '\t'.join((self.gene.name, str(self.is_present), str(self.is_simple)))
        out = '\n'.join([str(h) for h in self.hits])
        return(out)
