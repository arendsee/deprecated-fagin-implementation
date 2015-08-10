from lib.util import err
from lib.intervals import Interval

class Exonerate_hit:
    def __init__(self, row):
        self.name = row[0]

        if len(row) != 14:
            err("14 columns expected in exonerate output, {} detected" % len(row), file=sys.stderr)

        try:
            self.gene = Interval(contig=row[0],
                                start=int(row[1]),
                                stop=int(row[2]))

            self.target = Interval(contig=row[4],
                                start=int(row[5]),
                                stop=int(row[6]))
        except ValueError:
            err('start and stop values in exonerate file must be integers')

        try:
            (self.score,
             self.has_stop,
             self.has_frameshift,
             self.num_split_codons,
             self.num_intron,
             self.max_intron) = (int(x) for x in row[8:])
        except ValueError:
            err('score, has_stop, has_frameshift, num_split_codons, num_intron, max_intron columns in exonerate file must all be integers')

    def __str__(self):
        elements = (self.gene,
                    self.target,
                    self.score,
                    self.has_stop,
                    self.has_frameshift,
                    self.num_split_codons,
                    self.num_intron,
                    self.max_intron)
        out = '\t'.join((str(x) for x in elements))
        return(out)

class Exonerate:
    def __init__(self, _file, genome):
        self._file = _file

        # index GFF file by name
        self.gene_map = {g.name:g for g in genome.intervals()}

    def generator(self):
        # skip the header
        next(self._file)
        for line in self._file:
            row = line.split('\t')

            if not row[0] in self.gene_map:
                print("Gene '%s' in exonerate data missing in GFF file, skipping record", file=sys.stderr)
                continue

            yield Exonerate_hit(row=row)
