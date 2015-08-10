import collections
from lib.intervals import OrderedInterval, IntervalSet
from lib.util import Tabular, err

class Gene(OrderedInterval):
    def __init__(self, name, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = name

class Genome(Tabular, IntervalSet):
    def __init__(self, filename):
        Tabular.__init__(self, filename)
        IntervalSet.__init__(self, self.genes)

    def _load_rows(self, rows):
        try:
            self.genes = (Gene(name=i, contig=a, start=int(d), stop=int(e)) for a,b,c,d,e,f,g,h,i in rows)
        except IndexError:
            err('GFF formated files must have 9 columns')
        except ValueError:
            err("Start and stop positions must be integers")

    def _missing_input_error(self):
        err('GFF file is missing or unreadable')

    def __str__(self):
        rows = ((g.name, g.contig, str(g.start), str(g.stop)) for g in self.intervals())
        out = '\n'.join(['\t'.join(x) for x in rows])
        return(out)

