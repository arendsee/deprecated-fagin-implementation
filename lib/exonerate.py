from lib.util import err
from lib.intervals import Interval

class Hit:
    '''
    The base class for hit objects. Contains the query interval, the target interval, and the score.
    '''
    def __init__(self, row):
        self.name = row[0]

        try:
            self.gene = Interval(contig=row[0],
                                start=int(row[1]),
                                stop=int(row[2]))

            self.target = Interval(contig=row[4],
                                start=int(row[5]),
                                stop=int(row[6]))
        except ValueError:
            err('start and stop values in exonerate file must be integers')
        except IndexError:
            err('expected 9 columns')

        try:
            self.score = float(row[8])
        except ValueError:
            err('the score column (9) must be numeric')

    def __str__(self):
        elements = (self.gene, self.target, self.score)
        out = '\t'.join((str(x) for x in elements))
        return(out)

    def __eq__(self, other):
        return all((self.target == other.target,
                    self.score  == other.score,
                    self.gene   == other.gene))

    def __ne__(self, other):
        return not self.__eq__(other)

class IntronHit(Hit):
    def __init__(self, row):
        try:
            super().__init__(row[0:9])
        except IndexError:
            err('hit data must have at least 9 columns')

        try:
            (self.first_stop,
                self.has_frameshift,
                self.num_split_codons,
                self.num_intron,
                self.max_intron) = (int(x) for x in row[9:14])
        except ValueError:
            err('first_stop, has_frameshift, num_split_codons, num_intron, max_intron columns in exonerate file must all be integers')

        # try:
        #     self.intron_lengths

    def __str__(self):
        elements = (self.gene,
                    self.target,
                    self.score,
                    self.first_stop,
                    self.has_frameshift,
                    self.num_split_codons,
                    self.num_intron,
                    self.max_intron)
        out = '\t'.join((str(x) for x in elements))
        return(out)

class Exonerate:
    def __init__(self, _file):
        self._file = _file

    def generator(self):
        # skip the header
        header = next(self._file).split('\t')
        if(len(header) == 8):
            for line in self._file:
                row = line.split('\t')
                yield Hit(row=row)
        elif(len(header) == 14):
            for line in self._file:
                row = line.split('\t')
                yield IntronHit(row=row)
        else:
            err('Unrecognized hit input (incorrect number of columns)')
