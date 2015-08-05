import src.util as util

class NStrings(util.Tabular):
    def _assign_colnames(self, columns):
        try:
            self.chr = columns[0]
            self.start = tuple(int(i) for i in columns[1])
            self.length = tuple(int(i) for i in columns[2])
        except IndexError:
            err('NStrings file must have 3 columns: chr_name, start, and length')
        except ValueError:
            err('in Nstrings file, columns 2 and 3 must be counting numbers')
