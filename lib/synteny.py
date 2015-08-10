import itertools
from lib.intervals import MappedInterval, IntervalSet
from lib.util import Tabular, err

class Synteny(Tabular):
    def _load_rows(self, rows):
        '''
        Load the synteny file. This file must have the following columns:
            0. query scaffold
            1. query start
            2. query stop
            3. target scaffold
            4. target start
            5. target stop
            6. proportion identity (0 <= x <= 1)
            7. orientation (+/-)
        All start and stop locations are indexed from 0.
        '''
        try:
            # convert to appropriate types
            qint = [MappedInterval(contig=a, start=int(b), stop=int(c)) for a,b,c,d,e,f,g,h in rows]
            tint = [MappedInterval(contig=d, start=int(e), stop=int(f)) for a,b,c,d,e,f,g,h in rows]
            scores = [float(g) for a,b,c,d,e,f,g,h in rows]
        except IndexError:
            err('Synteny block file must have 8 columns')
        except ValueError:
            err('Columns 1,2,4,5 of the synteny file must be integers, column 6 must be numeric')

        for q,t,s in zip(qint, tint, scores):
            q.over, t.over = t, q
            t.score, q.score = s, s
        self.query, self.target = IntervalSet(qint), IntervalSet(tint)

    def _validate_data(self):
        if not all(0 <= b.score <= 1 for b in itertools.chain(self.query.intervals(), self.target.intervals())):
            err('In synteny file, proportion identity column (column 7) must be between 0 and 1')

    def _missing_input_error(self):
        err('Synteny file is missing or unreadable')

    def anchor_query(self, interval):
        return(self.query.anchor(interval))

    def anchor_target(self, interval):
        return(self.target.anchor(interval))

    def cut_low_score_pairs(self, minscore):
        for interval in self.query:
            if interval.score < minscore:
                interval.remove()

