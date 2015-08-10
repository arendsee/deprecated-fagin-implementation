import sys

def err(msg):
    sys.exit(msg)

class Tabular:
    def __init__(self, tab_data, validate=True):
        try:
            rows = tuple(s.split() for s in tab_data if s[0] != '#')
        except TypeError:
            self._missing_input_error()
        self._load_rows(rows)
        if validate:
            self._validate_row_lengths(rows)
            self._validate_data()

    def _validate_row_lengths(self, rows):
        lengths = [len(s) for s in rows]
        all_equal_length = all([i == lengths[0] for i in lengths])
        if not all_equal_length:
            msg = 'Unequal rows of unequal length in synteny file'

    def _validate_data(self):
        pass

    def _load_rows(self):
        raise NotImplemented

    def _missing_input_error(self):
        raise NotImplemented
