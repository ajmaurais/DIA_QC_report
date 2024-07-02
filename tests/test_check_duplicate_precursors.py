
import unittest
import pandas as pd

import setup_functions
from setup_functions import TEST_DIR

from DIA_QC_report.parse_data import check_duplicate_precursors

class TestCheckDuplicatePrecursorsFxn(unittest.TestCase):
    def setUp(self):
        project = 'Strap'
        fname = f'{TEST_DIR}/data/skyline_reports/{project}_by_protein_precursor_quality.tsv'
        self.precursors = pd.read_csv(fname, sep='\t')

    def test_invalid_mode_fails(self):
        pass

