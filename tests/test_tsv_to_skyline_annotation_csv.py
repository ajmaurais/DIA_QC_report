
import unittest
import os
import pandas as pd
import re

import setup_functions
from setup_functions import TEST_DIR

from DIA_QC_report.submodules.dia_db_utils import read_metadata
from DIA_QC_report.submodules.metadata import Dtype


class TestTsvToSkylineCsv(unittest.TestCase):

    META_TYPES = {'string_var': Dtype.STRING,
                  'bool_var': Dtype.BOOL,
                  'int_var': Dtype.INT,
                  'float_var': Dtype.FLOAT,
                  'na_var': Dtype.BOOL}

    SKY_TYPES = {'string_var': 'text',
                 'bool_var': 'true_false',
                 'int_var': 'number',
                 'float_var': 'number',
                 'na_var': 'text'}

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{TEST_DIR}/work/test_tsv_to_skyline_csv'
        setup_functions.make_work_dir(cls.work_dir, clear_dir=True)
        cls.metadata_tsv = f'{TEST_DIR}/data/metadata/Strap_multi_var_metadata.tsv'

        commands = ['tsv_to_skyline_annotation_csv', cls.metadata_tsv]
        cls.result = setup_functions.run_command(commands, cls.work_dir)


    def test_is_sucessful(self):
        self.assertEqual(0, self.result.returncode)


    def test_annotation_csv(self):
        annotation_csv = f'{self.work_dir}/sky_annotations.csv'
        self.assertTrue(os.path.isfile(annotation_csv))
        
        import pudb
        # pudb.set_trace()
        
        metadata_in, types_in = read_metadata(self.metadata_tsv)
        metadata_out, types_out = read_metadata(annotation_csv)

        # check that data types are correct
        self.assertDictEqual(self.META_TYPES, types_in)
        self.assertDictEqual(self.META_TYPES, types_out)


    def test_batch_file(self):
        annotation_bat = f'{self.work_dir}/sky_annotation_definitions.bat'
        self.assertTrue(os.path.isfile(annotation_bat))

        command_re = re.compile(r'^--annotation-name="(\w+)" --annotation-targets=replicate --annotation-type=([a-z_]+)$')

        import pudb
        # pudb.set_trace()

        with open(annotation_bat, 'r') as inF:
            lines = inF.readlines()

        for line in lines:
            self.assertTrue((match := command_re.search(line)) is not None)
            self.assertTrue(match[1] in self.SKY_TYPES)
            if self.SKY_TYPES[match[1]] != match[2]:
                print('FUCK!')
            self.assertEqual(self.SKY_TYPES[match[1]], match[2])

