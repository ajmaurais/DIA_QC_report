
import unittest
from unittest import mock
import os
import re

import setup_functions
from setup_functions import TEST_DIR

from DIA_QC_report.submodules.read_metadata import read_metadata
from DIA_QC_report.submodules.dtype import Dtype


class TestFileTypeBase(setup_functions.AbstractTestsBase):
    META_TYPES = {'string_var': Dtype.STRING,
                  'bool_var': Dtype.BOOL,
                  'int_var': Dtype.INT,
                  'float_var': Dtype.FLOAT,
                  'na_var': Dtype.NULL}

    SKY_TYPES = {'string_var': 'text',
                 'bool_var': 'true_false',
                 'int_var': 'number',
                 'float_var': 'number',
                 'na_var': 'text'}


    def __init__(self):
        self.work_dir = None
        self.metadata_file = None
        self.result = None


    def test_is_sucessful(self):
        self.assertEqual(0, self.result.returncode)


    @mock.patch('DIA_QC_report.submodules.read_metadata.LOGGER', mock.Mock())
    def test_annotation_csv(self):
        annotation_csv = f'{self.work_dir}/sky_annotations.csv'
        self.assertTrue(os.path.isfile(annotation_csv))

        _, types_in = read_metadata(self.metadata_file, exclude_null_from_skyline=True)
        _, types_out = read_metadata(annotation_csv, exclude_null_from_skyline=False)

        # check that data types are correct
        self.assertDictEqual(self.META_TYPES, types_in)
        self.assertDictEqual(self.META_TYPES, types_out)


    def test_batch_file(self):
        annotation_bat = f'{self.work_dir}/sky_annotation_definitions.bat'
        self.assertTrue(os.path.isfile(annotation_bat))

        command_re = re.compile(r'^--annotation-name="([\w\.\-]+)" --annotation-targets=replicate --annotation-type=([a-z_]+)$')

        with open(annotation_bat, 'r') as inF:
            lines = inF.readlines()

        for line in lines:
            self.assertTrue((match := command_re.search(line)) is not None)
            self.assertTrue(match[1] in self.SKY_TYPES)
            self.assertEqual(self.SKY_TYPES[match[1]], match[2])


class TestTsvToSkylineCsv(unittest.TestCase, TestFileTypeBase):
    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{TEST_DIR}/work/test_tsv_to_skyline_csv'
        setup_functions.make_work_dir(cls.work_dir, clear_dir=True)
        cls.metadata_file = f'{TEST_DIR}/data/metadata/Strap_multi_var_metadata.tsv'

        commands = ['metadata_convert', cls.metadata_file]
        cls.result = setup_functions.run_command(commands, cls.work_dir, prefix=__name__)


class TestCsvToSkylineCsv(unittest.TestCase, TestFileTypeBase):
    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{TEST_DIR}/work/test_csv_to_skyline_csv'
        setup_functions.make_work_dir(cls.work_dir, clear_dir=True)
        cls.metadata_file = f'{TEST_DIR}/data/metadata/HeLa_metadata.csv'

        commands = ['metadata_convert', cls.metadata_file]
        cls.result = setup_functions.run_command(commands, cls.work_dir, prefix=__name__)


class TestJsonToSkylineCsv(unittest.TestCase, TestFileTypeBase):
    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{TEST_DIR}/work/test_json_to_skyline_csv'
        setup_functions.make_work_dir(cls.work_dir, clear_dir=True)
        cls.metadata_file = f'{TEST_DIR}/data/metadata/HeLa_metadata.json'

        commands = ['metadata_convert', cls.metadata_file]
        cls.result = setup_functions.run_command(commands, cls.work_dir, prefix=__name__)


class TestProblamaticHeaders(unittest.TestCase, TestFileTypeBase):

    META_TYPES = {'cellLine': Dtype.STRING,
                  'sample.name': Dtype.STRING,
                  'experiment': Dtype.STRING,
                  'NCI7.std': Dtype.BOOL}

    SKY_TYPES = {'cellLine': 'text',
                 'sample.name': 'text',
                 'experiment': 'text',
                 'NCI7.std': 'true_false'}

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{TEST_DIR}/work/test_bad_header_tsv_to_skyline_csv'
        setup_functions.make_work_dir(cls.work_dir, clear_dir=True)
        cls.metadata_file = f'{TEST_DIR}/data/metadata/Sp3_metadata.json'

        commands = ['metadata_convert', cls.metadata_file]
        cls.result = setup_functions.run_command(commands, cls.work_dir, prefix=__name__)
