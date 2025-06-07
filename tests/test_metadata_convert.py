
import unittest
from unittest import mock
import os
import re

import setup_functions
from setup_functions import TEST_DIR

from DIA_QC_report.submodules.read_metadata import Metadata
from DIA_QC_report.submodules.dtype import Dtype
from DIA_QC_report import metadata_convert

PROG = 'dia_qc metadata_convert'

class TestOptions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.metadata_dir = f'{TEST_DIR}/data/metadata'
        cls.work_dir = f'{TEST_DIR}/work/test_metadata_convert_options'
        setup_functions.make_work_dir(cls.work_dir, clear_dir=True)


    def test_prefix_option(self):
        prefix = 'prefix'
        command = [f'--prefix={prefix}', f'{self.metadata_dir}/HeLa_metadata.json']
        result = setup_functions.run_main(
            metadata_convert._main, command, self.work_dir, prog=PROG
        )
        self.assertEqual(0, result.returncode, result.stderr)

        self.assertTrue(os.path.isfile(f'{self.work_dir}/{prefix}.annotations.csv'))
        self.assertTrue(os.path.isfile(f'{self.work_dir}/{prefix}.definitions.bat'))


class TestToFileBase(setup_functions.AbstractTestsBase):

    META_TYPES = {'string_var': Dtype.STRING,
                  'bool_var': Dtype.BOOL,
                  'int_var': Dtype.INT,
                  'float_var': Dtype.FLOAT,
                  'na_var': Dtype.NULL}


    def __init__(self):
        self.work_dir = None
        self.metadata_file = None
        self.result = None
        self.out_ext = None


    def test_is_sucessful(self):
        self.assertEqual(0, self.result.returncode, self.result.stderr)


    @mock.patch('DIA_QC_report.submodules.read_metadata.LOGGER', mock.Mock())
    def test_annotations(self):
        metadata = f'{self.work_dir}/{os.path.splitext(os.path.basename(self.metadata_file))[0]}.{self.out_ext}'
        self.assertTrue(os.path.isfile(metadata))

        writer = Metadata()
        self.assertTrue(writer.read(self.metadata_file, exclude_null_from_skyline=True))
        types_in = writer.types

        reader = Metadata()
        self.assertTrue(reader.read(metadata, exclude_null_from_skyline=False))
        types_out = reader.types

        # check that data types are correct
        self.assertDictEqual(self.META_TYPES, types_in)
        self.assertDictEqual(self.META_TYPES, types_out)

        # check that data dict is equal
        self.assertEqual(writer, reader)


class TestToCsv(unittest.TestCase, TestToFileBase):
    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{TEST_DIR}/work/test_metadata_to_csv'
        setup_functions.make_work_dir(cls.work_dir, clear_dir=True)
        cls.metadata_file = f'{TEST_DIR}/data/metadata/Strap_multi_var_metadata.tsv'
        cls.out_ext = 'csv'

        commands = ['--out=csv', cls.metadata_file]
        cls.result = setup_functions.run_main(
            metadata_convert._main, commands, cls.work_dir, prefix=__name__, prog=PROG
        )


class TestToTsv(unittest.TestCase, TestToFileBase):
    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{TEST_DIR}/work/test_metadata_to_tsv'
        setup_functions.make_work_dir(cls.work_dir, clear_dir=True)
        cls.metadata_file = f'{TEST_DIR}/data/metadata/Strap_missing_multi_var_metadata.json'
        cls.out_ext = 'tsv'

        commands = ['--out=tsv', cls.metadata_file]
        cls.result = setup_functions.run_main(
            metadata_convert._main, commands, cls.work_dir, prefix=__name__, prog=PROG
        )


class TestToJson(unittest.TestCase, TestToFileBase):
    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{TEST_DIR}/work/test_metadata_to_json'
        setup_functions.make_work_dir(cls.work_dir, clear_dir=True)
        cls.metadata_file = f'{TEST_DIR}/data/metadata/Strap_multi_var_metadata.tsv'
        cls.out_ext = 'json'

        commands = ['--out=json', cls.metadata_file]
        cls.result = setup_functions.run_main(
            metadata_convert._main, commands, cls.work_dir, prefix=__name__, prog=PROG
        )


class TestToSkylineBase(TestToFileBase):

    SKY_TYPES = {'string_var': 'text',
                 'bool_var': 'true_false',
                 'int_var': 'number',
                 'float_var': 'number',
                 'na_var': 'text'}

    def test_batch_file(self):
        annotation_bat = f'{self.work_dir}/{os.path.splitext(os.path.basename(self.metadata_file))[0]}.definitions.bat'
        self.assertTrue(os.path.isfile(annotation_bat))

        command_re = re.compile(r'^--annotation-name="([\w\.\-]+)" --annotation-targets=replicate --annotation-type=([a-z_]+)$')

        with open(annotation_bat, 'r') as inF:
            lines = inF.readlines()

        for line in lines:
            self.assertTrue((match := command_re.search(line)) is not None)
            self.assertTrue(match[1] in self.SKY_TYPES)
            self.assertEqual(self.SKY_TYPES[match[1]], match[2])


class TestTsvToSkylineCsv(unittest.TestCase, TestToSkylineBase):
    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{TEST_DIR}/work/test_tsv_to_skyline_csv'
        setup_functions.make_work_dir(cls.work_dir, clear_dir=True)
        cls.metadata_file = f'{TEST_DIR}/data/metadata/Strap_multi_var_metadata.tsv'
        cls.out_ext = 'annotations.csv'

        commands = ['--out=skyline', cls.metadata_file]
        cls.result = setup_functions.run_main(
            metadata_convert._main, commands, cls.work_dir, prefix=__name__, prog=PROG
        )


class TestCsvToSkylineCsv(unittest.TestCase, TestToSkylineBase):
    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{TEST_DIR}/work/test_csv_to_skyline_csv'
        setup_functions.make_work_dir(cls.work_dir, clear_dir=True)
        cls.metadata_file = f'{TEST_DIR}/data/metadata/HeLa_metadata.csv'
        cls.out_ext = 'annotations.csv'

        commands = ['--out=skyline', cls.metadata_file]
        cls.result = setup_functions.run_main(
            metadata_convert._main, commands, cls.work_dir, prefix=__name__, prog=PROG
        )


class TestJsonToSkylineCsv(unittest.TestCase, TestToSkylineBase):
    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{TEST_DIR}/work/test_json_to_skyline_csv'
        setup_functions.make_work_dir(cls.work_dir, clear_dir=True)
        cls.metadata_file = f'{TEST_DIR}/data/metadata/HeLa_metadata.json'
        cls.out_ext = 'annotations.csv'

        commands = ['--out=skyline', cls.metadata_file]
        cls.result = setup_functions.run_main(
            metadata_convert._main, commands, cls.work_dir, prefix=__name__, prog=PROG
        )


class TestProblamaticHeaders(unittest.TestCase, TestToSkylineBase):
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
        cls.out_ext = 'annotations.csv'

        commands = ['--out=skyline', cls.metadata_file]
        cls.result = setup_functions.run_main(
            metadata_convert._main, commands, cls.work_dir, prefix=__name__, prog=PROG
        )
