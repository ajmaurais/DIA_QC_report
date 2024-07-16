
import unittest
from unittest import mock
import argparse
import sys
import os
import sqlite3
import re

import setup_functions

import DIA_QC_report.submodules.dia_db_utils as db_utils
import DIA_QC_report.generate_batch_rmd as generate_batch_rmd


class TestMainFxns(unittest.TestCase):
    def test_remove_bc_tables_bc_options(self):
        directions = ('wide', 'long')
        tests = [(directions, option) for option in  ('44', '55', '66', '77')]
        tests += [(('wide',), f'{option}0') for option in range(4, 8)]
        tests += [(('long',), f'0{option}') for option in range(4, 8)]

        for test_directions, option in tests:
            tables = db_utils.parse_bitmask_options(option, directions,
                                                    generate_batch_rmd.PYTHON_METHOD_NAMES)

            for direction in test_directions:
                self.assertTrue(tables[direction]['batch_corrected'])

            with self.assertLogs(generate_batch_rmd.LOGGER) as cm:
                tables = generate_batch_rmd.remove_bc_tables(tables)

            self.assertTrue(len(test_directions) == len(cm.output))

            for direction, out in zip(test_directions, cm.output):
                self.assertEqual(out, f'WARNING:root:Batch corrected {direction} table not available when batch correction is skipped!')


    def test_remove_bc_tables_no_bc_options(self):
        directions = ('wide', 'long')
        tests = [(directions, option) for option in  ('00', '11', '22', '33')]
        tests += [(('wide',), f'{option}0') for option in range(4)]
        tests += [(('long',), f'0{option}') for option in range(4)]

        for test_directions, option in tests:
            tables = db_utils.parse_bitmask_options(option, directions,
                                                    generate_batch_rmd.PYTHON_METHOD_NAMES)

            for direction in test_directions:
                self.assertFalse(tables[direction]['batch_corrected'])

            with self.assertNoLogs(generate_batch_rmd.LOGGER) as cm:
                tables = generate_batch_rmd.remove_bc_tables(tables)


class TestMakeBatchRmd(unittest.TestCase):
    RENDER_RMD = False

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_generate_batch_rmd/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        # remove tables subdirectory in work_dir if necissary
        if os.path.isdir(f'{cls.work_dir}/tables'):
            for file in os.listdir(f'{cls.work_dir}/tables'):
                os.remove(f'{cls.work_dir}/tables/{file}')
            os.rmdir(f'{cls.work_dir}/tables')

        cls.parse_results = setup_functions.setup_multi_db(cls.data_dir,
                                                           cls.work_dir,
                                                           clear_dir=True)

        if not all(result.returncode == 0 for result in cls.parse_results):
            raise RuntimeError('Setup of test db failed!')

        normalize_command = ['dia_qc', 'normalize', cls.db_path]
        cls.normalize_result = setup_functions.run_command(normalize_command,
                                                           cls.work_dir,
                                                           prefix='normalize_db')

        cls.conn = None
        if cls.normalize_result.returncode == 0:
            if os.path.isfile(cls.db_path):
                cls.conn = sqlite3.connect(cls.db_path)


    @classmethod
    def tearDownClass(cls):
        if cls.conn is not None:
            cls.conn.close()


    def test_is_successful(self):
        self.assertTrue(self.conn is not None)
        self.assertTrue(db_utils.is_normalized(self.conn))

        rmd_name = 'basic_test'
        command = ['dia_qc', 'batch_rmd',
                   '-o', f'{rmd_name}.rmd', self.db_path]
        result = setup_functions.run_command(command, self.work_dir)

        self.assertEqual(result.returncode, 0)
        self.assertTrue(os.path.isfile(f'{self.work_dir}/{rmd_name}.rmd'))

        if self.RENDER_RMD:
            render_command = ['Rscript', '-e', f"rmarkdown::render('{rmd_name}.rmd')"]
            render_result = setup_functions.run_command(render_command, self.work_dir,
                                                        prefix=f'render_{rmd_name}')
            self.assertEqual(render_result.returncode, 0)

            for file in [f'{rmd_name}.html',
                         'proteins_batch_corrected_wide.tsv',
                         'precursors_batch_corrected_wide.tsv',
                         'metadata_wide.tsv']:
                self.assertTrue(os.path.isfile(f'{self.work_dir}/{file}'))


    def test_is_normalized_check(self):
        self.assertTrue(self.conn is not None)
        self.assertTrue(db_utils.is_normalized(self.conn))

        try:
            # set is_normalized metadata entry to False
            self.conn = db_utils.update_meta_value(self.conn, 'is_normalized', 'False')

            # make sure generate_batch_rmd fails
            command = ['dia_qc', 'batch_rmd', self.db_path]
            result = setup_functions.run_command(command, self.work_dir)
            self.assertEqual(result.returncode, 1)
            self.assertTrue('Database file it not normalized!' in result.stderr)

        finally:
            # reset is_normalized metadata entry to True
            self.conn = db_utils.update_meta_value(self.conn, 'is_normalized', 'True')


    @mock.patch('DIA_QC_report.submodules.dia_db_utils.LOGGER', mock.Mock())
    def test_single_batch_works(self):
        self.assertTrue(self.conn is not None)

        try:
            # set 'Strap' project replicates skipped
            self.assertTrue(db_utils.mark_reps_skipped(self.conn, projects=('Strap',)))

            rmd_name = 'single_batch'
            command = ['dia_qc', 'batch_rmd',
                       '--proteinTables=00', '--precursorTables=00',
                       '-o', f'{rmd_name}.rmd', self.db_path]
            result = setup_functions.run_command(command, self.work_dir)

            self.assertEqual(result.returncode, 0)
            self.assertTrue(os.path.isfile(f'{self.work_dir}/{rmd_name}.rmd'))
            self.assertTrue('WARNING: Only 1 project in replicates! Skipping batch correction.' in result.stderr)

            for flag in ['precursorTables', 'proteinTables']:
                command = ['dia_qc', 'batch_rmd', f'--{flag}=44', self.db_path]
                result = setup_functions.run_command(command, self.work_dir,
                                                     prefix='invalid_table_options')
                self.assertEqual(result.returncode, 0)

                name = re.sub('^p', 'P', flag.replace('Tables', ''))
                for d in ('long', 'wide'):
                    self.assertTrue(f'WARNING: {name} batch corrected {d} table not available when batch correction is skipped!' in result.stderr)

            if self.RENDER_RMD:
                render_command = ['Rscript', '-e', f"rmarkdown::render('{rmd_name}.rmd')"]
                render_result = setup_functions.run_command(render_command, self.work_dir,
                                                            prefix=f'render_{rmd_name}')
                self.assertEqual(render_result.returncode, 0)

        finally:
            # Set all replicates to be included
            db_utils.mark_all_reps_includced(self.conn)


    def test_controlKey_check(self):
        self.assertTrue(self.conn is not None)

        command = ['dia_qc', 'batch_rmd', '--addControlValue=A549', self.db_path]
        result = setup_functions.run_command(command, self.work_dir)
        self.assertEqual(result.returncode, 1)
        self.assertTrue('No control key specified!' in result.stderr)

        # make sure --skipTests option works
        command.insert(2, '--skipTests')
        result = setup_functions.run_command(command, self.work_dir,
                                             prefix='test_controlKey_check_skipTests')
        self.assertEqual(result.returncode, 0)


    def test_addControlValue_check(self):
        self.assertTrue(self.conn is not None)

        command = ['dia_qc', 'batch_rmd', '--controlKey=cellLine', self.db_path]
        result = setup_functions.run_command(command, self.work_dir)
        self.assertEqual(result.returncode, 1)
        self.assertTrue('No control value(s) specified!' in result.stderr)

        # make sure --skipTests option works
        command.insert(2, '--skipTests')
        result = setup_functions.run_command(command, self.work_dir,
                                             prefix='test_addControlValue_check_skipTests')
        self.assertEqual(result.returncode, 0)


class TestMissingMetadata(unittest.TestCase):
    TEST_PROJECT = 'Strap'
    RENDER_RMD = False

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_batch_rmd_missing_metadata/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        # remove tables subdirectory in work_dir if necissary
        if os.path.isdir(f'{cls.work_dir}/tables'):
            for file in os.listdir(f'{cls.work_dir}/tables'):
                os.remove(f'{cls.work_dir}/tables/{file}')
            os.rmdir(f'{cls.work_dir}/tables')

        cls.parse_result = setup_functions.setup_single_db(cls.data_dir,
                                                           cls.work_dir,
                                                           cls.TEST_PROJECT,
                                                           metadata_suffix='_missing_multi_var_metadata.tsv',
                                                           clear_dir=True)

        if cls.parse_result.returncode != 0:
            raise RuntimeError('Setup of test db failed!')

        normalize_command = ['dia_qc', 'normalize', cls.db_path]
        cls.normalize_result = setup_functions.run_command(normalize_command,
                                                           cls.work_dir,
                                                           prefix='normalize_db')

        if cls.normalize_result.returncode != 0:
            raise RuntimeError('Setup of test db failed!')

        cls.conn = None
        if os.path.isfile(cls.db_path):
            cls.conn = sqlite3.connect(cls.db_path)


    @classmethod
    def tearDownClass(cls):
        if cls.conn is not None:
            cls.conn.close()


    def test_is_successful(self):
        self.assertTrue(self.conn is not None)
        self.assertTrue(db_utils.is_normalized(self.conn))

        rmd_name = 'basic_test'
        command = ['dia_qc', 'batch_rmd',
                   '-c=string_var', '-c=bool_var', '-c=int_var', '-c=float_var',
                   '-o', f'{rmd_name}.rmd', self.db_path]
        result = setup_functions.run_command(command, self.work_dir)

        self.assertEqual(result.returncode, 0)
        self.assertTrue(os.path.isfile(f'{self.work_dir}/{rmd_name}.rmd'))
        self.assertTrue('WARNING: Only 1 project in replicates! Skipping batch correction.' in result.stderr)

        if self.RENDER_RMD:
            render_command = ['Rscript', '-e', f"rmarkdown::render('{rmd_name}.rmd')"]
            render_result = setup_functions.run_command(render_command, self.work_dir,
                                                        prefix=f'render_{rmd_name}')
            self.assertEqual(render_result.returncode, 0)

            for file in [f'{rmd_name}.html', 'metadata_wide.tsv']:
                self.assertTrue(os.path.isfile(f'{self.work_dir}/{file}'))


class TestBadMetadataHeaders(unittest.TestCase):
    RENDER_RMD = False

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_batch_rmd_bad_metadata_headers/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        # remove tables subdirectory in work_dir if necissary
        if os.path.isdir(f'{cls.work_dir}/tables'):
            for file in os.listdir(f'{cls.work_dir}/tables'):
                os.remove(f'{cls.work_dir}/tables/{file}')
            os.rmdir(f'{cls.work_dir}/tables')

        cls.parse_results = setup_functions.setup_multi_db(cls.data_dir,
                                                           cls.work_dir,
                                                           clear_dir=True,
                                                           metadata_suffix='_bad_headers_metadata.tsv')

        if not all(result.returncode == 0 for result in cls.parse_results):
            raise RuntimeError('Setup of test db failed!')

        normalize_command = ['dia_qc', 'normalize', cls.db_path]
        cls.normalize_result = setup_functions.run_command(normalize_command,
                                                           cls.work_dir,
                                                           prefix='normalize_db')

        cls.conn = None
        if cls.normalize_result.returncode == 0:
            if os.path.isfile(cls.db_path):
                cls.conn = sqlite3.connect(cls.db_path)


    @classmethod
    def tearDownClass(cls):
        if cls.conn is not None:
            cls.conn.close()


    def test_space_headers(self):
        self.assertTrue(self.conn is not None)
        self.assertTrue(db_utils.is_normalized(self.conn))

        rmd_name = 'space_header_test'
        command = ['dia_qc', 'batch_rmd',
                   '-c', 'string var', '-c', 'bool var', '-c', 'int var', '-c', 'float var',
                   '--controlKey', 'string var',
                   '--addControlValue', 'NCI7 7-pool', '--addControlValue', 'NCI7 4-pool',
                   '-o', f'{rmd_name}.rmd', self.db_path]
        result = setup_functions.run_command(command, self.work_dir)

        self.assertEqual(result.returncode, 0)
        self.assertTrue(os.path.isfile(f'{self.work_dir}/{rmd_name}.rmd'))

        if self.RENDER_RMD:
            render_command = ['Rscript', '-e', f"rmarkdown::render('{rmd_name}.rmd')"]
            render_result = setup_functions.run_command(render_command, self.work_dir,
                                                        prefix=f'render_{rmd_name}')
            self.assertEqual(render_result.returncode, 0)

            for file in [f'{rmd_name}.html',
                         'proteins_batch_corrected_wide.tsv',
                         'precursors_batch_corrected_wide.tsv',
                         'metadata_wide.tsv']:
                self.assertTrue(os.path.isfile(f'{self.work_dir}/{file}'))


    def test_symbol_headers(self):
        self.assertTrue(self.conn is not None)
        self.assertTrue(db_utils.is_normalized(self.conn))

        rmd_name = 'symbol_header_test'
        command = ['dia_qc', 'batch_rmd',
                   '-c', 'string var', '-c', 'This_is-a@bad~header$ (alphanumeric ONLY please!)',
                   '--controlKey', 'This_is-a@bad~header$ (alphanumeric ONLY please!)',
                   '--addControlValue', 'NCI7 7-pool', '--addControlValue', 'NCI7 4-pool',
                   '-o', f'{rmd_name}.rmd', self.db_path]
        result = setup_functions.run_command(command, self.work_dir)

        self.assertEqual(result.returncode, 0)
        self.assertTrue(os.path.isfile(f'{self.work_dir}/{rmd_name}.rmd'))

        if self.RENDER_RMD:
            render_command = ['Rscript', '-e', f"rmarkdown::render('{rmd_name}.rmd')"]
            render_result = setup_functions.run_command(render_command, self.work_dir,
                                                        prefix=f'render_{rmd_name}')
            self.assertEqual(render_result.returncode, 0)

            for file in [f'{rmd_name}.html',
                         'proteins_batch_corrected_wide.tsv',
                         'precursors_batch_corrected_wide.tsv',
                         'metadata_wide.tsv']:
                self.assertTrue(os.path.isfile(f'{self.work_dir}/{file}'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Tests for generate_batch_rmd')
    parser.add_argument('-r', '--render', action='store_true', default=False,
                        help='Also test if rmd file can be rendered?')
    args, unittest_args = parser.parse_known_args()

    TestMakeBatchRmd.RENDER_RMD = args.render
    TestMissingMetadata.RENDER_RMD = args.render
    TestBadMetadataHeaders.RENDER_RMD = args.render

    unittest_args.insert(0, __file__)
    unittest.main(argv=unittest_args, verbosity=2)
