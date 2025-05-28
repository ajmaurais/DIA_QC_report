
import unittest
import os
import sqlite3
import re

import setup_functions

import DIA_QC_report.submodules.dia_db_utils as db_utils
from DIA_QC_report import generate_batch_rmd
from DIA_QC_report.submodules.normalization import NORMALIZATION_METHODS


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


class TestNormMethodsKnown(unittest.TestCase):
    def test_all_norm_methods_known(self):
        self.assertEqual(len(generate_batch_rmd.NORMALIZATION_METHOD_NAMES),
                         len(NORMALIZATION_METHODS))

        for method in NORMALIZATION_METHODS:
            self.assertTrue(method in generate_batch_rmd.NORMALIZATION_METHOD_NAMES)


class TestMakeBatchRmd(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.render_rmd = os.getenv('RENDER_RMD', 'False').lower() == 'true'

        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_generate_batch_rmd/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        # remove html artifact subdirectory in work_dir if necissary
        for d in ('basic_test_files', 'single_batch_files'):
            setup_functions.remove_test_dir(f'{cls.work_dir}/{d}', recursive=True)

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

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(os.path.isfile(f'{self.work_dir}/{rmd_name}.rmd'))

        if self.render_rmd:
            render_command = ['Rscript', '-e', f"rmarkdown::render('{rmd_name}.rmd')"]
            render_result = setup_functions.run_command(render_command, self.work_dir,
                                                        prefix=f'render_{rmd_name}')
            self.assertEqual(render_result.returncode, 0)

            self.assertFalse(os.path.isdir(f'{self.work_dir}/plots'))
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
            db_utils.update_meta_value(self.conn, db_utils.IS_NORMALIZED, 'False')

            # make sure generate_batch_rmd fails
            command = ['dia_qc', 'batch_rmd', self.db_path]
            result = setup_functions.run_command(command, self.work_dir)
            self.assertEqual(result.returncode, 1)
            self.assertTrue('Database file it not normalized!' in result.stderr)

        finally:
            # reset is_normalized metadata entry to True
            db_utils.update_meta_value(self.conn, db_utils.IS_NORMALIZED, 'True')


    def test_norm_method_found_unknown_value(self):
        self.assertTrue(self.conn is not None)

        current_method = db_utils.get_meta_value(self.conn, db_utils.PROTEIN_NORM_METHOD)
        try:
            # set is_normalized metadata entry to False
            db_utils.update_meta_value(self.conn, db_utils.PROTEIN_NORM_METHOD, 'Nothing')

            # make sure generate_batch_rmd fails
            command = ['dia_qc', 'batch_rmd', self.db_path]
            result = setup_functions.run_command(command, self.work_dir)
            self.assertEqual(result.returncode, 1)
            self.assertTrue("Unknown normalization method: 'Nothing'" in result.stderr)

        finally:
            # reset protein_normalization_method to original method
            db_utils.update_meta_value(self.conn, db_utils.PROTEIN_NORM_METHOD, current_method)


    def test_norm_method_found_missing_key(self):
        self.assertTrue(self.conn is not None)

        cur = self.conn.cursor()
        cur.execute("SELECT * FROM metadata WHERE key == ?;", (db_utils.PRECURSOR_NORM_METHOD,))
        current_entry = cur.fetchall()
        self.assertEqual(len(current_entry), 1)

        try:
            # delete current entry
            cur.execute("DELETE FROM metadata WHERE key == ?", (db_utils.PRECURSOR_NORM_METHOD,))
            self.conn.commit()

            # make sure generate_batch_rmd fails
            command = ['dia_qc', 'batch_rmd', self.db_path]
            result = setup_functions.run_command(command, self.work_dir)
            self.assertEqual(result.returncode, 1)
            self.assertTrue('Missing normalization methods in metadata table!' in result.stderr)

        finally:
            # reset original precursor_normalization_method entry
            cur = self.conn.cursor()
            cur.executemany('INSERT INTO metadata (key, value) VALUES (?, ?);', current_entry)
            self.conn.commit()


    def test_single_batch_works(self):
        self.assertTrue(self.conn is not None)

        try:
            # set 'Strap' project replicates skipped
            self.assertTrue(db_utils.mark_reps_skipped(self.conn, projects=('Strap',), quiet=True))

            rmd_name = 'single_batch'
            command = ['dia_qc', 'batch_rmd',
                       '--proteinTables=00', '--precursorTables=00', '--metadataTables=00',
                       '-o', f'{rmd_name}.rmd', self.db_path]
            result = setup_functions.run_command(command, self.work_dir)

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(os.path.isfile(f'{self.work_dir}/{rmd_name}.rmd'))
            self.assertTrue('WARNING: Only 1 project in replicates! Skipping batch correction.' in result.stderr)

            for flag in ['precursorTables', 'proteinTables']:
                command = ['dia_qc', 'batch_rmd', f'--{flag}=44', self.db_path]
                result = setup_functions.run_command(command, self.work_dir,
                                                     prefix='invalid_table_options')
                self.assertEqual(result.returncode, 0, result.stderr)

                name = re.sub('^p', 'P', flag.replace('Tables', ''))
                for d in ('long', 'wide'):
                    self.assertTrue(f'WARNING: {name} batch corrected {d} table not available when batch correction is skipped!' in result.stderr)

            if self.render_rmd:
                render_command = ['Rscript', '-e', f"rmarkdown::render('{rmd_name}.rmd')"]
                render_result = setup_functions.run_command(render_command, self.work_dir,
                                                            prefix=f'render_{rmd_name}')
                self.assertEqual(render_result.returncode, 0)

        finally:
            # Set all replicates to be included
            db_utils.mark_all_reps_included(self.conn, quiet=True)


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
        self.assertEqual(result.returncode, 0, result.stderr)


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
        self.assertEqual(result.returncode, 0, result.stderr)


class TestPDFReport(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.render_rmd = os.getenv('RENDER_RMD', 'False').lower() == 'true'

        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_generate_batch_rmd_pdf/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        setup_functions.remove_test_dir(f'{cls.work_dir}/pdf_test_files', recursive=True)
        cls.parse_results = setup_functions.setup_multi_db(cls.data_dir,
                                                           cls.work_dir,
                                                           clear_dir=True)

        if not all(result.returncode == 0 for result in cls.parse_results):
            raise RuntimeError('Setup of test db failed!')

        normalize_command = ['dia_qc', 'normalize', cls.db_path]
        cls.normalize_result = setup_functions.run_command(normalize_command,
                                                           cls.work_dir,
                                                           prefix='normalize_db')

        if cls.normalize_result.returncode != 0:
            raise RuntimeError('Setup of test db failed!')


    def test_pdf_report(self):
        self.assertTrue(self.normalize_result.returncode == 0)
        self.assertTrue(os.path.isfile(self.db_path))

        rmd_name = 'pdf_test'
        command = ['dia_qc', 'batch_rmd', '--savePlots=pdf',
                   '--proteinTables=77', '--precursorTables=77', '--metadataTables=11',
                   '-o', f'{rmd_name}.rmd', self.db_path]
        result = setup_functions.run_command(command, self.work_dir)

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(os.path.isfile(f'{self.work_dir}/{rmd_name}.rmd'))

        if self.render_rmd:
            render_command = ['Rscript', '-e',
                              f"rmarkdown::render('{rmd_name}.rmd', output_format='pdf_document', params=list(save_plots=FALSE, write_tables=FALSE))"]
            render_result = setup_functions.run_command(render_command, self.work_dir,
                                                        prefix=f'render_{rmd_name}')
            self.assertEqual(render_result.returncode, 0, render_result.stderr)
            self.assertTrue(os.path.isfile(f'{self.work_dir}/{rmd_name}.pdf'))
            self.assertFalse(os.path.isdir(f'{self.work_dir}/plots'))

            # Check that no .tsv files are in the work_dir
            for file in os.listdir(self.work_dir):
                self.assertFalse(file.endswith('.tsv'), f'Found unexpected .tsv file: {file}')


class TestInteractive(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.render_rmd = os.getenv('RENDER_RMD', 'False').lower() == 'true'

        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_generate_batch_rmd_interactive/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        # remove plots subdirectory in work_dir if necissary
        for d in ('plots', 'savePlots_files'):
            setup_functions.remove_test_dir(f'{cls.work_dir}/{d}', recursive=True)

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


    def test_savePlots(self):
        self.assertTrue(self.conn is not None)
        self.assertTrue(db_utils.is_normalized(self.conn))

        rmd_name = 'savePlots'
        command = ['dia_qc', 'batch_rmd', '--savePlots=pdf',
                   '--proteinTables=00', '--precursorTables=00', '--metadataTables=00',
                   '-o', f'{rmd_name}.rmd', self.db_path]
        result = setup_functions.run_command(command, self.work_dir)

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(os.path.isfile(f'{self.work_dir}/{rmd_name}.rmd'))

        if self.render_rmd:
            render_command = ['Rscript', '-e', f"rmarkdown::render('{rmd_name}.rmd')"]
            render_result = setup_functions.run_command(render_command, self.work_dir,
                                                        prefix=f'render_{rmd_name}')
            self.assertEqual(render_result.returncode, 0, render_result.stderr)

            self.assertTrue(os.path.isdir(f'{self.work_dir}/plots'))
            for file in [f'{rmd_name}.html',
                         'plots/cv_dist.pdf',
                         'plots/precursor_normalization.tiff',
                         'plots/precursor_pca.pdf',
                         'plots/protein_pca.pdf']:
                self.assertTrue(os.path.isfile(f'{self.work_dir}/{file}'))


    def area_dist_interactive(self, rmd_path):
        self.assertTrue(os.path.isfile(rmd_path))
        with open(rmd_path, 'r') as inF:
            text = inF.read()

        not_i, _ = generate_batch_rmd.precursor_norm_plot(None, interactive=False)
        if not_i in text:
            self.assertFalse('ggiraph::girafe' in not_i)
            return False

        i, _ = generate_batch_rmd.precursor_norm_plot(None, interactive=True)
        if i in text:
            self.assertTrue('ggiraph::girafe' in i)
            return True

        return None


    def pca_interactive(self, rmd_path):
        self.assertTrue(os.path.isfile(rmd_path))
        with open(rmd_path, 'r') as inF:
            text = inF.read()

        not_i_prec, _ = generate_batch_rmd.pca_plot('precursor', interactive=False)
        not_i_prot, _ = generate_batch_rmd.pca_plot('protein', interactive=False)
        if not_i_prec in text and not_i_prot in text:
            self.assertFalse('ggiraph::girafe' in not_i_prec)
            self.assertFalse('ggiraph::girafe' in not_i_prot)
            return False

        i_prec, _ = generate_batch_rmd.pca_plot('precursor', interactive=True)
        i_prot, _ = generate_batch_rmd.pca_plot('protein', interactive=True)
        if i_prec in text and i_prot in text:
            self.assertTrue('ggiraph::girafe' in i_prec)
            self.assertTrue('ggiraph::girafe' in i_prot)
            return True

        return None


    def test_interactive_options(self):
        self.assertTrue(self.conn is not None)
        self.assertTrue(db_utils.is_normalized(self.conn))

        for option in ('0', '1', '2', '3'):
            rmd_name = f'test_interactive_{option}.rmd'
            command = ['dia_qc', 'batch_rmd', f'--interactive={option}',
                       '--proteinTables=00', '--precursorTables=00', '--metadataTables=00',
                       '-o', rmd_name, self.db_path]
            result = setup_functions.run_command(command, self.work_dir)
            self.assertEqual(result.returncode, 0, result.stderr)

            plots = db_utils.parse_bitmask_options(option, ('plots',), ('area_dist', 'pca'))

            self.assertEqual(self.area_dist_interactive(f'{self.work_dir}{rmd_name}'),
                             plots['plots']['area_dist'])
            self.assertEqual(self.pca_interactive(f'{self.work_dir}/{rmd_name}'),
                             plots['plots']['pca'])


class TestMissingMetadata(unittest.TestCase):
    TEST_PROJECT = 'Strap'

    @classmethod
    def setUpClass(cls):
        cls.render_rmd = os.getenv('RENDER_RMD', 'False').lower() == 'true'

        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_batch_rmd_missing_metadata/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        # remove subdirectoies in work_dir if necissary
        for d in ('tables', 'basic_test_files'):
            setup_functions.remove_test_dir(f'{cls.work_dir}/{d}', recursive=True)

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

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(os.path.isfile(f'{self.work_dir}/{rmd_name}.rmd'))
        self.assertTrue('WARNING: Only 1 project in replicates! Skipping batch correction.' in result.stderr)

        if self.render_rmd:
            render_command = ['Rscript', '-e', f"rmarkdown::render('{rmd_name}.rmd')"]
            render_result = setup_functions.run_command(render_command, self.work_dir,
                                                        prefix=f'render_{rmd_name}')
            self.assertEqual(render_result.returncode, 0, render_result.stderr)

            for file in [f'{rmd_name}.html', 'metadata_wide.tsv']:
                self.assertTrue(os.path.isfile(f'{self.work_dir}/{file}'))


class TestBadMetadataHeaders(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.render_rmd = os.getenv('RENDER_RMD', 'False').lower() == 'true'

        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_batch_rmd_bad_metadata_headers/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

		# remove tables subdirectory in work_dir if necissary
        setup_functions.remove_test_dir(f'{cls.work_dir}/tables')

        # remove html artifact subdirectory in work_dir if necissary
        for d in ('space_header_test_files', 'symbol_header_test_files'):
            setup_functions.remove_test_dir(f'{cls.work_dir}/{d}', recursive=True)

        setup_functions.remove_test_dir(cls.work_dir, recursive=True)
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

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(os.path.isfile(f'{self.work_dir}/{rmd_name}.rmd'))

        if self.render_rmd:
            render_command = ['Rscript', '-e', f"rmarkdown::render('{rmd_name}.rmd')"]
            render_result = setup_functions.run_command(render_command, self.work_dir,
                                                        prefix=f'render_{rmd_name}')
            self.assertEqual(render_result.returncode, 0, render_result.stderr)

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

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(os.path.isfile(f'{self.work_dir}/{rmd_name}.rmd'))

        if self.render_rmd:
            render_command = ['Rscript', '-e', f"rmarkdown::render('{rmd_name}.rmd')"]
            render_result = setup_functions.run_command(render_command, self.work_dir,
                                                        prefix=f'render_{rmd_name}')
            self.assertEqual(render_result.returncode, 0, render_result.stderr)

            for file in [f'{rmd_name}.html',
                         'proteins_batch_corrected_wide.tsv',
                         'precursors_batch_corrected_wide.tsv',
                         'metadata_wide.tsv']:
                self.assertTrue(os.path.isfile(f'{self.work_dir}/{file}'))
