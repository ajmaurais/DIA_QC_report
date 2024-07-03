
import unittest
import os
import sqlite3

import pandas as pd
from numpy import log2

import setup_functions

import DIA_QC_report.submodules.dia_db_utils as db_utils


class CommonTests(setup_functions.AbstractTestsBase):
    def __init__(self):
        self.conn = None
        self.precursor_method = None
        self.protein_method = None


    @classmethod
    def tearDownClass(cls):
        if cls.conn is not None:
            cls.conn.close()


    def test_is_normalized(self):
        self.assertIsNotNone(self.conn)
        self.assertTrue(db_utils.is_normalized(self.conn))
        self.assertEqual(db_utils.get_meta_value(self.conn, 'precursor_normalization_method'),
                         self.precursor_method)
        self.assertEqual(db_utils.get_meta_value(self.conn, 'protein_normalization_method'),
                         self.protein_method)


    def test_precursor_medians_equal(self):
        self.assertIsNotNone(self.conn)

        query = ''' SELECT
                r.replicate,
                p.modifiedSequence,
                p.precursorCharge,
                p.normalizedArea
            FROM precursors p
            LEFT JOIN replicates r ON r.id == p.replicateId
            WHERE p.normalizedArea IS NOT NULL; '''
        df = pd.read_sql(query, self.conn)

        self.assertTrue(not any(df['normalizedArea'] == 0))

        df['log2NormArea'] = log2(df['normalizedArea'] + 1)
        medians = df.groupby('replicate')['log2NormArea'].median()

        self.assertTrue(all(medians.iloc[0] == medians))


class TestSingleProject(CommonTests):
    TEST_PROJECT = 'Sp3'

    @classmethod
    def run_commands(cls, normalize_command):
        cls.parse_result = setup_functions.setup_single_db(cls.data_dir,
                                                           cls.work_dir,
                                                           cls.TEST_PROJECT,
                                                           clear_dir=True)

        if cls.parse_result.returncode != 0:
            raise RuntimeError('Setup of test db failed!')

        cls.normalize_result = setup_functions.run_command(normalize_command,
                                                           cls.work_dir,
                                                           prefix='normalize_single_proj')

        cls.conn = None
        if cls.normalize_result.returncode == 0:
            if os.path.isfile(cls.db_path):
                cls.conn = sqlite3.connect(cls.db_path)


    def test_is_successful(self):
        self.assertEqual(self.parse_result.returncode, 0)
        self.assertEqual(self.normalize_result.returncode, 0)


class TestMedianSingle(unittest.TestCase, TestSingleProject):
    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_normalize_db_median_single_project/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        cls.precursor_method = 'median'
        cls.protein_method = 'median'

        normalize_command = ['normalize_db', '-m=median', cls.db_path]
        cls.run_commands(normalize_command)


    def test_protein_medians_equal(self):
        self.assertIsNotNone(self.conn)

        query = ''' SELECT
                r.replicate,
                q.proteinId,
                q.normalizedAbundance
            FROM proteinQuants q
            LEFT JOIN replicates r ON r.id == q.replicateId
            WHERE q.normalizedAbundance IS NOT NULL; '''
        df = pd.read_sql(query, self.conn)

        self.assertTrue(not any(df['normalizedAbundance'] == 0))

        df['log2NormAbundance'] = log2(df['normalizedAbundance'] + 1)
        medians = df.groupby('replicate')['log2NormAbundance'].median()

        self.assertTrue(all(medians.iloc[0] == medians))


class TestDirectLFQSingle(unittest.TestCase, TestSingleProject):
    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_normalize_db_DirectLFQ_single_project/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        cls.precursor_method = 'median'
        cls.protein_method = 'DirectLFQ'

        normalize_command = ['normalize_db', '-m=DirectLFQ', cls.db_path]
        cls.run_commands(normalize_command)


    def test_keep_missing_fails(self):
        command = ['normalize_db', '-m=DirectLFQ', '--keepMissing', self.db_path]

        result = setup_functions.run_command(command, self.work_dir,
                                             prefix='test_keepMissing')
        self.assertEqual(result.returncode, 1)
        self.assertTrue('--keepMissing option not compatible with DirectLFQ' in result.stderr)


class TestMultiProject(CommonTests):
    @classmethod
    def run_commands(cls, normalize_command):
        if any(result.returncode != 0 for result in cls.parse_results):
            raise RuntimeError('Setup of test db failed!')

        cls.normalize_result = setup_functions.run_command(normalize_command,
                                                           cls.work_dir,
                                                           prefix='normalize_multi_proj')

        cls.conn = None
        if cls.normalize_result.returncode == 0:
            if os.path.isfile(cls.db_path):
                cls.conn = sqlite3.connect(cls.db_path)


    def test_is_successful(self):
        self.assertTrue(all(result.returncode == 0 for result in self.parse_results))
        self.assertEqual(self.normalize_result.returncode, 0)


    def test_no_new_protein_quants(self):
        '''
        Check that there aren't any protein quants that have NA unnormalized abundances,
        but finite normalized abundances. This should never happen. If it does there is
        something wrong with the peptideToProtein table.
        '''
        self.assertIsNotNone(self.conn)

        cur = self.conn.cursor()
        cur.execute('''
            SELECT
                replicateId, proteinId,
                abundance, normalizedAbundance
            FROM proteinQuants
            WHERE abundance IS NULL and normalizedAbundance is NOT NULL;''')
        self.assertEqual(0, len(cur.fetchall()))


    def test_precursor_db_update_sucessful(self):
        self.assertTrue(self.conn is not None)

        keep_cols = ['replicateId', 'ion', 'normalizedArea']

        def df_to_dict(df):
            ret = dict()
            for row in df.itertuples():
                key = f'{row.replicateId}_{row.ion}'
                ret[key] = row.normalizedArea

            return ret

        query = ''' SELECT
                replicateId,
                peptideId,
                modifiedSequence,
                precursorCharge,
                normalizedArea
            FROM precursors; '''
        db_df = pd.read_sql(query, self.conn)
        db_df['ion'] = db_df.apply(lambda x: f"{x['modifiedSequence']}_{x['precursorCharge']}", axis=1)
        db_df = db_df[keep_cols]

        gt_df = pd.read_csv(f'{self.data_dir}/intermediate_files/multi_project_ion_df.tsv', sep='\t')
        gt_df = gt_df[keep_cols]

        df_norm_area = df_to_dict(gt_df)
        db_norm_area = df_to_dict(db_df)

        df_keys = set(df_norm_area.keys())
        db_keys = set(db_norm_area.keys())

        self.assertTrue(df_keys <= db_keys)

        for ion in df_norm_area:
            if pd.isna(df_norm_area[ion]):
                self.assertTrue(pd.isna(db_norm_area[ion]))
            else:
                self.assertAlmostEqual(df_norm_area[ion], db_norm_area[ion], places=5)


    def test_protein_db_update_sucessful(self):
        self.assertTrue(self.conn is not None)

        keep_cols = ['replicateId', 'proteinId', 'normalizedAbundance']

        def df_to_dict(df):
            ret = dict()
            for row in df.itertuples():
                key = f'{row.replicateId}_{row.proteinId}'
                ret[key] = row.normalizedAbundance

            return ret

        query = ''' SELECT
                replicateId,
                proteinId,
                normalizedAbundance
            FROM proteinQuants; '''
        db_df = pd.read_sql(query, self.conn)
        db_df = db_df[keep_cols]

        gt_df = pd.read_csv(f'{self.data_dir}/intermediate_files/multi_project_{self.protein_method}_protein_df.tsv', sep='\t')
        gt_df = gt_df[keep_cols]

        df_norm_area = df_to_dict(gt_df)
        db_norm_area = df_to_dict(db_df)

        df_keys = set(df_norm_area.keys())
        db_keys = set(db_norm_area.keys())

        self.assertTrue(df_keys <= db_keys)

        for protein in df_norm_area:
            if pd.isna(df_norm_area[protein]):
                self.assertTrue(pd.isna(db_norm_area[protein]))
            else:
                self.assertAlmostEqual(df_norm_area[protein], db_norm_area[protein], places=5)


class TestMedianMulti(unittest.TestCase, TestMultiProject):
    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_normalize_db_median_multi_project/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'
        cls.parse_results = setup_functions.setup_multi_db(cls.data_dir, cls.work_dir)

        cls.precursor_method = 'median'
        cls.protein_method = 'median'

        normalize_command = ['normalize_db', '-m=median', cls.db_path]
        cls.run_commands(normalize_command)


class TestDirectLFQMulti(unittest.TestCase, TestMultiProject):
    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_normalize_db_DirectLFQ_multi_project/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'
        cls.parse_results = setup_functions.setup_multi_db(cls.data_dir, cls.work_dir)

        cls.precursor_method = 'median'
        cls.protein_method = 'DirectLFQ'

        normalize_command = ['normalize_db', '-m=DirectLFQ', cls.db_path]
        cls.run_commands(normalize_command)


    def test_keep_missing_fails(self):
        command = ['normalize_db', '-m=DirectLFQ', '--keepMissing', self.db_path]

        result = setup_functions.run_command(command, self.work_dir,
                                             prefix='test_keepMissing')
        self.assertEqual(result.returncode, 1)
        self.assertTrue('--keepMissing option not compatible with DirectLFQ' in result.stderr)


class TestAllPrecursorsMissing(unittest.TestCase):
    TEST_PROJECT = 'GFP'

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_normalize_db_missing_precursors/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        command = ['parse_data', f'--projectName={cls.TEST_PROJECT}',
                   f'{cls.data_dir}/skyline_reports/{cls.TEST_PROJECT}_replicate_quality.tsv',
                   f'{cls.data_dir}/skyline_reports/{cls.TEST_PROJECT}_precursor_quality.tsv']

        setup_functions.make_work_dir(cls.work_dir, True)
        cls.parse_result = setup_functions.run_command(command, cls.work_dir,
                                                       prefix='parse_missing_precursors')

        if cls.parse_result.returncode != 0:
            raise RuntimeError('Setup of test db failed!')

        cls.conn = None
        if os.path.isfile(cls.db_path):
            cls.conn = sqlite3.connect(cls.db_path)


    @classmethod
    def tearDownClass(cls):
        if cls.conn is not None:
            cls.conn.close()


    def test_median_normalization(self):
        self.assertIsNotNone(self.conn)

        norm_levels = ['precursor_normalization_method', 'protein_normalization_method']
        cur = self.conn.cursor()
        for level in norm_levels:
            cur.execute('DELETE FROM metadata WHERE key == ?', (level,))
        cur.execute("UPDATE metadata SET value == FALSE WHERE key == 'is_normalized'")
        self.conn.commit()

        self.assertFalse(db_utils.is_normalized(self.conn))
        with self.assertLogs(db_utils.LOGGER, level='ERROR') as cm:
            for level in norm_levels:
                self.assertIsNone(db_utils.get_meta_value(self.conn, level))

        normalize_command = ['normalize_db', '-m=median', self.db_path]
        normalize_result = setup_functions.run_command(normalize_command,
                                                       self.work_dir,
                                                       prefix='normalize_db_median')

        self.assertTrue(db_utils.is_normalized(self.conn))
        for level in norm_levels:
            self.assertEqual(db_utils.get_meta_value(self.conn, level), 'median')


    def test_DirectLFQ_normalization(self):
        normalize_command = ['normalize_db', '-m=DirectLFQ', self.db_path]

        normalize_result = setup_functions.run_command(normalize_command,
                                                       self.work_dir,
                                                       prefix='normalize_db_DirectLFQ')

        self.assertEqual(normalize_result.returncode, 1)
        self.assertTrue('ERROR: Can not perform DirectLFQ normalization with 0 precursors!' in normalize_result.stderr)
