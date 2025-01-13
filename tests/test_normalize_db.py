
import unittest
import os
import sqlite3
from math import isclose

import pandas as pd
from numpy import log2

import setup_functions

import DIA_QC_report.submodules.dia_db_utils as db_utils


def precursor_medians_equal(conn, epislon=1e-6):
    query = ''' SELECT
            r.replicate,
            p.modifiedSequence,
            p.precursorCharge,
            p.normalizedArea
        FROM precursors p
        LEFT JOIN replicates r ON r.id == p.replicateId
        WHERE p.normalizedArea IS NOT NULL; '''
    df = pd.read_sql(query, conn)

    df['log2NormArea'] = log2(df['normalizedArea'] + 1)
    medians = df.groupby('replicate')['log2NormArea'].median()

    return max(abs(medians.iloc[0] - median) for median in medians) <= epislon


def protein_medians_eqeual(conn, epislon=1e-6):
    query = ''' SELECT
            r.replicate,
            q.proteinId,
            q.normalizedAbundance
        FROM proteinQuants q
        LEFT JOIN replicates r ON r.id == q.replicateId
        WHERE q.normalizedAbundance IS NOT NULL; '''
    df = pd.read_sql(query, conn)

    df['log2NormAbundance'] = log2(df['normalizedAbundance'] + 1)
    medians = df.groupby('replicate')['log2NormAbundance'].median()

    return max(abs(medians.iloc[0] - median) for median in medians) <= epislon


class CommonTests(setup_functions.AbstractTestsBase):
    def __init__(self):
        self.conn = None
        self.precursor_method = None
        self.protein_method = None
        self.command = None


    @classmethod
    def tearDownClass(cls):
        if cls.conn is not None:
            cls.conn.close()


    def test_is_normalized(self):
        self.assertIsNotNone(self.conn)
        self.assertTrue(db_utils.is_normalized(self.conn))
        self.assertEqual(db_utils.get_meta_value(self.conn, db_utils.PRECURSOR_NORM_METHOD),
                         self.precursor_method)
        self.assertEqual(db_utils.get_meta_value(self.conn, db_utils.PROTEIN_NORM_METHOD),
                         self.protein_method)


    def test_precursor_medians_equal(self):
        self.assertIsNotNone(self.conn)
        self.assertTrue(precursor_medians_equal(self.conn))


    def test_command_log_updated(self):
        self.assertIsNotNone(self.conn)

        last_command = db_utils.get_last_command(self.conn)
        last_command[0] = os.path.basename(last_command[0])

        self.assertEqual(last_command, self.command)


class TestSingleProject(CommonTests):
    TEST_PROJECT = 'Sp3'

    @classmethod
    def run_commands(cls):
        cls.parse_result = setup_functions.setup_single_db(cls.data_dir,
                                                           cls.work_dir,
                                                           cls.TEST_PROJECT,
                                                           clear_dir=True)

        if cls.parse_result.returncode != 0:
            raise RuntimeError('Setup of test db failed!')

        cls.normalize_result = setup_functions.run_command(cls.command,
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
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_normalize_db_median_single_project'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        cls.precursor_method = 'median'
        cls.protein_method = 'median'

        cls.command = ['dia_qc', 'normalize', '-m=median', cls.db_path]
        cls.run_commands()


    def test_protein_medians_equal(self):
        self.assertIsNotNone(self.conn)
        self.assertTrue(protein_medians_eqeual(self.conn))


class TestDirectLFQSingle(unittest.TestCase, TestSingleProject):
    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_normalize_db_DirectLFQ_single_project'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        cls.precursor_method = 'median'
        cls.protein_method = 'DirectLFQ'

        cls.command = ['dia_qc', 'normalize', '-m=DirectLFQ', cls.db_path]
        cls.run_commands()


    def test_keep_missing_fails(self):
        command = ['dia_qc', 'normalize', '-m=DirectLFQ', '--keepMissing', self.db_path]

        result = setup_functions.run_command(command, self.work_dir,
                                             prefix='test_keepMissing')
        self.assertEqual(result.returncode, 1)
        self.assertTrue('--keepMissing option not compatible with DirectLFQ' in result.stderr)


class TestMedianDiaNN(unittest.TestCase, CommonTests):
    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_normalize_db_median_DiaNN_keep_missing'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        cls.precursor_method = 'median'
        cls.protein_method = 'median'

        # set up test db
        setup_functions.make_work_dir(cls.work_dir, True)
        parse_command = ['dia_qc', 'parse', '--projectName=Sp3',
                         '-m', f'{cls.data_dir}/metadata/Sp3_metadata.json',
                         f'{cls.data_dir}/skyline_reports/Sp3_replicate_quality.tsv',
                         f'{cls.data_dir}/skyline_reports/Sp3_DiaNN_precursor_quality.tsv']
        cls.parse_result = setup_functions.run_command(parse_command, cls.work_dir, prefix='parse')

        if cls.parse_result.returncode != 0:
            raise RuntimeError('Setup of test db failed!')

        cls.command = ['dia_qc', 'normalize', '-m=median', '--keepMissing', cls.db_path]
        cls.normalize_result = setup_functions.run_command(cls.command, cls.work_dir,
                                                           prefix='normalize')

        cls.conn = None
        if cls.normalize_result.returncode == 0:
            if os.path.isfile(cls.db_path):
                cls.conn = sqlite3.connect(cls.db_path)


    def test_NAs_empty(self):
        '''
        Test that a protein's quantity is NULL when all its precursor quantities are NULL.
        '''
        self.assertIsNotNone(self.conn)
        self.assertTrue(db_utils.is_normalized(self.conn))

        df_pre = pd.read_sql('''
            SELECT
                p.replicateId, ptp.proteinId, p.peptideId, p.precursorCharge,
                p.totalAreaFragment as area, p.normalizedArea
            FROM precursors p
            LEFT JOIN peptideToProtein ptp ON ptp.peptideId = p.peptideId;''', self.conn)

        self.assertSeriesEqual(df_pre['area'].isna(), df_pre['normalizedArea'].isna(),
                               check_names=False)

        protein_keys = ['replicateId', 'proteinId']
        na_precursors = df_pre.groupby(protein_keys)['area'].apply(lambda x: all(pd.isna(x)))
        na_precursors.name = 'all_na'

        df_prot = pd.read_sql('''
            SELECT
                q.replicateId, q.proteinId, q.abundance, q.normalizedAbundance
            FROM proteinQuants q
            LEFT JOIN proteins prot ON prot.proteinId == q.proteinId;''', self.conn)

        self.assertSeriesEqual(df_prot['abundance'].isna(),
                               df_prot['normalizedAbundance'].isna(),
                               check_names=False)

        df_prot = df_prot.set_index(protein_keys).join(na_precursors)
        self.assertFalse(any(df_prot['all_na'].isna()))
        self.assertSeriesEqual(df_prot['all_na'], df_prot['abundance'].isna(),
                               check_names=False)


class TestMultiProject(CommonTests):
    @classmethod
    def run_commands(cls):
        if any(result.returncode != 0 for result in cls.parse_results):
            raise RuntimeError('Setup of test db failed!')

        cls.normalize_result = setup_functions.run_command(cls.command,
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
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_normalize_db_median_multi_project'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'
        cls.parse_results = setup_functions.setup_multi_db(cls.data_dir, cls.work_dir)

        cls.precursor_method = 'median'
        cls.protein_method = 'median'

        cls.command = ['dia_qc', 'normalize', '-m=median', cls.db_path]
        cls.run_commands()


class TestDirectLFQMulti(unittest.TestCase, TestMultiProject):
    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_normalize_db_DirectLFQ_multi_project'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'
        cls.parse_results = setup_functions.setup_multi_db(cls.data_dir, cls.work_dir)

        cls.precursor_method = 'median'
        cls.protein_method = 'DirectLFQ'

        cls.command = ['dia_qc', 'normalize', '-m=DirectLFQ', cls.db_path]
        cls.run_commands()


    def test_keep_missing_fails(self):
        command = ['dia_qc', 'normalize', '-m=DirectLFQ', '--keepMissing', self.db_path]

        result = setup_functions.run_command(command, self.work_dir,
                                             prefix='test_keepMissing')
        self.assertEqual(result.returncode, 1)
        self.assertTrue('--keepMissing option not compatible with DirectLFQ' in result.stderr)


class TestAllPrecursorsMissing(unittest.TestCase):
    TEST_PROJECT = 'GPF'

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_normalize_db_missing_precursors'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        command = ['dia_qc', 'parse', f'--projectName={cls.TEST_PROJECT}',
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


    def test_median_normalization_na_rm(self):
        self.assertIsNotNone(self.conn)

        norm_levels = [db_utils.PRECURSOR_NORM_METHOD, db_utils.PROTEIN_NORM_METHOD]
        cur = self.conn.cursor()
        for level in norm_levels:
            cur.execute('DELETE FROM metadata WHERE key == ?', (level,))
        cur.execute("UPDATE metadata SET value == FALSE WHERE key == ?", (db_utils.IS_NORMALIZED, ))
        self.conn.commit()

        self.assertFalse(db_utils.is_normalized(self.conn))
        with self.assertLogs(db_utils.LOGGER, level='ERROR') as cm:
            for level in norm_levels:
                self.assertIsNone(db_utils.get_meta_value(self.conn, level))

        normalize_command = ['dia_qc', 'normalize', '-m=median', self.db_path]
        normalize_result = setup_functions.run_command(normalize_command,
                                                       self.work_dir,
                                                       prefix='normalize_db_median_na_rm')

        self.assertEqual(normalize_result.returncode, 0)
        self.assertTrue(db_utils.is_normalized(self.conn))
        for level in norm_levels:
            self.assertEqual(db_utils.get_meta_value(self.conn, level), 'median')

        # make sure all normalized precursor areas are NULL
        cur = self.conn.cursor()
        cur.execute('SELECT * FROM precursors WHERE normalizedArea IS NOT NULL')
        self.assertEqual(len(cur.fetchall()), 0)

        # make sure all normalized protein quants are NULL
        cur = self.conn.cursor()
        cur.execute('SELECT * FROM proteinQuants WHERE normalizedAbundance IS NOT NULL')
        self.assertEqual(len(cur.fetchall()), 0)


    def test_median_normalization_keep_na(self):
        self.assertIsNotNone(self.conn)

        norm_levels = [db_utils.PRECURSOR_NORM_METHOD, db_utils.PROTEIN_NORM_METHOD]
        cur = self.conn.cursor()
        for level in norm_levels:
            cur.execute('DELETE FROM metadata WHERE key == ?', (level,))
        cur.execute("UPDATE metadata SET value == FALSE WHERE key == ?", (db_utils.IS_NORMALIZED,))
        self.conn.commit()

        self.assertFalse(db_utils.is_normalized(self.conn))
        with self.assertLogs(db_utils.LOGGER, level='ERROR') as cm:
            for level in norm_levels:
                self.assertIsNone(db_utils.get_meta_value(self.conn, level))

        normalize_command = ['dia_qc', 'normalize', '--keepMissing', '-m=median', self.db_path]
        normalize_result = setup_functions.run_command(normalize_command,
                                                       self.work_dir,
                                                       prefix='normalize_db_median_keep_na')

        self.assertEqual(normalize_result.returncode, 0)
        self.assertTrue(db_utils.is_normalized(self.conn))
        for level in norm_levels:
            self.assertEqual(db_utils.get_meta_value(self.conn, level), 'median')

        # make sure all normalized precursor areas are not NULL
        cur = self.conn.cursor()
        cur.execute('SELECT totalAreaFragment, normalizedArea FROM precursors WHERE normalizedArea IS NULL')
        na_precursors = [norm_area for area, norm_area in cur.fetchall() if not (area is None and norm_area is None)]
        self.assertEqual(len(na_precursors), 0)

        # make sure all normalized protein quants are not NULL
        cur = self.conn.cursor()
        cur.execute('SELECT abundance, normalizedAbundance FROM proteinQuants WHERE normalizedAbundance IS NULL')
        na_proteins = cur.fetchall()
        na_proteins = [norm_abundance for abundance, norm_abundance in cur.fetchall()
                       if not (abundance is None and norm_abundance is None)]
        self.assertEqual(len(na_proteins), 0)

        # make sure protein and precursor medians are equal
        self.assertTrue(precursor_medians_equal(self.conn))
        self.assertTrue(protein_medians_eqeual(self.conn))


    def test_DirectLFQ_normalization(self):
        normalize_command = ['dia_qc', 'normalize', '-m=DirectLFQ', self.db_path]

        normalize_result = setup_functions.run_command(normalize_command,
                                                       self.work_dir,
                                                       prefix='normalize_db_DirectLFQ')

        self.assertEqual(normalize_result.returncode, 1)
        self.assertTrue('ERROR: Can not perform DirectLFQ normalization with 0 precursors!' in normalize_result.stderr)
