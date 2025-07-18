
import unittest
import os
import sqlite3
from collections import Counter
import pandas as pd

import setup_functions

import DIA_QC_report.submodules.normalization as normalization
import DIA_QC_report.submodules.dia_db_utils as db_utils
from DIA_QC_report.submodules.skyline_reports import PrecursorReport
from DIA_QC_report import parse_data


class TestNormalizationBase(setup_functions.AbstractTestsBase):
    def __init__(self):
        self.work_dir = None
        self.db_path = None
        self.data_dir = None
        self.parse_result = None
        self.conn = None


    @classmethod
    def tearDownClass(cls):
        if cls.conn is not None:
            cls.conn.close()


    def check_medians_equal(self, df, quant_col):
        norm_col = normalization.cammel_case('normalized', quant_col)
        log2_norm_col = normalization.cammel_case('log2', norm_col)

        medians = df.groupby('replicateId')[log2_norm_col].median()
        for m in medians:
            self.assertAlmostEqual(medians.iloc[0], m)


    def test_abc_instantiation_fails(self):
        self.assertIsNotNone(self.conn)

        with self.assertRaises(TypeError):
            normalization.NormalizationManagerBase(self.conn)


    def test_median_normalization_na_rm(self):
        self.assertIsNotNone(self.conn)

        manager = normalization.MedianNormalizer(self.conn)
        with self.assertLogs(normalization.LOGGER) as cm:
            self.assertTrue(manager.normalize())

        precursors, proteins = manager.get_long_tables()

        self.check_medians_equal(precursors, 'area')
        self.check_medians_equal(proteins, 'abundance')


    def test_median_normalization_keep_missing(self):
        self.assertIsNotNone(self.conn)

        manager = normalization.MedianNormalizer(self.conn, keep_na=True)
        with self.assertLogs(normalization.LOGGER) as cm:
            self.assertTrue(manager.normalize())

        precursors, proteins = manager.get_long_tables()

        self.check_medians_equal(precursors, 'area')
        self.check_medians_equal(proteins, 'abundance')


    def test_directlfq_normalization(self):
        self.assertIsNotNone(self.conn)

        manager = normalization.DirectlfqNormalizer(self.conn)
        with self.assertLogs(normalization.LOGGER) as cm:
            self.assertTrue(manager.normalize())

        precursors, proteins = manager.get_long_tables()

        self.check_medians_equal(precursors, 'area')


    def test_get_long_tables(self):
        self.assertIsNotNone(self.conn)

        manager = normalization.MedianNormalizer(self.conn)
        with self.assertLogs(normalization.LOGGER) as cm:
            self.assertTrue(manager.normalize())

        precursors, proteins = manager.get_long_tables(use_db_ids=False)

        self.assertTrue('replicate' in precursors.columns)
        self.assertTrue('replicate' in proteins.columns)
        self.assertTrue('protein' in proteins.columns)
        self.assertEqual(proteins.shape[1], 8)
        self.assertEqual(precursors.shape[1], 9)

        precursors, proteins = manager.get_long_tables(use_db_ids=True)

        self.assertFalse('replicate' in precursors.columns)
        self.assertFalse('replicate' in proteins.columns)
        self.assertFalse('protein' in proteins.columns)
        self.assertEqual(proteins.shape[1], 6)
        self.assertEqual(precursors.shape[1], 8)


    def test_all_reps_skipped(self):
        self.assertIsNotNone(self.conn)

        test_conn = sqlite3.connect(':memory:')
        self.conn.backup(test_conn)
        manager = normalization.MedianNormalizer(test_conn)

        try:
            cur = test_conn.cursor()
            cur.execute('UPDATE replicates SET includeRep = FALSE;')
            test_conn.commit()

            with self.assertLogs(normalization.LOGGER, level='ERROR') as cm:
                self.assertFalse(manager._read_precursors())
            self.assertTrue(any('All replicates in database have been excluded!' in entry
                                for entry in cm.output))

            with self.assertLogs(normalization.LOGGER, level='ERROR') as cm:
                self.assertFalse(manager.normalize())
            self.assertTrue(any('All replicates in database have been excluded!' in entry
                                for entry in cm.output))

        finally:
            test_conn.close()


class TestSingleNormalization(unittest.TestCase, TestNormalizationBase):
    TEST_PROJECT = 'Sp3'

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_normalization_functions/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        cls.parse_result = setup_functions.setup_single_db(cls.data_dir,
                                                           cls.work_dir,
                                                           cls.TEST_PROJECT,
                                                           clear_dir=True)

        cls.conn = None
        if cls.parse_result.returncode == 0:
            if os.path.isfile(cls.db_path):
                cls.conn = sqlite3.connect(cls.db_path)


    def test_read_precursors(self):
        self.assertIsNotNone(self.conn)

        manager = normalization.MedianNormalizer(self.conn)

        with self.assertNoLogs(normalization.LOGGER, level='WARNING') as cm:
            self.assertTrue(manager._read_precursors())


class TestMedianDiaNN(unittest.TestCase, TestNormalizationBase):
    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_normalization_functions_DiaNN'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        # set up test db
        setup_functions.make_work_dir(cls.work_dir, True)
        parse_command = ['--projectName=Sp3', '-m', f'{cls.data_dir}/metadata/Sp3_metadata.json',
                         f'{cls.data_dir}/skyline_reports/Sp3_replicate_quality.tsv',
                         f'{cls.data_dir}/skyline_reports/Sp3_DiaNN_precursor_quality.tsv']
        cls.parse_result = setup_functions.run_main(
            parse_data._main, parse_command, cls.work_dir, prefix='parse', prog='dia_qc parse'
        )

        cls.conn = None
        if cls.parse_result.returncode == 0:
            if os.path.isfile(cls.db_path):
                cls.conn = sqlite3.connect(cls.db_path)


    def test_na_precursors(self):
        self.assertIsNotNone(self.conn)

        manager = normalization.MedianNormalizer(self.conn, keep_na=True)

        with self.assertNoLogs(normalization.LOGGER, level='WARNING') as cm:
            manager.normalize()

        df_pre = manager.precursors.copy()
        self.assertSeriesEqual(df_pre['area'].isna(), df_pre['normalizedArea'].isna(),
                               check_names=False)

        # add proteinID column back to precursor df
        df_ptp = pd.read_sql('SELECT peptideId, proteinId FROM peptideToProtein;', self.conn)
        df_pre = df_pre.set_index('peptideId')
        df_ptp = df_ptp.set_index('peptideId')
        df_pre = df_pre.join(df_ptp).reset_index()

        protein_keys = ['replicateId', 'proteinId']
        df_pre = df_pre.set_index(protein_keys)
        na_precursors = df_pre.groupby(protein_keys)['area'].apply(lambda x: all(pd.isna(x)))
        na_precursors.name = 'all_na'

        df_prot = manager.proteins.copy()
        self.assertSeriesEqual(df_prot['abundance'].isna(),
                               df_prot['normalizedAbundance'].isna(),
                               check_names=False)

        df_prot = df_prot.set_index(protein_keys).join(na_precursors)
        self.assertFalse(any(df_prot['all_na'].isna()))
        self.assertSeriesEqual(df_prot['all_na'], df_prot['abundance'].isna(),
                               check_names=False)


class TestAllPrecursorsMissing(unittest.TestCase, TestNormalizationBase):
    TEST_PROJECT = 'GPF'

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_normalization_functions_missing/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        command = [f'--projectName={cls.TEST_PROJECT}',
                   f'{cls.data_dir}/skyline_reports/{cls.TEST_PROJECT}_replicate_quality.tsv',
                   f'{cls.data_dir}/skyline_reports/{cls.TEST_PROJECT}_precursor_quality.tsv']

        setup_functions.make_work_dir(cls.work_dir, True)
        cls.parse_result = setup_functions.run_main(
            parse_data._main, command, cls.work_dir, prefix='parse_missing_precursors', prog='dia_qc parse'
        )

        cls.conn = None
        if cls.parse_result.returncode == 0:
            if os.path.isfile(cls.db_path):
                cls.conn = sqlite3.connect(cls.db_path)


    def test_directlfq_normalization(self):
        self.assertIsNotNone(self.conn)

        manager = normalization.DirectlfqNormalizer(self.conn)
        with self.assertLogs(normalization.LOGGER) as cm:
            self.assertFalse(manager.normalize())

        self.assertTrue('Can not perform DirectLFQ normalization with 0 precursors!' in cm.output[-1])


    def test_read_precursors_na_rm(self):
        self.assertIsNotNone(self.conn)

        # read ground truth data
        df = pd.read_csv(f'{self.data_dir}/skyline_reports/{self.TEST_PROJECT}_precursor_quality.tsv', sep='\t')
        cols = {col.skyline_name: col.name for col in PrecursorReport().required_columns()}
        df = df.rename(columns=cols)
        df = df[cols.values()]
        df = df[['replicateName', 'modifiedSequence', 'precursorCharge',
                 'totalAreaFragment']].drop_duplicates()
        df_w = df.pivot(index=['modifiedSequence', 'precursorCharge'],
                        columns='replicateName', values='totalAreaFragment')

        prec_counts = df_w.apply(lambda x: any(x.isnull().values), axis=1)
        na_counts = Counter(prec_counts)
        prec_counts = prec_counts[prec_counts == False]

        n_missing = na_counts[True]
        n_not_missing = na_counts[False]
        n_precursors = sum(na_counts.values())

        manager = normalization.MedianNormalizer(self.conn)
        with self.assertLogs(normalization.LOGGER, level='WARNING') as cm:
            self.assertTrue(manager._read_precursors())

        # check that log shows correct number of missing values
        self.assertTrue(f'Removed {n_missing} of {n_precursors} precursors with missing values.' in l for l in cm.output)
        self.assertTrue(f'{n_not_missing} precursors without missing values remain.' in l for l in cm.output)

        # check that db and ground truth precursor dfs match
        # make gt dict
        df_gt = df.set_index(['modifiedSequence', 'precursorCharge'])
        df_gt = df_gt[df_gt.index.isin(prec_counts.index)]
        df_gt = df_gt.reset_index()
        df_gt = df_gt.set_index(['replicateName', 'modifiedSequence', 'precursorCharge'])
        df_gt = {row.Index: row.totalAreaFragment for row in df_gt.itertuples()}

        # make db dict
        cur = self.conn.cursor()
        cur.execute('SELECT id, replicate FROM replicates;')
        rep_ids = {x[0]: x[1] for x in cur.fetchall()}

        df_db = manager.precursors
        df_db['replicateName'] = df_db['replicateId'].apply(lambda x: rep_ids[x])
        df_db = df_db.set_index(['replicateName', 'modifiedSequence', 'precursorCharge'])
        df_db = {row.Index: row.area for row in df_db.itertuples()}

        self.assertEqual(len(df_db), len(df_gt))
        self.assertDictEqual(df_db, df_gt)


class TestMultiNormalization(unittest.TestCase, TestNormalizationBase):
    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_normalization_functions_multi/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        cls.parse_results = setup_functions.setup_multi_db(cls.data_dir, cls.work_dir)

        cls.conn = None
        if all(result.returncode == 0 for result in cls.parse_results):
            if os.path.isfile(cls.db_path):
                cls.conn = sqlite3.connect(cls.db_path)


    def test_read_precursors_na_rm(self):
        self.assertIsNotNone(self.conn)

        # read ground truth data
        dfs = [f'{self.data_dir}/skyline_reports/Sp3_by_protein_precursor_quality.tsv',
               f'{self.data_dir}/skyline_reports/Strap_by_protein_precursor_quality.tsv']
        cols = {col.skyline_name: col.name for col in PrecursorReport().required_columns()}
        for i in range(len(dfs)):
            dfs[i] = pd.read_csv(dfs[i], sep='\t')
            dfs[i] = dfs[i].rename(columns=cols)
            dfs[i] = dfs[i][cols.values()]
            dfs[i] = dfs[i][['replicateName', 'modifiedSequence', 'precursorCharge',
                             'totalAreaFragment']].drop_duplicates()
        df = pd.concat(dfs)
        df_w = df.pivot(index=['modifiedSequence', 'precursorCharge'],
                        columns='replicateName', values='totalAreaFragment')

        prec_counts = df_w.apply(lambda x: any(x.isnull().values), axis=1)
        na_counts = Counter(prec_counts)
        prec_counts = prec_counts[prec_counts == False]

        n_missing = na_counts[True]
        n_not_missing = na_counts[False]
        n_precursors = sum(na_counts.values())

        manager = normalization.MedianNormalizer(self.conn)
        with self.assertLogs(normalization.LOGGER, level='WARNING') as cm:
            self.assertTrue(manager._read_precursors())

        # check that log shows correct number of missing values
        self.assertTrue(f'Removed {n_missing} of {n_precursors} precursors with missing values.' in l for l in cm.output)
        self.assertTrue(f'{n_not_missing} precursors without missing values remain.' in l for l in cm.output)

        # check that db and ground truth precursor dfs match
        # make gt dict
        df_gt = df.set_index(['modifiedSequence', 'precursorCharge'])
        df_gt = df_gt[df_gt.index.isin(prec_counts.index)]
        df_gt = df_gt.reset_index()
        df_gt = df_gt.set_index(['replicateName', 'modifiedSequence', 'precursorCharge'])
        df_gt = {row.Index: row.totalAreaFragment for row in df_gt.itertuples()}

        # make db dict
        cur = self.conn.cursor()
        cur.execute('SELECT id, replicate FROM replicates;')
        rep_ids = {x[0]: x[1] for x in cur.fetchall()}

        df_db = manager.precursors
        df_db['replicateName'] = df_db['replicateId'].apply(lambda x: rep_ids[x])
        df_db = df_db.set_index(['replicateName', 'modifiedSequence', 'precursorCharge'])
        df_db = {row.Index: row.area for row in df_db.itertuples()}

        self.assertEqual(len(df_db), len(df_gt))
        self.assertDictEqual(df_db, df_gt)
