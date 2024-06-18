
import unittest
import os
import sqlite3
import pandas as pd
from numpy import log2

import setup_functions

import DIA_QC_report.submodules.dia_db_utils as db_utils

class TestSingleProject(unittest.TestCase):
    TEST_PROJECT = 'Sp3'

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_normalize_single_project/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        cls.parse_result = setup_functions.setup_single_db(cls.data_dir,
                                                           cls.work_dir,
                                                           cls.TEST_PROJECT,
                                                           clear_dir=True)

        if cls.parse_result.returncode != 0:
            raise RuntimeError('Setup of test db failed!')

        normalize_command = ['normalize_db', cls.db_path]
        cls.normalize_result = setup_functions.run_command(normalize_command,
                                                           cls.work_dir,
                                                           prefix='normalize_single_proj')

        cls.conn = None
        if cls.normalize_result.returncode == 0:
            if os.path.isfile(cls.db_path):
                cls.conn = sqlite3.connect(cls.db_path)


    @classmethod
    def tearDownClass(cls):
        if cls.conn is not None:
            cls.conn.close()


    def test_is_successful(self):
        self.assertEqual(self.parse_result.returncode, 0)
        self.assertEqual(self.normalize_result.returncode, 0)


    def test_is_normalized(self):
        self.assertTrue(self.conn is not None)
        self.assertTrue(db_utils.is_normalized(self.conn))
        self.assertEqual(db_utils.get_meta_value(self.conn, 'precursor_normalization_method'),
                         'median')
        self.assertEqual(db_utils.get_meta_value(self.conn, 'protein_normalization_method'),
                        'DirectLFQ')


    def test_precursor_medians_equal(self):
        self.assertTrue(self.conn is not None)

        query = ''' SELECT
                r.replicate,
                p.modifiedSequence,
                p.precursorCharge,
                p.normalizedArea
            FROM precursors p
            LEFT JOIN replicates r ON r.id == p.replicateId; '''
        df = pd.read_sql(query, self.conn)

        self.assertTrue(not any(df['normalizedArea'] == 0))

        df['log2NormArea'] = log2(df['normalizedArea'])
        medians = df.groupby('replicate')['log2NormArea'].median()

        self.assertTrue(all(medians.iloc[0] == medians))


class TestMultiProject(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_normalize_multi_project/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'
        cls.parse_results = setup_functions.setup_multi_db(cls.data_dir, cls.work_dir)
        
        if any(result.returncode != 0 for result in cls.parse_results):
            raise RuntimeError('Setup of test db failed!')

        normalize_command = ['normalize_db', cls.db_path]
        cls.normalize_result = setup_functions.run_command(normalize_command,
                                                           cls.work_dir,
                                                           prefix='normalize_multi_proj')

        cls.conn = None
        if cls.normalize_result.returncode == 0:
            if os.path.isfile(cls.db_path):
                cls.conn = sqlite3.connect(cls.db_path)


    @classmethod
    def tearDownClass(cls):
        if cls.conn is not None:
            cls.conn.close()


    def test_is_successful(self):
        self.assertTrue(all(result.returncode == 0 for result in self.parse_results))
        self.assertEqual(self.normalize_result.returncode, 0)


    def test_is_normalized(self):
        self.assertTrue(self.conn is not None)
        self.assertTrue(db_utils.is_normalized(self.conn))
        self.assertEqual(db_utils.get_meta_value(self.conn, 'precursor_normalization_method'),
                         'median')
        self.assertEqual(db_utils.get_meta_value(self.conn, 'protein_normalization_method'),
                        'DirectLFQ')


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

        df['log2NormArea'] = log2(df['normalizedArea'])
        medians = df.groupby('replicate')['log2NormArea'].median()

        self.assertTrue(all(medians.iloc[0] == medians))


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


if __name__ == '__main__':
    unittest.main()
