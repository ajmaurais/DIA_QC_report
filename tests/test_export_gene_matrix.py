
import unittest
import os
import re

import pandas as pd

import setup_functions

from DIA_QC_report import parse_data
from DIA_QC_report import normalize_db
from DIA_QC_report import export_gene_matrix

class TestExportGeneMatrix(unittest.TestCase):
    TEST_PROJECT = 'Strap'

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_export_gene_matrix/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'
        cls.gene_id_path = f'{cls.data_dir}/metadata/prhuman2gene_2023_05_24_subset.csv'
        cls.prog = 'dia_qc export_gene_matrix'

        parse_result = setup_functions.setup_single_db(
            cls.data_dir, cls.work_dir, cls.TEST_PROJECT, clear_dir=True, group_by_gene=True
        )
        if parse_result.returncode != 0:
            raise RuntimeError('Setup of test db failed!')


    def assertIsFile(self, path):
        if not os.path.isfile(path):
            raise AssertionError(f'File does not exist: {path}')


    def test_is_successful(self):
        prefix='test_is_sucessful'
        command = [f'--prefix={prefix}', self.gene_id_path, self.db_path]
        result = setup_functions.run_main(
            export_gene_matrix._main, command, self.work_dir, prog=self.prog
        )
        self.assertEqual(result.returncode, 0, result.stderr)

        # make sure expected reports are generated
        for file in ('precursors', 'proteins'):
            self.assertIsFile(f'{self.work_dir}/{prefix}_{file}_unnormalized.tsv')


    def test_gene_table_tsv(self):
        prefix='test_tsv'
        command = [f'--prefix={prefix}',
                   f'{self.data_dir}/metadata/prhuman2gene_2023_05_24_subset.tsv',
                   self.db_path]
        result = setup_functions.run_main(
            export_gene_matrix._main, command, self.work_dir, prog=self.prog
        )
        self.assertEqual(result.returncode, 0, result.stderr)

        # make sure expected reports are generated
        for file in ('precursors', 'proteins'):
            self.assertIsFile(f'{self.work_dir}/{prefix}_{file}_unnormalized.tsv')


    def test_use_gene_hash_option(self):
        # make sure id table without hash fails
        command = ['--addGeneUuid', self.gene_id_path, self.db_path]
        result = setup_functions.run_main(
            export_gene_matrix._main, command, self.work_dir,
            prefix='failed_gene_id_hash_table', prog=self.prog
        )
        self.assertEqual(result.returncode, 1)

        # make sure correct table succeeds
        prefix='test_uuid'
        gene_id_hash_path = f'{self.data_dir}/metadata/prhuman2gene_gene_uuid_2023_05_24_subset.csv'
        command = [f'--prefix={prefix}', '--addGeneUuid', gene_id_hash_path, self.db_path]
        result = setup_functions.run_main(
            export_gene_matrix._main, command, self.work_dir, prog=self.prog
        )
        self.assertEqual(result.returncode, 0, result.stderr)

        # make sure expected reports are generated
        for file in ('precursors', 'proteins'):
            report_path = f'{self.work_dir}/{prefix}_{file}_unnormalized.tsv'
            self.assertIsFile(report_path)

            # make sure gene_uuid column is in reports
            df = pd.read_csv(report_path, sep='\t')
            self.assertTrue('gene_uuid' in df.columns)


class TestGroupMethod(unittest.TestCase):
    TEST_PROJECT = 'Strap'

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_export_gene_matrix_group_by'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data'
        cls.gene_id_path = f'{cls.data_dir}/metadata/prhuman2gene_2023_05_24_subset.csv'
        cls.prog = 'dia_qc export_gene_matrix'

        setup_functions.make_work_dir(cls.work_dir, clear_dir=True)
        argv = [
            '--groupBy=gene',
            f'{cls.data_dir}/skyline_reports/PDC_NF_PNET_Validation_Study_replicate_quality.tsv',
            f'{cls.data_dir}/skyline_reports/PDC_NF_PNET_Validation_Study_precursor_quality.tsv'
        ]
        parse_result = setup_functions.run_main(
            parse_data._main, argv, cls.work_dir, prog='dia_qc parse', prefix='parse'
        )
        if parse_result.returncode != 0:
            raise RuntimeError('Setup of test db failed!')

        normalize_result = setup_functions.run_main(
            normalize_db._main, [cls.db_path], cls.work_dir,
            prog='dia_qc normalize', prefix='normalize'
        )
        if normalize_result.returncode != 0:
            raise RuntimeError('Normalization of test db failed!')


    def test_split(self):
        prefix = 'test_split'
        command = [f'--prefix={prefix}', '--groupMethod=split', self.gene_id_path, self.db_path]
        result = setup_functions.run_main(
            export_gene_matrix._main, command, self.work_dir, prog=self.prog
        )
        self.assertEqual(result.returncode, 0, result.stderr)

        protein_id_cols = ['Gene', 'GeneGroup']
        precursor_id_cols = protein_id_cols + ['modifiedSequence', 'precursorCharge']

        for file in ('unnormalized', 'normalized'):
            df_path = f'{self.work_dir}/{prefix}_precursors_{file}.tsv'
            with self.subTest(f'Check {prefix}_precursors_{file}'):
                self.assertTrue(os.path.isfile(df_path), f'File does not exist: {df_path}')

                df = pd.read_csv(df_path, sep='\t')
                self.assertEqual(len(df.index), len(df[precursor_id_cols].drop_duplicates().index))
                df = df[precursor_id_cols].drop_duplicates()
                df.loc[:, 'precursor'] = df.modifiedSequence.astype(str).str.cat(df.precursorCharge.astype(str), sep='_')
                df = df.drop(columns=['modifiedSequence', 'precursorCharge'])

                for group, group_df in df.groupby(['GeneGroup', 'precursor']):
                    genes = re.split(r';\s*', group[0])
                    self.assertEqual(len(genes), len(group_df.index))

        for file in ('unnormalized', 'normalized'):
            df_path = f'{self.work_dir}/{prefix}_proteins_{file}.tsv'
            with self.subTest(f'Check {prefix}_proteins_{file}'):
                self.assertTrue(os.path.isfile(df_path), f'File does not exist: {df_path}')

                df = pd.read_csv(df_path, sep='\t')
                self.assertEqual(len(df.index), len(df[protein_id_cols].drop_duplicates()))
                df = df[protein_id_cols].drop_duplicates()

                for group, group_df in df.groupby('GeneGroup'):
                    genes = re.split(r';\s*', group)
                    self.assertEqual(len(genes), len(group_df.index))


    def test_drop(self):
        prefix = 'test_drop'
        command = [f'--prefix={prefix}', '--groupMethod=drop', self.gene_id_path, self.db_path]
        result = setup_functions.run_main(
            export_gene_matrix._main, command, self.work_dir, prog=self.prog
        )
        self.assertEqual(result.returncode, 0, result.stderr)

        protein_id_cols = ['Gene', 'GeneGroup']
        precursor_id_cols = protein_id_cols + ['modifiedSequence', 'precursorCharge']

        for level in ('precursor', 'protein'):
            for file in ('unnormalized', 'normalized'):
                df_path = f'{self.work_dir}/{prefix}_{level}s_{file}.tsv'
                with self.subTest(f'Check {prefix}_{level}s_{file}'):
                    self.assertTrue(os.path.isfile(df_path), f'File does not exist: {df_path}')

                    df = pd.read_csv(df_path, sep='\t')
                    if level == 'precursor':
                        self.assertEqual(len(df.index), len(df[precursor_id_cols].drop_duplicates().index))
                    else:
                        self.assertEqual(len(df.index), len(df[protein_id_cols].drop_duplicates().index))

                    gene_groups = df.GeneGroup.drop_duplicates().to_list()
                    for group in gene_groups:
                        self.assertEqual(1, len(re.split(r';\s*', group)))


if __name__ == '__main__':
    unittest.main()