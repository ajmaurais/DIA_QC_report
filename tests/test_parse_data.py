
import unittest
from unittest import mock
import os
import sqlite3
import json
import re
from collections import Counter
import pandas as pd

import setup_functions

from DIA_QC_report.submodules.metadata import Dtype
from DIA_QC_report.submodules.dia_db_utils import is_normalized
from DIA_QC_report.submodules.dia_db_utils import update_meta_value, get_meta_value

TEST_DIR = os.path.dirname(os.path.abspath(__file__))


def get_db_metadata(db_path):
    if not os.path.isfile(db_path):
        raise RuntimeError(f"'{db_path}' does not exist!")

    conn = sqlite3.connect(db_path)

    # get replicates
    query = 'SELECT replicateId, replicate, project FROM replicates;'
    cur = conn.cursor()
    cur.execute(query)
    replicates = {row[0]: {'replicate': row[1], 'project': row[2]} for row in cur.fetchall()}

    # get metadata
    query = 'SELECT * FROM sampleMetadata;'
    cur = conn.cursor()
    cur.execute(query)

    # pivot metadata wider
    metadata = dict()
    for row in cur.fetchall():
        if row[0] not in metadata:
            metadata[row[0]] = list()
        metadata[row[0]].append((row[1], row[2]))

    # get metadata types
    query = 'SELECT annotationKey, annotationType FROM sampleMetadataTypes;'
    cur = conn.cursor()
    cur.execute(query)
    types = {row[0]: row[1] for row in cur.fetchall()}

    conn.close()

    ret = {}
    for replicateId, data in replicates.items():
        ret[data['replicate']] = {'project': data['project']} | dict(metadata[replicateId])

    return ret, types


def load_metadata():
    metadata = dict()
    for project in ('Sp3', 'Strap'):
        with open(f'{TEST_DIR}/data/metadata/{project}_metadata.json', 'r') as inF:
            metadata[project] = json.load(inF)

    metadata_table = {}
    for project, replicates in metadata.items():
        for replicate, annotations in replicates.items():
            metadata_table[replicate] = {'project': project}
            for variable, value in annotations.items():
                metadata_table[replicate][variable] = value

    return metadata_table


METADATA = load_metadata()
METADATA_TYPES = {'cellLine': 'STRING',
                  'sample.name': 'STRING',
                  'experiment': 'STRING',
                  'NCI7.std': 'BOOL'}


class TestMetadata(unittest.TestCase):
    TEST_PROJECT = 'Sp3'

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{TEST_DIR}/work/test_metadata'
        setup_functions.make_work_dir(cls.work_dir, clear_dir=True)


    @staticmethod
    def setup_command(project, ext, meta_suffix='_metadata'):
        db_path = f'test_{ext}{meta_suffix}.db3'
        command = ['parse_data',
                   f'--projectName={project}', '--overwriteMode=overwrite',
                   '-m', f'{TEST_DIR}/data/metadata/{project}{meta_suffix}.{ext}',
                   '-o', db_path,
                   f'{TEST_DIR}/data/skyline_reports/{project}_replicate_quality.tsv',
                   f'{TEST_DIR}/data/skyline_reports/{project}_by_protein_precursor_quality.tsv']
        return command, db_path


    def test_tsv(self):
        command, db_path = self.setup_command(self.TEST_PROJECT, 'tsv')
        result = setup_functions.run_command(command, self.work_dir)

        # test command was successful
        self.assertEqual(0, result.returncode)

        # get metadata from test db
        db_metadata, db_metadata_types = get_db_metadata(f'{self.work_dir}/{db_path}')

        # test all replicates are in test db
        project_reps = [rep for rep, data in METADATA.items() if data['project'] == self.TEST_PROJECT]
        self.assertEqual(len(project_reps), len(db_metadata))
        for rep in project_reps:
            self.assertTrue(rep in db_metadata)

        # test all annotations are correct
        for rep in project_reps:
            self.assertDictEqual(db_metadata[rep], METADATA[rep])

        # test metadata types are correct
        for variable in db_metadata_types:
            self.assertEqual(db_metadata_types[variable], METADATA_TYPES[variable])


    def test_json(self):
        command, db_path = self.setup_command(self.TEST_PROJECT, 'json')
        result = setup_functions.run_command(command, self.work_dir)

        # test command was successful
        self.assertEqual(0, result.returncode)

        # get metadata from test db
        db_metadata, db_metadata_types = get_db_metadata(f'{self.work_dir}/{db_path}')

        # test all replicates are in test db
        project_reps = [rep for rep, data in METADATA.items() if data['project'] == self.TEST_PROJECT]
        self.assertEqual(len(project_reps), len(db_metadata))
        for rep in project_reps:
            self.assertTrue(rep in db_metadata)

        # test all annotations are correct
        for rep in project_reps:
            self.assertDictEqual(db_metadata[rep], METADATA[rep])

        # test metadata types are correct
        for variable in db_metadata_types:
            self.assertEqual(db_metadata_types[variable], METADATA_TYPES[variable])


    def test_skyline_csv(self):
        command, db_path = self.setup_command(self.TEST_PROJECT, 'csv', meta_suffix='_annotations')
        result = setup_functions.run_command(command, self.work_dir)

        # test command was successful
        self.assertEqual(0, result.returncode)

        # Log entry specific for skyline annotations csv
        self.assertTrue('Found Skyline annotations csv' in result.stderr)

        # get metadata from test db
        db_metadata, db_metadata_types = get_db_metadata(f'{self.work_dir}/{db_path}')

        # test all replicates are in test db
        project_reps = [rep for rep, data in METADATA.items() if data['project'] == self.TEST_PROJECT]
        self.assertEqual(len(project_reps), len(db_metadata))
        for rep in project_reps:
            self.assertTrue(rep in db_metadata)

        # test all annotations are correct
        for rep in project_reps:
            self.assertDictEqual(db_metadata[rep], METADATA[rep])

        # test metadata types are correct
        for variable in db_metadata_types:
            self.assertEqual(db_metadata_types[variable], METADATA_TYPES[variable])


    def test_nomal_csv(self):
        command, db_path = self.setup_command(self.TEST_PROJECT, 'csv')
        result = setup_functions.run_command(command, self.work_dir)

        # test command was successful
        self.assertEqual(0, result.returncode)

        # get metadata from test db
        db_metadata, db_metadata_types = get_db_metadata(f'{self.work_dir}/{db_path}')

        # test all replicates are in test db
        project_reps = [rep for rep, data in METADATA.items() if data['project'] == self.TEST_PROJECT]
        self.assertEqual(len(project_reps), len(db_metadata))
        for rep in project_reps:
            self.assertTrue(rep in db_metadata)

        # test all annotations are correct
        for rep in project_reps:
            self.assertDictEqual(db_metadata[rep], METADATA[rep])

        # test metadata types are correct
        for variable in db_metadata_types:
            self.assertEqual(db_metadata_types[variable], METADATA_TYPES[variable])


class TestInferDtypes(unittest.TestCase):
    TEST_PROJECT = 'Strap'

    METADATA_TYPES = {'string_var': 'STRING',
                      'bool_var': 'BOOL',
                      'int_var': 'INT',
                      'float_var': 'FLOAT',
                      'na_var': 'BOOL'}

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{TEST_DIR}/work/test_infer_types'
        setup_functions.make_work_dir(cls.work_dir, clear_dir=True)


    @staticmethod
    def setup_command(project, ext, meta_suffix='_multi_var_metadata'):
        db_path = f'test_{ext}{meta_suffix}.db3'
        command = ['parse_data',
                   f'--projectName={project}', '--overwriteMode=overwrite',
                   '-m', f'{TEST_DIR}/data/metadata/{project}{meta_suffix}.{ext}',
                   '-o', db_path,
                   f'{TEST_DIR}/data/skyline_reports/{project}_replicate_quality.tsv',
                   f'{TEST_DIR}/data/skyline_reports/{project}_by_protein_precursor_quality.tsv']

        return command, db_path


    def test_infer_dtypes(self):
        command, db_path = self.setup_command(self.TEST_PROJECT, 'tsv')
        result = setup_functions.run_command(command, self.work_dir)

        # test command was successful
        self.assertEqual(0, result.returncode)

        # get metadata from test db
        db_metadata, db_metadata_types = get_db_metadata(f'{self.work_dir}/{db_path}')

        self.assertDictEqual(self.METADATA_TYPES, db_metadata_types)


class TestMultiProject(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{TEST_DIR}/work/test_multi_project/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{TEST_DIR}/data/'
        cls.results = setup_functions.setup_multi_db(cls.data_dir, cls.work_dir)

        cls.conn = None
        if any(result.returncode == 0 for result in cls.results):
            if os.path.isfile(cls.db_path):
                cls.conn = sqlite3.connect(cls.db_path)


    @classmethod
    def tearDownClass(cls):
        if cls.conn is not None:
            cls.conn.close()


    def test_is_successful(self):
        for result in self.results:
            self.assertEqual(0, result.returncode)


    @mock.patch('DIA_QC_report.submodules.dia_db_utils.LOGGER', mock.Mock())
    def test_is_not_normalized(self):
        self.assertTrue(self.conn is not None)
        self.assertFalse(is_normalized(self.conn))


    def test_replicates(self):
        self.assertTrue(self.conn is not None)

        cur = self.conn.cursor()
        cur.execute('SELECT replicate, project FROM replicates;')
        db_reps = {row[0]: row[1] for row in cur.fetchall()}

        ground_truth_reps = {replicate: data['project'] for replicate, data in METADATA.items()}
        self.assertDictEqual(ground_truth_reps, db_reps)


    def test_metadata(self):
        # get metadata from test db
        db_metadata, db_metadata_types = get_db_metadata(self.db_path)

        self.assertDictEqual(METADATA_TYPES, db_metadata_types)
        self.assertDictEqual(METADATA, db_metadata)


    def test_protein_names(self):
        self.assertTrue(self.conn is not None)

        # read ground truth proteins
        df_combined = None
        for project in ('Sp3', 'Strap'):
            df = pd.read_csv(f'{TEST_DIR}/data/skyline_reports/{project}_by_protein_precursor_quality.tsv',
                             sep='\t')
            df['project'] = project
            df = df[['ProteinName', 'project']].drop_duplicates()
            if df_combined is None:
                df_combined = df
            else:
                df_combined = pd.concat([df_combined, df])
        proteins = {(row.ProteinName, row.project) for row in df_combined.itertuples()}

        query = '''
            SELECT DISTINCT prot.name, r.project FROM proteinQuants q
            LEFT JOIN proteins prot ON prot.proteinID == q.proteinId
            LEFT JOIN replicates r ON r.replicateId == q.replicateId;
        '''

        cur = self.conn.cursor()
        cur.execute(query)
        db_proteins = {(row[0], row[1]) for row in cur.fetchall()}

        self.assertEqual(proteins, db_proteins)


class TestMultiProjectStepped(unittest.TestCase):
    PROJECT_1 = 'Strap'
    PROJECT_2 = 'Sp3'
    DATA_DIR = f'{TEST_DIR}/data'

    def test_peptideToProtein(self):
        '''
        Test that no new peptide to protein mappings are added
        when skyline documents have different protein parsimpny settings.
        '''

        # setup variables
        work_dir = f'{TEST_DIR}/work/test_peptideToProtein'
        db_path = f'{work_dir}/data.db3'

        # get report protein mappings
        def get_report_mappings(path):
            df = pd.read_csv(path, sep='\t')
            ret = set()
            for row in df[['ProteinName', 'ModifiedSequence']].drop_duplicates().itertuples():
                ret.add((row.ProteinName, row.ModifiedSequence))
            return ret
        proj_1_map = get_report_mappings(f'{self.DATA_DIR}/skyline_reports/{self.PROJECT_1}'
                                          '_by_protein_precursor_quality.tsv')
        proj_2_map = get_report_mappings(f'{self.DATA_DIR}/invalid_reports/{self.PROJECT_2}'
                                          '_by_protein_not_minimal_precursor_quality.tsv')

        # add PROJECT_1
        first_result = setup_functions.setup_single_db(self.DATA_DIR, work_dir, self.PROJECT_1,
                                                       clear_dir=True, output_prefix='add_project_1')

        # get initial peptide to protein mappings for Strap project
        self.assertEqual(first_result.returncode, 0)
        self.assertTrue(os.path.isfile(db_path))
        conn = sqlite3.connect(db_path)
        cur = conn.cursor()
        query = '''SELECT prot.name, p.modifiedSequence FROM precursors p
                   LEFT JOIN peptideToProtein ptp ON ptp.peptideId == p.peptideId
                   LEFT JOIN proteins prot ON ptp.proteinId == prot.proteinId
                   LEFT JOIN replicates r ON r.replicateId == p.replicateId
                   WHERE r.project == ?; '''
        cur.execute(query, (self.PROJECT_1,))
        db_proj_1_map = {(x[0], x[1]) for x in cur.fetchall()}
        self.assertTrue(db_proj_1_map == proj_1_map)

        command = ['parse_data', '--overwriteMode=append', f'--projectName={self.PROJECT_2}',
                   '-m', f'{self.DATA_DIR}/metadata/{self.PROJECT_2}_metadata.tsv',
                   f'{self.DATA_DIR}/skyline_reports/{self.PROJECT_2}_replicate_quality.tsv',
                   f'{self.DATA_DIR}/invalid_reports/{self.PROJECT_2}_by_protein_not_minimal_precursor_quality.tsv']
        second_result = setup_functions.run_command(command, work_dir, prefix='add_project_2')
        self.assertEqual(second_result.returncode, 0)

        cur = conn.cursor()
        cur.execute(query, (self.PROJECT_1,))
        new_db_proj_1_map = {(x[0], x[1]) for x in cur.fetchall()}
        self.assertEqual(db_proj_1_map, new_db_proj_1_map)
        self.assertEqual(proj_1_map, new_db_proj_1_map)

        cur = conn.cursor()
        cur.execute(query, (self.PROJECT_2,))
        db_proj_2_map = {(x[0], x[1]) for x in cur.fetchall()}
        self.assertEqual(proj_2_map, db_proj_2_map)

        conn.close()


    def test_group_by_gene(self):
        ''' Test --groupBy=gene option '''

        # setup variables
        work_dir = f'{TEST_DIR}/work/test_group_by_gene'
        db_path = f'{work_dir}/data.db3'
        gene_id_path = f'{self.DATA_DIR}/metadata/prhuman2gene_2023_05_24_subset.csv'

        # read protein genes in precursor reports
        gene_groups = dict()
        for project in [self.PROJECT_1, self.PROJECT_2]:
            report_file = f'{self.DATA_DIR}/skyline_reports/{project}_by_gene_precursor_quality.tsv'
            gene_groups[project] = []
            for row in pd.read_csv(report_file, sep='\t').itertuples():
                if pd.isna(row.ProteinGene):
                    gene_groups[project].append(row.ProteinName)
                else:
                    gene_groups[project].append(row.ProteinGene)

            gene_groups[project] = Counter(gene_groups[project])

        # add PROJECT_1
        first_result = setup_functions.setup_single_db(self.DATA_DIR, work_dir, self.PROJECT_1,
                                                       output_prefix='add_project_1',
                                                       clear_dir=True, group_by_gene=True)
        self.assertEqual(first_result.returncode, 0)

        # initialize db connection
        self.assertTrue(os.path.isfile(db_path))
        conn = sqlite3.connect(db_path)
        self.assertEqual('gene', get_meta_value(conn, 'group_precursors_by'))

        # Test adding PROJECT_2 with by_protein grouping fails
        second_result_by_protein = setup_functions.setup_single_db(self.DATA_DIR, work_dir, self.PROJECT_2,
                                                                   output_prefix='add_project_2_by_protein',
                                                                   overwrite_mode='append',
                                                                   clear_dir=False, group_by_gene=False)
        self.assertEqual(second_result_by_protein.returncode, 1)

        # add PROJECT_2
        second_result = setup_functions.setup_single_db(self.DATA_DIR, work_dir, self.PROJECT_2,
                                                        output_prefix='add_project_2',
                                                        overwrite_mode='append',
                                                        clear_dir=False, group_by_gene=True)
        self.assertEqual(second_result.returncode, 0)

        # check that gene groups in reports match groups in database
        for project in [self.PROJECT_1, self.PROJECT_2]:
            cur = conn.cursor()
            query = '''SELECT prot.name FROM precursors p
                       LEFT JOIN peptideToProtein ptp ON ptp.peptideId == p.peptideId
                       LEFT JOIN proteins prot ON ptp.proteinId == prot.proteinId
                       LEFT JOIN replicates r ON r.replicateId == p.replicateId
                       WHERE r.project == ?; '''
            cur.execute(query, (project,))
            db_gene_groups = Counter([x[0] for x in cur.fetchall()])
            self.assertDictEqual(gene_groups[project], db_gene_groups)

        # make sure make_gene_matrix fails if grouped by protein
        conn = update_meta_value(conn, 'group_precursors_by', 'protein')
        bad_matrix_result = setup_functions.run_command(['make_gene_matrix', gene_id_path, db_path],
                                                        work_dir, prefix='failed_matrix')
        self.assertEqual(bad_matrix_result.returncode, 1)
        self.assertTrue('Precursors in database must be grouped by gene!' in bad_matrix_result.stderr)
        conn = update_meta_value(conn, 'group_precursors_by', 'gene')


class TestDuplicatePrecursorsOption(unittest.TestCase):
    TEST_PROJECT = 'Sp3'
    WORK_DIR = f'{TEST_DIR}/work/test_duplicate_precursors'
    DB_PATH = f'{WORK_DIR}/data.db3'
    PRECURSOR_REPPRT=f'{TEST_DIR}/data/invalid_reports/Sp3_by_protein_duplicate_precursor_quality.tsv'

    @classmethod
    def setUpClass(cls):
        setup_functions.make_work_dir(cls.WORK_DIR, clear_dir=True)

    @staticmethod
    def setup_command(project, option, precursor_rep):
        db_name = f'test_{option}_option.db3'
        command = ['parse_data',
                   f'--projectName={project}', '--overwriteMode=overwrite',
                   '-m', f'{TEST_DIR}/data/metadata/{project}_metadata.tsv',
                   '-o', db_name,
                   '--duplicatePrecursors', option,
                   f'{TEST_DIR}/data/skyline_reports/{project}_replicate_quality.tsv',
                   f'{TEST_DIR}/data/invalid_reports/{project}_by_protein_{precursor_rep}_precursor_quality.tsv']
        return command, db_name


    @staticmethod
    def read_precursor_report(path, use_user_set_total):
        df = pd.read_csv(path, sep='\t')

        data = dict()
        for (rep, seq, charge), group in df.groupby(['ReplicateName',
                                                     'ModifiedSequence',
                                                     'PrecursorCharge']):
            precursor = f'{seq}_{charge}'
            if precursor not in data:
                data[precursor] = dict()

            areas = group['TotalAreaFragment'].to_list()
            if all(areas[0] == area for area in areas):
                data[precursor][rep] = areas[0]
                continue

            for row in group.itertuples():
                if not use_user_set_total:
                    assert rep not in data[precursor]
                    data[precursor][rep] = row.TotalAreaFragment
                    break

                if row.UserSetTotal:
                    assert rep not in data[precursor]
                    data[precursor][rep] = row.TotalAreaFragment
                    break

        return data


    @staticmethod
    def read_db_precursos(conn):
        query = ''' SELECT r.replicate, p.modifiedSequence, p.precursorCharge, p.totalAreaFragment
                    FROM precursors p
                    LEFT JOIN replicates r ON r.replicateId == p.replicateId; '''
        df = pd.read_sql(query, conn)
        data = dict()
        for row in df.itertuples():
            precursor = f'{row.modifiedSequence}_{row.precursorCharge}'
            if precursor not in data:
                data[precursor] = {}
            data[precursor][row.replicate] = row.totalAreaFragment

        return data


    def test_error_option_fails(self):
        command, db_name = self.setup_command(self.TEST_PROJECT, 'e', 'duplicate')
        result = setup_functions.run_command(command, self.WORK_DIR)
        self.assertEqual(1, result.returncode)


    def test_invalid_no_user_set_report_fails(self):
        command, db_name = self.setup_command(self.TEST_PROJECT, 'm', 'invalid_no_user_set')
        result = setup_functions.run_command(command, self.WORK_DIR)
        self.assertEqual(1, result.returncode)
        self.assertTrue(re.search(r'ERROR: [0-9]+ precursor groups have no user set peak boundaries!',
                                  result.stderr))


    def test_invalid_other_diff_report_fails(self):
        command, db_name = self.setup_command(self.TEST_PROJECT, 'm', 'invalid_other_diff')
        result = setup_functions.run_command(command, self.WORK_DIR)
        self.assertEqual(1, result.returncode)
        self.assertTrue(re.search(r'ERROR: There are [0-9]+ non-unique precursor areas!',
                                  result.stderr))


    def test_use_user_set_total_option(self):
        command, db_name = self.setup_command(self.TEST_PROJECT, 'm', 'duplicate')
        result = setup_functions.run_command(command, self.WORK_DIR)
        self.assertEqual(0, result.returncode)
        self.assertTrue('There are 9 non-unique precursors!' in result.stderr)
        self.assertTrue('After selecting precursors with user set peak boundaries, ' in result.stderr)

        data = self.read_precursor_report(self.PRECURSOR_REPPRT, True)

        # read db precursors
        self.assertTrue(os.path.isfile(f'{self.WORK_DIR}/{db_name}'))
        conn = sqlite3.connect(f'{self.WORK_DIR}/{db_name}')
        db_data = self.read_db_precursos(conn)
        conn.close()

        self.assertEqual(len(data), len(db_data))
        for precursor in data:
            self.assertEqual(len(data[precursor]), len(db_data[precursor]))
            for rep in data[precursor]:
                self.assertAlmostEqual(data[precursor][rep], db_data[precursor][rep])


    def test_first_option(self):
        command, db_name = self.setup_command(self.TEST_PROJECT, 'f', 'duplicate')
        result = setup_functions.run_command(command, self.WORK_DIR)
        self.assertEqual(0, result.returncode)
        self.assertTrue('There are 9 non-unique precursors!' in result.stderr)

        data = self.read_precursor_report(self.PRECURSOR_REPPRT, False)

        # read db precursors
        self.assertTrue(os.path.isfile(f'{self.WORK_DIR}/{db_name}'))
        conn = sqlite3.connect(f'{self.WORK_DIR}/{db_name}')
        db_data = self.read_db_precursos(conn)
        conn.close()

        self.assertEqual(len(data), len(db_data))
        for precursor in data:
            self.assertEqual(len(data[precursor]), len(db_data[precursor]))
            for rep in data[precursor]:
                self.assertAlmostEqual(data[precursor][rep], db_data[precursor][rep])


if __name__ == '__main__':
    unittest.main()
