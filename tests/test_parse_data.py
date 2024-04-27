
import unittest
from unittest import mock
import sqlite3
import os
import json

import setup_functions

from DIA_QC_report.submodules.metadata import Dtype
from DIA_QC_report.submodules.dia_db_utils import is_normalized

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

        # test command was sucessful
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

        # test command was sucessful
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

        # test command was sucessful
        self.assertEqual(0, result.returncode)

        # Log entry specific for skyline annotations csv
        self.assertTrue('Found Skyline annotations csv' in result.stderr.decode('utf-8'))

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

        # test command was sucessful
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

        # test command was sucessful
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


if __name__ == '__main__':
    unittest.main()
