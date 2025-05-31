
import unittest
from unittest import mock
import os
import pandas as pd

from DIA_QC_report import validate_pipeline_params
from DIA_QC_report.submodules.panorama import url_exists, have_internet

import setup_functions


class TestGenerateSchemaUrl(unittest.TestCase):
    def test_branch(self):
        if not have_internet():
            self.skipTest("No internet connection available for testing.")

        url = validate_pipeline_params.generate_schema_url(
            repo='ajmaurais/nf-skyline-dia-ms',
            revision='main', filename = 'main.nf'
        )
        self.assertTrue(url_exists(url))


    def test_tag(self):
        if not have_internet():
            self.skipTest("No internet connection available for testing.")

        url = validate_pipeline_params.generate_schema_url(
            repo='mriffle/nf-skyline-dia-ms',
            revision='v2.6.1', filename = 'main.nf'
        )
        self.assertTrue(url_exists(url))


    def test_commit_hash(self):
        if not have_internet():
            self.skipTest("No internet connection available for testing.")

        url = validate_pipeline_params.generate_schema_url(
            repo='ajmaurais/nf-skyline-dia-ms',
            revision='45c69b6bdb8111fb14c140a1275a1a11df47a69d',
            filename = 'main.nf'
        )
        self.assertTrue(url_exists(url))


class TestParseInputFiles(unittest.TestCase):
    def setUp(self):
        self.data_dir = f'{setup_functions.TEST_DIR}/data/'
        self.local_file_dir = f'{self.data_dir}/validate_pipeline_params/fake_local_ms_files/'

        self.projects = ('Strap', 'Sp3')
        self.project_replicates = {}
        for project in self.projects:
            df = pd.read_csv(f'{self.data_dir}/skyline_reports/{project}_replicate_quality.tsv', sep='\t')
            self.project_replicates[project] = df['Replicate'].tolist()


    def test_flat_local_mzml(self):
        test_project = 'Sp3'
        multi_batch_mode, files = validate_pipeline_params.parse_input_files(
            input_files=f'{self.local_file_dir}/{test_project}',
            file_glob='*.mzML'
        )

        self.assertFalse(multi_batch_mode)
        self.assertEqual(len(files), len(self.project_replicates[test_project]))
        self.assertTrue(all(os.path.splitext(file)[1] == '.mzML' for file in files))
        self.assertEqual(
            set(os.path.splitext(file)[0] for file in files),
            set(self.project_replicates[test_project])
        )


    def test_list_local_mzml(self):
        multi_batch_mode, files = validate_pipeline_params.parse_input_files(
            input_files=[f'{self.local_file_dir}/{project}' for project in self.projects],
            file_glob='*.mzML'
        )

        target_reps = set(self.project_replicates['Strap'] + self.project_replicates['Sp3'])
        self.assertFalse(multi_batch_mode)
        self.assertEqual(len(files), len(target_reps))
        self.assertTrue(all(os.path.splitext(file)[1] == '.mzML' for file in files))
        self.assertEqual(
            set(os.path.splitext(file)[0] for file in files),
            set(target_reps)
        )


    def test_multi_batch_local_mzml(self):
        multi_batch_mode, files = validate_pipeline_params.parse_input_files(
            input_files={project: f'{self.local_file_dir}/{project}' for project in self.projects},
            file_glob='*.mzML'
        )

        self.assertIsInstance(files, dict)
        self.assertTrue(multi_batch_mode)
        for project, reps in files.items():
            self.assertTrue(all(os.path.splitext(file)[1] == '.mzML' for file in reps))
            self.assertIn(project, self.projects)
            self.assertEqual(
                set(os.path.splitext(file)[0] for file in reps),
                set(self.project_replicates[project])
            )