
import unittest
from unittest import mock
import os
from functools import partial
import re
import pandas as pd

from DIA_QC_report import validate_pipeline_params
from DIA_QC_report.submodules.panorama import url_exists, have_internet
from DIA_QC_report.submodules.panorama import PANORA_PUBLIC_KEY

import setup_functions

STRAP_URL = 'https://panoramaweb.org/_webdav/ICPC/NCI-7%20Joint%20Project/NCI-7%20Data%20Harmonization/LFQ-Analyses/USA-UW/%40files/RawFiles/S-Trap/'
SP3_URL = 'https://panoramaweb.org/_webdav/ICPC/NCI-7%20Joint%20Project/NCI-7%20Data%20Harmonization/LFQ-Analyses/USA-UW/%40files/RawFiles/SP3/'
PUBLIC_URL = 'https://panoramaweb.org/_webdav/Panorama%20Public/2024/Thermo%20Fisher%20Research%20and%20Development%20-%202024_Stellar_Instrument_Platform/Label%20Free%20-%20E.%20coli/%40files/RawFiles/ReplicatesSmall/'

class TestAddParams(unittest.TestCase):
    def test_add_params(self):
        lhs = {
            "fasta": None,
            "search_engine": "encyclopedia",
            "skyline": {"skip": False, "doc_name": "final"},
        }
        rhs = {
            "dir":       {"b1": "/path/b1", "b2": "/path/b2"},
            "qc_report": {"color_vars": ["a", "b", "c"]},
            "fasta": "db.fasta",
            "search_engine": "diann",
        }
        target = {
            'dir': {'b1': '/path/b1', 'b2': '/path/b2'},
            'qc_report': {'color_vars': ['a', 'b', 'c']},
            'fasta': 'db.fasta',
            'search_engine': 'diann',
            'skyline': {'skip': False, 'doc_name': 'final'}
        }
        result = validate_pipeline_params.add_params(lhs, rhs)

        self.assertIsInstance(result, dict)
        self.assertEqual(set(result.keys()), set(target.keys()))
        for key in target:
            if isinstance(target[key], dict):
                self.assertIsInstance(result[key], dict)
                self.assertEqual(set(result[key].keys()), set(target[key].keys()))
                for subkey in target[key]:
                    self.assertEqual(result[key][subkey], target[key][subkey])
            else:
                self.assertEqual(result[key], target[key])


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

        self.projects = ['Strap', 'Sp3', 'ReplicatesSmall']
        self.project_replicates = {}
        for project in ['Strap', 'Sp3']:
            df = pd.read_csv(f'{self.data_dir}/skyline_reports/{project}_replicate_quality.tsv', sep='\t')
            self.project_replicates[project] = df['Replicate'].tolist()

        self.project_replicates['ReplicatesSmall'] = [
            '240419_P1_Neo_60SPD_EColiHeLa_Small_Assay_Rep1',
            '240419_P1_Neo_60SPD_EColiHeLa_Small_Assay_Rep2',
            '240419_P1_Neo_60SPD_EColiHeLa_Small_Assay_Rep3',
            '240419_P1_Neo_60SPD_EColiHeLa_Small_Assay_Rep4',
            '240419_P1_Neo_60SPD_EColiHeLa_Small_Assay_Rep5',
            '240419_P1_Neo_60SPD_EColiHeLa_Small_Assay_Rep6',
            '240419_P1_Neo_60SPD_EColiHeLa_Small_Assay_Rep7',
            '240419_P1_Neo_60SPD_EColiHeLa_Small_Assay_Rep8'
        ]


    def do_test(self, input_files, target_files, **kwargs):
        with self.assertNoLogs(validate_pipeline_params.LOGGER, 'ERROR'):
            multi_batch_mode, files = validate_pipeline_params.parse_input_files(
                input_files=input_files,
                **kwargs
            )

        self.assertEqual(multi_batch_mode, isinstance(input_files, dict))

        if multi_batch_mode:
            self.assertIsInstance(files, dict)
            for project, reps in files.items():
                self.assertIn(project, target_files)
                self.assertEqual(set(reps), set(target_files[project]))
                self.assertEqual(len(reps), len(target_files[project]))
        else:
            self.assertIsInstance(files, list)
            self.assertEqual(set(files), set(target_files))
            self.assertEqual(len(files), len(target_files))


    def _target_files(self, project, ext):
        return [f'{rep}.{ext}' for rep in self.project_replicates[project]]


    def _target_files_from_url(self, url, ext, api_key=None):
        if url == SP3_URL:
            return self._target_files('Sp3', ext)
        elif url == STRAP_URL:
            return self._target_files('Strap', ext)
        elif url == PUBLIC_URL:
            return self._target_files('ReplicatesSmall', ext)
        else:
            raise ValueError(f'Unexpected URL: {url}')


    def test_flat_local(self):
        input_files = f'{self.local_file_dir}/Sp3'
        for ext in ['mzML', 'raw']:
            target_files = self._target_files('Sp3', ext)
            self.do_test(input_files, target_files, file_glob=f'*-8mz-ovlp-400to1000-*.{ext}')


    def test_list_local(self):
        input_files = [f'{self.local_file_dir}/{project}' for project in ('Sp3', 'Strap')]
        for ext in ['mzML', 'raw']:
            target_files = [f'{file}.{ext}' for file in self.project_replicates['Strap'] + self.project_replicates['Sp3']]
            self.do_test(input_files, target_files, file_glob=f'*-8mz-ovlp-400to1000-*.{ext}')


    def test_multi_batch_local_mzml(self):
        input_files = {project: f'{self.local_file_dir}/{project}' for project in ('Sp3', 'Strap')}
        for ext in ['mzML', 'raw']:
            target_files = {project: [f'{rep}.{ext}' for rep in replicates]
                            for project, replicates in self.project_replicates.items()}
            self.do_test(input_files, target_files, file_glob=f'*-8mz-ovlp-400to1000-*.{ext}')


    def test_flat_panorama(self):
        input_files = SP3_URL
        target_files = [f'{rep}.raw' for rep in self.project_replicates['Sp3']]
        with mock.patch('DIA_QC_report.validate_pipeline_params.list_panorama_files',
                        return_value=target_files) as mock_list_panorama:
            self.do_test(input_files, target_files, file_glob='*-8mz-ovlp-400to1000-*.raw')

        mock_list_panorama.assert_called_once_with(input_files, api_key=None)


    def test_list_panorama(self):
        ext = 'raw'
        input_files = [SP3_URL, PUBLIC_URL]
        target_files = self._target_files('Sp3', ext=ext) + self._target_files('ReplicatesSmall', ext=ext)
        with mock.patch('DIA_QC_report.validate_pipeline_params.list_panorama_files',
                        side_effect=partial(self._target_files_from_url, ext=ext)) as mock_list_panorama:
            self.do_test(input_files, target_files, file_regex=r'.+?(-8mz-ovlp-400to1000-|Small_Assay).+?\.%s$' % ext)

        mock_list_panorama.assert_has_calls([mock.call(arg, api_key=None) for arg in input_files])


    def test_list_mixed(self):
        ext = 'raw'
        input_files = [SP3_URL, f'{self.local_file_dir}/ReplicatesSmall']
        target_files = self._target_files('Sp3', ext=ext) + self._target_files('ReplicatesSmall', ext=ext)
        with mock.patch('DIA_QC_report.validate_pipeline_params.list_panorama_files',
                        side_effect=partial(self._target_files_from_url, ext=ext)) as mock_list_panorama:
            self.do_test(input_files, target_files, file_regex=r'.+?(-8mz-ovlp-400to1000-|Small_Assay).+?\.%s$' % ext)

        mock_list_panorama.assert_called_once_with(input_files[0], api_key=None)


    def test_multi_batch_panorama(self):
        ext = 'raw'
        input_files = {'Strap': STRAP_URL, 'Sp3': SP3_URL }
        target_files = {
            'Strap': self._target_files('Strap', ext),
            'Sp3': self._target_files('Sp3', ext)
        }
        with mock.patch('DIA_QC_report.validate_pipeline_params.list_panorama_files',
                        side_effect=partial(self._target_files_from_url, ext=ext)) as mock_list_panorama:
            self.do_test(input_files, target_files, file_glob=f'*-8mz-ovlp-400to1000-*.{ext}')

        mock_list_panorama.assert_has_calls([
            mock.call(input_files['Strap'], api_key=None),
            mock.call(input_files['Sp3'], api_key=None)
        ])


    def test_multi_batch_mixed(self):
        ext = 'mzML'
        input_files = {
            'Strap': STRAP_URL,
            'Sp3': f'{self.local_file_dir}/Sp3'
        }
        target_files = {
            'Strap': self._target_files('Strap', ext),
            'Sp3': self._target_files('Sp3', ext),
        }
        with mock.patch('DIA_QC_report.validate_pipeline_params.list_panorama_files',
                        side_effect=partial(self._target_files_from_url, ext=ext)) as mock_list_panorama:
            self.do_test(
                input_files, target_files,
                api_key=PANORA_PUBLIC_KEY,
                file_glob=f'*-8mz-ovlp-400to1000-*.{ext}'
            )

        mock_list_panorama.assert_called_once_with(input_files['Strap'], api_key=PANORA_PUBLIC_KEY)


    def test_multi_batch_mixed_list(self):
        ext = 'raw'
        input_files = {
            'all_remote': [PUBLIC_URL, STRAP_URL],
            'all_local': [f'{self.local_file_dir}/Sp3', f'{self.local_file_dir}/Strap'],
            'mixed': [f'{self.local_file_dir}/Strap', PUBLIC_URL, SP3_URL]
        }
        target_files = {
            'all_remote': self._target_files('Strap', ext) + self._target_files('ReplicatesSmall', ext),
            'all_local': self._target_files('Strap', ext) + self._target_files('Sp3', ext),
            'mixed': self._target_files('Sp3', ext) + self._target_files('Strap', ext) + self._target_files('ReplicatesSmall', ext)
        }
        with mock.patch('DIA_QC_report.validate_pipeline_params.list_panorama_files',
                        side_effect=partial(self._target_files_from_url, ext=ext)) as mock_list_panorama:
            self.do_test(
                input_files, target_files,
                file_regex=r'.+?(-8mz-ovlp-400to1000-|Small_Assay).+?\.%s$' % ext,
            )