
import unittest
from unittest import mock
import os
import re
from functools import partial
from shutil import which
import random
import json

import pandas as pd

from DIA_QC_report import validate_pipeline_params
from DIA_QC_report.submodules.panorama import url_exists, have_internet
from DIA_QC_report.submodules.panorama import PANORA_PUBLIC_KEY
from DIA_QC_report.submodules.read_metadata import Metadata

import setup_functions
from test_panorama_functions import PUBLIC_FILE

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
        result = validate_pipeline_params.merge_params(lhs, rhs)

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


class TestGetConfigPath(unittest.TestCase):
    def setUp(self):
        self.rhs = {
            "dir":       {"b1": "/path/b1", "b2": "/path/b2"},
            "qc_report": {"color_vars": ["a", "b", "c"]},
            "fasta": "db.fasta",
            "search_engine": "diann",
        }

    def test_single(self):
        target = 'diann'
        result = validate_pipeline_params._get_config_path(self.rhs, ['search_engine'])
        self.assertEqual(result, target)


    def test_double(self):
        target = '/path/b1'
        result = validate_pipeline_params._get_config_path(
            self.rhs, ['dir', 'b1']
        )
        self.assertEqual(result, target)


    def test_double_missing(self):
        target = None
        result = validate_pipeline_params._get_config_path(
            self.rhs, ['dir', 'b3']
        )
        self.assertEqual(result, target)
        result = validate_pipeline_params._get_config_path(
            self.rhs, ['not_a_key', 'b1']
        )
        self.assertEqual(result, target)


class TestGenerateSchemaUrl(unittest.TestCase):
    def test_branch(self):
        if not have_internet():
            self.skipTest("No internet connection available for testing.")

        url = validate_pipeline_params.generate_git_url(
            repo='ajmaurais/nf-skyline-dia-ms',
            revision='main', filename = 'main.nf'
        )
        self.assertTrue(url_exists(url))


    def test_tag(self):
        if not have_internet():
            self.skipTest("No internet connection available for testing.")

        url = validate_pipeline_params.generate_git_url(
            repo='mriffle/nf-skyline-dia-ms',
            revision='v2.6.1', filename = 'main.nf'
        )
        self.assertTrue(url_exists(url))


    def test_commit_hash(self):
        if not have_internet():
            self.skipTest("No internet connection available for testing.")

        url = validate_pipeline_params.generate_git_url(
            repo='ajmaurais/nf-skyline-dia-ms',
            revision='45c69b6bdb8111fb14c140a1275a1a11df47a69d',
            filename = 'main.nf'
        )
        self.assertTrue(url_exists(url))


class TestGetInputFileText(unittest.TestCase):
    def test_get_input_file_text(self):
        if not have_internet():
            self.skipTest("No internet connection available for testing.")

        text = validate_pipeline_params.get_input_file_text(PUBLIC_FILE, api_key=PANORA_PUBLIC_KEY)
        with open(f'{setup_functions.TEST_DIR}/data/validate_pipeline_params/panorama/Clean-SMTG-B1-1410-Oct2022-3qtrans-Meta.csv', 'r') as f:
            target_text = f.read()

        self.assertEqual(text.splitlines(), target_text.splitlines())


class TestGetApiKeyFromNextflow(unittest.TestCase, setup_functions.AbstractTestsBase):
    def _mock_process_run(self, return_value):
        mock = unittest.mock.MagicMock()
        mock.returncode = 0
        mock.stdout = return_value
        return mock


    def test_get_api_key_from_nextflow_secrets(self):
        with mock.patch('DIA_QC_report.validate_pipeline_params.subprocess.run',
                        return_value=self._mock_process_run(PANORA_PUBLIC_KEY)) as mock_run, \
             mock.patch('DIA_QC_report.validate_pipeline_params.which',
                        return_value='nextflow') as mock_which:
            api_key = validate_pipeline_params.get_api_key_from_nextflow_secrets()

        self.assertEqual(api_key, PANORA_PUBLIC_KEY)
        mock_run.assert_called_once_with(
            ['nextflow', 'secrets', 'get', 'PANORAMA_API_KEY'],
            capture_output=True, text=True
        )


    def test_no_key_in_secrets_manager(self):
        with self.assertLogs(validate_pipeline_params.LOGGER, 'ERROR') as cm:
            with mock.patch('DIA_QC_report.validate_pipeline_params.subprocess.run',
                            return_value=self._mock_process_run('null')), \
                mock.patch('DIA_QC_report.validate_pipeline_params.which', return_value='nextflow'):
                 api_key = validate_pipeline_params.get_api_key_from_nextflow_secrets()

        self.assertIsNone(api_key)
        self.assertInLog("No Panorama API key found in Nextflow secrets.", cm)



class TestParseInputFiles(unittest.TestCase):
    def setUp(self):
        self.data_dir = f'{setup_functions.TEST_DIR}/data/'
        self.local_file_dir = f'{self.data_dir}/validate_pipeline_params/mock_local_ms_files/'

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
                        side_effect=partial(self._target_files_from_url, ext=ext)):
            self.do_test(
                input_files, target_files,
                file_regex=r'.+?(-8mz-ovlp-400to1000-|Small_Assay).+?\.%s$' % ext,
            )


    def test_multi_batch_invalid_directory(self):
        ext = 'raw'
        input_files = {
            'all_remote': [PUBLIC_URL, STRAP_URL],
            'all_local': [f'{self.local_file_dir}/Sp3', f'{self.local_file_dir}/Strap'],
            'mixed': [f'{self.local_file_dir}/Strap', PUBLIC_URL, SP3_URL]
        }
        target_files = {
            'all_remote': self._target_files('Strap', ext),
            'all_local': self._target_files('Strap', ext) + self._target_files('Sp3', ext),
            'mixed': self._target_files('Sp3', ext) + self._target_files('Strap', ext)
        }
        with mock.patch('DIA_QC_report.validate_pipeline_params.list_panorama_files',
                        side_effect=partial(self._target_files_from_url, ext=ext)):
            with self.assertLogs(validate_pipeline_params.LOGGER, 'WARNING') as cm:
                use_multi_batch_mode, files = validate_pipeline_params.parse_input_files(
                    input_files, file_glob=f'*-8mz-ovlp-400to1000-*.{ext}', strict=False
                )

        n_files = sum(len(files) for files in target_files.values())
        self.assertTrue(any(
            f"No MS files could be found in 2 directories. Continuing with {n_files} available files." in m
            for m in cm.output
        ), cm.output)

        self.assertTrue(use_multi_batch_mode)
        for project in target_files:
            self.assertIn(project, files)
            self.assertEqual(set(target_files[project]), set(files[project]),
                             f"Files for project '{project}' do not match expected files.")


    def test_multi_batch_invalid_directory_strict(self):
        ext = 'raw'
        input_files = {
            'all_remote': [PUBLIC_URL, STRAP_URL],
            'all_local': [f'{self.local_file_dir}/Sp3', f'{self.local_file_dir}/Strap'],
            'mixed': [f'{self.local_file_dir}/Strap', PUBLIC_URL, SP3_URL]
        }

        with mock.patch('DIA_QC_report.validate_pipeline_params.list_panorama_files',
                        side_effect=partial(self._target_files_from_url, ext=ext)):
            with self.assertRaises(SystemExit):
                with self.assertLogs(validate_pipeline_params.LOGGER, 'ERROR') as cm:
                    validate_pipeline_params.parse_input_files(
                        input_files, file_glob=f'*-8mz-ovlp-400to1000-*.{ext}', strict=True,
                    )


    def test_list_invalid_directory(self):
        ext = 'raw'
        input_files = [PUBLIC_URL, STRAP_URL]
        target_files = self._target_files('Strap', ext)

        with mock.patch('DIA_QC_report.validate_pipeline_params.list_panorama_files',
                        side_effect=partial(self._target_files_from_url, ext=ext)):
            with self.assertLogs(validate_pipeline_params.LOGGER, 'WARNING') as cm:
                use_batch_mode, files = validate_pipeline_params.parse_input_files(
                    input_files, file_glob=f'*-8mz-ovlp-400to1000-*.{ext}', strict=False
                )

        self.assertTrue(any(
            f"No MS files could be found in 1 directory. Continuing with {len(target_files)} available files." in m
            for m in cm.output
        ), cm.output)

        self.assertFalse(use_batch_mode)
        self.assertEqual(set(target_files), set(files))
        self.assertEqual(len(files), len(target_files))


    def test_flat_invalid_directory_strict(self):
        ext = 'raw'
        input_files = [PUBLIC_URL, STRAP_URL]

        with mock.patch('DIA_QC_report.validate_pipeline_params.list_panorama_files',
                        side_effect=partial(self._target_files_from_url, ext=ext)):
            with self.assertRaises(SystemExit):
                with self.assertLogs(validate_pipeline_params.LOGGER, 'WARNING') as cm:
                    validate_pipeline_params.parse_input_files(
                        input_files, file_glob=f'*-8mz-ovlp-400to1000-*.{ext}', strict=True
                    )


class ValidateMetadata(unittest.TestCase, setup_functions.AbstractTestsBase):
    def setUp(self):
        self.data_dir = f'{setup_functions.TEST_DIR}/data'

        self.projects = ('Strap', 'Sp3')
        self.project_replicates = {}
        for project in self.projects:
            df = pd.read_csv(f'{self.data_dir}/skyline_reports/{project}_replicate_quality.tsv', sep='\t')
            self.project_replicates[project] = df['FileName'].tolist()

        metadata_reader = Metadata()
        metadata_file = f'{self.data_dir}/metadata/Sp3_Strap_combined_metadata.tsv'
        if not metadata_reader.read(metadata_file):
            raise FileNotFoundError(f"Metadata file '{metadata_file}' not found or could not be read.")

        self.combined_metadata = metadata_reader.df
        self.metadata_types = metadata_reader.types


class TestValidateMetadataParams(ValidateMetadata):
    def test_missing_metadata_file(self):
        color_vars=['cellLine', 'experiment']

        with self.assertLogs(validate_pipeline_params.LOGGER, 'ERROR') as cm:
            success = validate_pipeline_params.validate_metadata(
                self.project_replicates, None, None, color_vars=color_vars
            )

        self.assertFalse(success)
        self.assertInLog('Replicate metadata is None, but metadata parameters are specified.', cm)
        for var in color_vars:
            self.assertInLog(f"Without a metadata file, parameter color_vars = '{var}' will cause an error.", cm)


    def test_all_null_metadata_var(self):
        reader = Metadata()
        self.assertTrue(reader.read(f'{self.data_dir}/metadata/Strap_missing_multi_var_metadata.json'),
                        'Failed to read metadata file')
        reps = {None: self.project_replicates['Strap']}

        with self.assertLogs(validate_pipeline_params.LOGGER, 'WARN') as cm:
            success = validate_pipeline_params.validate_metadata(
                reps, reader.df, reader.types, color_vars=['na_var']
            )
        self.assertTrue(success)
        self.assertInLog("All values for metadata variable 'na_var' from 'color_vars' parameter are NULL.", cm)


    def test_missing_color_var(self):
        with self.assertLogs(validate_pipeline_params.LOGGER, 'ERROR') as cm:
            success = validate_pipeline_params.validate_metadata(
                self.project_replicates, self.combined_metadata, self.metadata_types,
                color_vars=['missing_color_var']
            )

        self.assertFalse(success)
        self.assertInLog("Metadata variable 'missing_color_var' from 'color_vars' parameter not found in metadata.", cm)


    def test_missing_batch_var(self):
        with self.assertLogs(validate_pipeline_params.LOGGER, 'ERROR') as cm:
            success = validate_pipeline_params.validate_metadata(
                self.project_replicates, self.combined_metadata, self.metadata_types,
                batch1='notABatch'
            )

        self.assertFalse(success)
        self.assertInLog("Metadata variable 'notABatch' from 'batch1' parameter not found in metadata.", cm)


    def test_only_control_key(self):
        with self.assertLogs(validate_pipeline_params.LOGGER, 'ERROR') as cm:
            success = validate_pipeline_params.validate_metadata(
                self.project_replicates, self.combined_metadata, self.metadata_types,
                control_key='cellLine'
            )
        self.assertFalse(success)
        self.assertInLog("Both 'control_key' and 'color_vars' must be specified or both omitted.", cm)


    def test_missing_control_key(self):
        with self.assertLogs(validate_pipeline_params.LOGGER, 'ERROR') as cm:
            success = validate_pipeline_params.validate_metadata(
                self.project_replicates, self.combined_metadata, self.metadata_types,
                control_key='notAVar', control_values=['T47D', 'HeLa']
            )
        self.assertFalse(success)
        self.assertInLog("Metadata variable 'notAVar' from 'control_key' parameter not found in metadata.", cm)


    def test_missing_control_values(self):
        with self.assertLogs(validate_pipeline_params.LOGGER, 'ERROR') as cm:
            success = validate_pipeline_params.validate_metadata(
                self.project_replicates, self.combined_metadata, self.metadata_types,
                control_key='cellLine', control_values=['notACellLine']
            )
        self.assertFalse(success)
        self.assertInLog("Control value 'notACellLine' for key 'cellLine' not found in metadata.", cm)


class TestValidateMetadataReps(ValidateMetadata):
    def test_duplicate_replicate(self):
        random.seed(40)
        duplicate_rep = random.choice(self.project_replicates['Sp3'])
        self.project_replicates['Strap'].append(duplicate_rep)

        with self.assertLogs(validate_pipeline_params.LOGGER, 'ERROR') as cm:
            success = validate_pipeline_params.validate_metadata(
                self.project_replicates, self.combined_metadata, self.metadata_types
            )
        self.assertInLog(f"Replicate '{duplicate_rep}' is duplicated 1 time in quant_spectra_dir.", cm)
        self.assertFalse(success)


    def test_duplicate_replicate_flat(self):
        random.seed(12)
        duplicate_rep = random.choice(self.project_replicates['Strap'])
        project_reps_flat = {None: self.project_replicates['Strap']}

        n_duplicates = random.randint(1, 5)
        for _ in range(n_duplicates):
            project_reps_flat[None].append(duplicate_rep)

        with self.assertLogs(validate_pipeline_params.LOGGER, 'ERROR') as cm:
            success = validate_pipeline_params.validate_metadata(
                project_reps_flat, self.combined_metadata, self.metadata_types
            )
        self.assertInLog(f"Replicate '{duplicate_rep}' is duplicated {n_duplicates} times in quant_spectra_dir.", cm)
        self.assertFalse(success)


    def test_extra_metadata_reps(self):
        project_reps_flat = {None: self.project_replicates['Strap']}

        with self.assertLogs(validate_pipeline_params.LOGGER, 'WARN') as cm:
            success = validate_pipeline_params.validate_metadata(
                project_reps_flat, self.combined_metadata, self.metadata_types
            )
        self.assertInLog("There are 20 replicates in the metadata that are not in quant_spectra_dir", cm)
        self.assertTrue(success)


    def test_missing_metadata_rep(self):
        random.seed(7)
        remove_reps = set(random.sample(self.combined_metadata['Replicate'].drop_duplicates().to_list(), 2))
        metadata_df = self.combined_metadata.loc[~self.combined_metadata['Replicate'].isin(remove_reps)]
        metadata_df = metadata_df.reset_index(drop=True)

        # test with strict=True
        with self.assertLogs(validate_pipeline_params.LOGGER, 'ERROR') as cm:
            success = validate_pipeline_params.validate_metadata(
                self.project_replicates, metadata_df, self.metadata_types,
                strict=True
            )
        for rep in remove_reps:
            self.assertInLog(f":Replicate '{rep}' in quant_spectra_dir is not present in the metadata.", cm)
        self.assertFalse(success)

        # test with strict=False
        with self.assertLogs(validate_pipeline_params.LOGGER, 'WARNING') as cm:
            success = validate_pipeline_params.validate_metadata(
                self.project_replicates, metadata_df, self.metadata_types,
                strict=False
            )
        for rep in remove_reps:
            self.assertInLog(f":Replicate '{rep}' in quant_spectra_dir is not present in the metadata.", cm)
        self.assertTrue(success)


    def test_duplicate_metadata_value(self):
        random.seed(12)
        i = random.randint(0, len(self.combined_metadata.index) - 1)
        new_row = self.combined_metadata.loc[i].to_list()
        new_row[2] = 'NewValue'
        metadata_df = self.combined_metadata
        metadata_df.loc[len(metadata_df.index)] = new_row

        # test with strict=True
        with self.assertLogs(validate_pipeline_params.LOGGER, 'ERROR') as cm:
            success = validate_pipeline_params.validate_metadata(
                self.project_replicates, metadata_df, self.metadata_types,
                strict=True
            )
        self.assertInLog(f"Replicate '{new_row[0]}' already has metadata key '{new_row[1]}'. Overwriting with value '{new_row[2]}'.", cm)
        self.assertInLog('Duplicate metadata keys found in replicates.', cm)
        self.assertFalse(success)

        # test with strict=False
        with self.assertLogs(validate_pipeline_params.LOGGER, 'WARN') as cm:
            success = validate_pipeline_params.validate_metadata(
                self.project_replicates, metadata_df, self.metadata_types,
                strict=False
            )
        self.assertInLog(f"Replicate '{new_row[0]}' already has metadata key '{new_row[1]}'. Overwriting with value '{new_row[2]}'.", cm)
        self.assertNotInLog('Duplicate metadata keys found in replicates.', cm)
        self.assertTrue(success)


    def test_missing_metadata_value(self):
        random.seed(7)

        color_vars = ['cellLine', 'experiment']
        metadata_df = self.combined_metadata
        remove_index = random.sample(
            metadata_df[metadata_df['annotationKey'].isin(color_vars)].index.to_list(),
            random.randint(1, min(len(metadata_df.index), 5))
        )
        removed_reps = metadata_df.loc[remove_index].to_dict()
        metadata_df = metadata_df.drop(remove_index).reset_index(drop=True)

        # test with strict=True
        with self.assertLogs(validate_pipeline_params.LOGGER, 'ERROR') as cm:
            success = validate_pipeline_params.validate_metadata(
                self.project_replicates, metadata_df, self.metadata_types,
                color_vars=color_vars, strict=True
            )
        self.assertFalse(success)
        for i in remove_index:
            self.assertInLog(
                "Replicate '{}' is missing metadata key '{}' from 'color_vars' parameter.".format(
                    removed_reps['Replicate'][i],
                    removed_reps['annotationKey'][i]),
                cm
            )

        # test with strict=False
        with self.assertLogs(validate_pipeline_params.LOGGER, 'WARNING') as cm:
            success = validate_pipeline_params.validate_metadata(
                self.project_replicates, metadata_df, self.metadata_types,
                color_vars=color_vars, strict=False
            )
        self.assertTrue(success)
        for i in remove_index:
            self.assertInLog(
                "Replicate '{}' is missing metadata key '{}' from 'color_vars' parameter.".format(
                    removed_reps['Replicate'][i],
                    removed_reps['annotationKey'][i]),
                cm
            )


    def test_na_metadata_value(self):
        reader = Metadata()
        self.assertTrue(reader.read(f'{self.data_dir}/metadata/Strap_missing_multi_var_metadata.json'),
                        "Failed to read metadata file for testing.")
        reps = {None: self.project_replicates['Strap']}
        color_vars=['int_var', 'float_var', 'bool_var']
        missing_counts = {}
        for var in color_vars:
            missing_counts[var] = sum(reader.df[reader.df['annotationKey'] == var]['annotationValue'].apply(
                lambda x: pd.isna(x) or x is None or x == ''
            ))

        with self.assertLogs(validate_pipeline_params.LOGGER, 'WARNING') as cm:
            success = validate_pipeline_params.validate_metadata(
                reps, reader.df, reader.types, color_vars=color_vars
            )
        self.assertTrue(success)

        for var, n_missing in missing_counts.items():
            if n_missing > 0:
                self.assertInLog(
                    f"{str(reader.types[var])} variable '{var}' is missing in {n_missing} of {len(reps[None])} replicates.",
                    cm
                )
            else:
                self.assertNotInLog(f"{str(reader.types[var])} variable '{var}' is missing in", cm)


class TestWriteReorts(ValidateMetadata):
    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_write_parameter_validation_reports'
        setup_functions.make_work_dir(cls.work_dir, clear_dir=True)


    def test_write_json_metadata_report(self):
        with self.assertLogs(validate_pipeline_params.LOGGER, 'INFO') as cm:
            success = validate_pipeline_params.validate_metadata(
                self.project_replicates, self.combined_metadata, self.metadata_types,
                color_vars=['cellLine', 'experiment', 'NCI7std'],
                control_key='cellLine', control_values=['T47D', 'HeLa'],
                write_metadata_report=True, report_dir=self.work_dir, report_format='json'
            )

        self.assertTrue(success, cm.output)
        self.assertInLog(f'Metadata report written to {self.work_dir}', cm)
        self.assertTrue(
            os.path.isfile(f'{self.work_dir}/metadata_validation_report.json'), "Metadata report file was not created."
        )

        with open(f'{self.work_dir}/metadata_validation_report.json', 'r') as inF:
            report = json.load(inF)

        expected_rows = {
            ('cellLine', 'color_vars'): {'type': 'STRING'},
            ('experiment', 'color_vars'): {'type': 'STRING'},
            ('NCI7std', 'color_vars'): {'type': 'BOOL'},
            ('cellLine', 'control_key'): {'type': 'STRING'}
        }
        for row in expected_rows:
            expected_rows[row]['variable'] = row[0]
            expected_rows[row]['parameter'] = row[1]
            expected_rows[row]['missing_in_n_replicates'] = 0
            expected_rows[row]['found_in_n_replicates'] = sum(len(reps) for reps in self.project_replicates.values())

        report_dict = {(row['variable'], row['parameter']): row for row in report}

        for key, row in expected_rows.items():
            self.assertIn(key, report_dict, f"Missing expected row for {key} in report.")
            self.assertDictEqual(report_dict[key], row, f"Row for {key} does not match expected values.")

        self.assertEqual(len(report), 4)


    def test_write_tsv_metadata_report(self):
        with self.assertLogs(validate_pipeline_params.LOGGER, 'INFO') as cm:
            success = validate_pipeline_params.validate_metadata(
                self.project_replicates, self.combined_metadata, self.metadata_types,
                color_vars=['cellLine', 'experiment', 'NCI7std'],
                control_key='cellLine', control_values=['T47D', 'HeLa'],
                write_metadata_report=True, report_dir=self.work_dir, report_format='tsv'
            )

        self.assertTrue(success, cm.output)
        self.assertInLog(f'Metadata report written to {self.work_dir}', cm)
        self.assertTrue(
            os.path.isfile(f'{self.work_dir}/metadata_validation_report.tsv'), "Metadata report file was not created."
        )

        with open(f'{self.work_dir}/metadata_validation_report.tsv', 'r') as inF:
            lines = inF.readlines()

        self.assertEqual(len(lines), 5, "Metadata report file does not have the expected number of lines.")


    def test_write_json_replicate_report(self):
        with self.assertLogs(validate_pipeline_params.LOGGER, 'INFO') as cm:
            success = validate_pipeline_params.validate_metadata(
                self.project_replicates, self.combined_metadata, self.metadata_types,
                color_vars=['cellLine', 'experiment', 'NCI7std'],
                write_replicate_report=True, report_dir=self.work_dir, report_format='json'
            )

        self.assertTrue(success, cm.output)
        self.assertInLog(f'Replicate report written to {self.work_dir}', cm)
        self.assertTrue(
            os.path.isfile(f'{self.work_dir}/replicate_validation_report.json'), "Replicate report file was not created."
        )


    def test_write_tsv_replicate_report(self):
        with self.assertLogs(validate_pipeline_params.LOGGER, 'INFO') as cm:
            success = validate_pipeline_params.validate_metadata(
                self.project_replicates, self.combined_metadata, self.metadata_types,
                color_vars=['cellLine', 'experiment', 'NCI7std'],
                write_replicate_report=True, report_dir=self.work_dir, report_format='tsv'
            )

        self.assertTrue(success, cm.output)
        self.assertInLog(f'Replicate report written to {self.work_dir}', cm)
        self.assertTrue(
            os.path.isfile(f'{self.work_dir}/replicate_validation_report.tsv'), "Replicate report file was not created."
        )

        report = pd.read_csv(f'{self.work_dir}/replicate_validation_report.tsv', sep='\t')
        self.assertEqual(len(report.index), sum(len(reps) for reps in self.project_replicates.values()),
                         "Replicate report does not have the expected number of rows.")


    def test_reserved_replicate_report_header(self):
        replicates = {}
        batch = 'Strap'
        for rep in self.project_replicates[batch]:
            replicates[rep] = validate_pipeline_params.Replicate(rep, batch=batch)
            replicates[rep].metadata['ParameterBatch'] = batch
            replicates[rep].metadata['ParameterBatch_1'] = batch

        with self.assertLogs(validate_pipeline_params.LOGGER, 'WARNING') as cm:
            validate_pipeline_params._write_replicate_report(
                replicates, output_path=f'{self.work_dir}/test_reserved_report_header_report.tsv',
            )

        self.assertInLog("'ParameterBatch' is a reserved column header name. It will be changed to 'ParameterBatch_2' in the report.", cm)


class TestRemoveNoneFromParamDict(unittest.TestCase):
    def test_remove_none_values(self):
        params = { 'param1': 'value1', 'param2': None, 'param3': 'value3', 'param4': None }
        expected = { 'param1': 'value1', 'param3': 'value3' }
        result = validate_pipeline_params.remove_none_from_param_dict(params)
        self.assertDictEqual(result, expected)


    def test_remove_none_values_in_branch(self):
        params   = {'param1': 'value1',
                    'param2': None,
                    'param3': 'value3',
                    'param4': {'n1': 'v', 'n2': None}}
        expected = {'param1': 'value1',
                    'param3': 'value3',
                    'param4': {'n1': 'v'} }
        result = validate_pipeline_params.remove_none_from_param_dict(params)
        self.assertDictEqual(result, expected)


    def test_remove_empty_branch(self):
        params   = {'param1': 'value1',
                    'param2': None,
                    'param3': 'value3',
                    'param4': {'n1': None, 'n2': None} }
        expected = {'param1': 'value1',
                    'param3': 'value3' }
        result = validate_pipeline_params.remove_none_from_param_dict(params)
        self.assertDictEqual(result, expected)


    def test_cascading_branch_removal(self):
        params   = {'param1': 'value1',
                    'param2': None,
                    'param3': 'value3',
                    'param4': {'n1': {'inner': None}}}
        expected = {'param1': 'value1',
                    'param3': 'value3' }
        result = validate_pipeline_params.remove_none_from_param_dict(params)
        self.assertDictEqual(result, expected)


    def test_list_unchanged(self):
        params   = {'param1': 'value1',
                    'param2': None,
                    'param3': 'value3',
                    'param4': ['a', None, 'b']}
        expected = {'param1': 'value1',
                    'param3': 'value3',
                    'param4': ['a', None, 'b']}
        result = validate_pipeline_params.remove_none_from_param_dict(params)
        self.assertDictEqual(result, expected)


    def test_primitive_returns_unchanged(self):
        self.assertEqual(validate_pipeline_params.remove_none_from_param_dict(42), 42)
        self.assertEqual(validate_pipeline_params.remove_none_from_param_dict("string"), "string")
        self.assertEqual(validate_pipeline_params.remove_none_from_param_dict(3.14), 3.14)
        self.assertEqual(validate_pipeline_params.remove_none_from_param_dict(True), True)


class TestValidateConfigFiles(unittest.TestCase, setup_functions.AbstractTestsBase):
    def setUp(self):
        self.data_dir = f'{setup_functions.TEST_DIR}/data/'
        self.local_schema_path = f'{self.data_dir}/validate_pipeline_params/config_files/nextflow_schema.json'
        self.local_config_path = f'{self.data_dir}/validate_pipeline_params/config_files/nextflow.config'


    def test_bad_config_url(self):
        with self.assertLogs(validate_pipeline_params.LOGGER, 'ERROR') as cm:
            successs, data = validate_pipeline_params.validate_config_files(
                'https://example.com/bad_config.config', self.local_schema_path
            )
        self.assertFalse(successs, "Config validation should fail for a bad URL.")
        self.assertInLog('Error downloading pipeline config file', cm)


    def test_valid_config(self):
        config_path = f'{self.data_dir}/validate_pipeline_params/config_files/pdc_pipeline.config'
        success, data = validate_pipeline_params.validate_config_files(
            [self.local_config_path, config_path], self.local_schema_path
        )
        self.assertTrue(success, "Config validation failed for a valid config file.")
