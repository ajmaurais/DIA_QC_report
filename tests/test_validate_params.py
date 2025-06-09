
import unittest
from unittest import mock
from types import SimpleNamespace
import os
import json

import pandas as pd

from DIA_QC_report.submodules.read_metadata import Metadata
from DIA_QC_report import validate_pipeline_params as vpp
from DIA_QC_report.submodules.pipeline_config import PipelineConfig
from DIA_QC_report.submodules.panorama import have_internet
from DIA_QC_report.submodules.panorama import PANORAMA_PUBLIC_KEY

import setup_functions

MOCK_PANORAMA = True

STRAP_URL = 'https://panoramaweb.org/_webdav/ICPC/NCI-7%20Joint%20Project/NCI-7%20Data%20Harmonization/LFQ-Analyses/USA-UW/%40files/RawFiles/S-Trap/'
SP3_URL = 'https://panoramaweb.org/_webdav/ICPC/NCI-7%20Joint%20Project/NCI-7%20Data%20Harmonization/LFQ-Analyses/USA-UW/%40files/RawFiles/SP3/'
PUBLIC_URL = 'https://panoramaweb.org/_webdav/Panorama%20Public/2024/Thermo%20Fisher%20Research%20and%20Development%20-%202024_Stellar_Instrument_Platform/Label%20Free%20-%20E.%20coli/%40files/RawFiles/ReplicatesSmall/'

METADATE_CONFIG_TO_ARG_PARAMS = {
    'qc_report.color_vars': 'color_vars', 'qc_report.batch1': 'batch1', 'qc_report.batch2': 'batch2',
    'qc_report.control_key': 'control_key', 'qc_report.control_values': 'control_values'
}

def mock_list_list_panorama_files(url, **kwargs):
    if url == STRAP_URL:
        return os.listdir(setup_functions.TEST_DIR + '/data/validate_pipeline_params/mock_local_ms_files/Strap')
    elif url == SP3_URL:
        return os.listdir(setup_functions.TEST_DIR + '/data/validate_pipeline_params/mock_local_ms_files/Sp3')
    elif url == PUBLIC_URL:
        return os.listdir(setup_functions.TEST_DIR + '/data/validate_pipeline_params/mock_local_ms_files/ReplicatesSmall')
    else:
        raise ValueError(f'Unknown URL: {url}')


def setup_test_config(base_config, output_path=None, add_params=None):
    config = PipelineConfig(file=base_config)

    def add_attribute(root, path, value):
        node = root
        for key in path[:-1]:
            if isinstance(node, SimpleNamespace):
                if not hasattr(node, key):
                    setattr(node, key, SimpleNamespace())
                node = getattr(node, key)
            elif isinstance(node, dict):
                if key not in node:
                    node[key] = SimpleNamespace()
                node = node[key]
            else:
                raise TypeError(f'cannot descend into {type(node).__name__}')

        last = path[-1]
        if isinstance(node, SimpleNamespace):
            setattr(node, last, value)
        elif isinstance(node, dict):
            node[last] = value
        else:
            raise TypeError(f'cannot set value on {type(node).__name__}')

    for key, value in add_params.items():
        path = key.split('.')
        add_attribute(config.params, path, value)

    if output_path is None:
        return config
    config.write(output_path)


class TestValidateSetup(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data_dir = f'{setup_functions.TEST_DIR}/data'
        cls.local_ms_files = f'{cls.data_dir}/validate_pipeline_params/mock_local_ms_files'
        cls.local_db = f'{cls.data_dir}/validate_pipeline_params/mock_local_db'

        cls.projects = ('Strap', 'Sp3')
        cls.project_replicates = {}
        for project in cls.projects:
            df = pd.read_csv(f'{cls.data_dir}/skyline_reports/{project}_replicate_quality.tsv', sep='\t')
            cls.project_replicates[project] = df['Replicate'].tolist()

        metadata_reader = Metadata()
        metadata_file = f'{cls.data_dir}/metadata/Sp3_Strap_combined_metadata.tsv'
        if not metadata_reader.read(metadata_file):
            raise FileNotFoundError(f"Metadata file '{metadata_file}' not found or could not be read.")

        cls.combined_metadata = metadata_reader.df
        cls.metadata_types = metadata_reader.types


    def check_metadata_report(
        self, report_path,
        color_vars=None, control_key=None, control_values=None, batch1=None, batch2=None
    ):
        meta_params = []
        if color_vars:
            meta_params = [(var, 'color_vars') for var in color_vars]
        if control_key:
            meta_params.append((control_key, 'control_key'))
        if batch1:
            meta_params.append((batch1, 'batch1'))
        if batch2:
            meta_params.append((batch2, 'batch2'))

        self.assertTrue(os.path.exists(report_path), f'Metadata report does not exist: {report_path}')
        with open(report_path, 'r') as f:
            report = json.load(f)

        self.assertEqual(len(report), len(meta_params))
        for row in report:
            for key in ['variable', 'parameter', 'type']:
                self.assertIn(key, row, f'Missing key in metadata report row: {key}')
            param_key = (row['variable'], row['parameter'])
            self.assertIn(param_key, meta_params, f'Unexpected metadata parameter: {param_key}')


    def check_replicate_report(self, report_path, projects):
        single_batch = False
        if isinstance(projects, str):
            projects = [projects]
            single_batch = True

        self.assertTrue(os.path.exists(report_path), f'Replicate report does not exist: {report_path}')
        with open(report_path, 'r') as f:
            report = json.load(f)

        report_grouped = {}
        for row in report:
            batch = row.pop('ParameterBatch')
            if batch not in report_grouped:
                report_grouped[batch] = []
            report_grouped[batch].append(row)

        self.assertEqual(len(report_grouped), len(projects), 'Unexpected number of batches in report')
        for project in projects:
            if single_batch:
                self.assertIn(None, report_grouped, f'Report should have only 1 batch')
            else:
                self.assertIn(project, report_grouped, f'Missing batch: {project}')

            self.assertEqual(
                len(report_grouped[None if single_batch else project]), len(self.project_replicates[project]),
                f'Unexpected number of replicates for project {project} in report'
            )
            report_project_reps = {row['Replicate'] for row in report_grouped[None if single_batch else project]}
            self.assertEqual(
                report_project_reps, set(self.project_replicates[project]),
                f'Mismatched replicates for project {project} in report'
            )


class TestConfig(TestValidateSetup):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.config_dir = f'{cls.data_dir}/validate_pipeline_params/config_files'
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_validate_pipeline_config'
        cls.prog = 'dia_qc validate'
        cls.common_test_args = [
            'config', '--report-format=json',
            '--schema', f'{cls.config_dir}/nextflow_schema.json'
        ]
        cls.common_test_args += ['--api-key', PANORAMA_PUBLIC_KEY] if MOCK_PANORAMA else []
        setup_functions.make_work_dir(cls.work_dir, clear_dir=True)


    def test_all_local_flat(self):
        project = 'Strap'
        meta_vars = {
            'qc_report.color_vars': ['experiment', 'cellLine'],
            'qc_report.control_key': 'cellLine',
            'qc_report.control_values': ['A549', 'CCRF-CEM', 'COLO-205', 'H226', 'H23', 'HeLa',
                                         'NCI7 4-pool', 'NCI7 7-pool', 'RPMI-8228', 'T47D']
        }
        test_config = f'{self.work_dir}/test_all_local_flat.config'
        test_prefix = 'test_all_local_flat'
        args = self.common_test_args + [test_config, '--report-prefix', f'{test_prefix}_']

        setup_test_config(
            f'{self.config_dir}/template.config', output_path=test_config,
            add_params={
                'quant_spectra_dir': f'{self.local_ms_files}/{project}',
                'quant_spectra_glob': '*400to1000*.raw',
                'chromatogram_library_spectra_dir': f'{self.local_ms_files}/{project}',
                'chromatogram_library_spectra_glob': '*-Lib.raw',
                'fasta': f'{self.local_db}/uniprot_human_jan2021_yeastENO1.fasta',
                'spectral_library': f'{self.local_db}/uniprot_human_jan2021_yeastENO1.dlib',
                'replicate_metadata': f'{self.data_dir}/metadata/Strap_metadata.tsv'
            } | meta_vars
        )
        result = setup_functions.run_main(
            vpp._main, args, self.work_dir, prog=self.prog
        )
        self.assertEqual(result.returncode, 0, result.stderr)

        with self.subTest('Check metadata report'):
            self.check_metadata_report(
                f'{self.work_dir}/{test_prefix}_metadata_validation_report.json',
                **{METADATE_CONFIG_TO_ARG_PARAMS[k]: v for k, v in meta_vars.items()}
            )

        with self.subTest('Check replicate report'):
            self.check_replicate_report(
                f'{self.work_dir}/{test_prefix}_replicate_validation_report.json', project
            )


    def test_all_local_batched(self):
        projects = ['Strap', 'Sp3']
        meta_vars = {
            'qc_report.color_vars': ['experiment', 'cellLine']
        }
        test_prefix = 'test_all_local_batched'
        test_config = f'{self.work_dir}/test_all_local_multi_batch.config'
        args = self.common_test_args + [test_config, '--report-prefix', f'{test_prefix}_']

        setup_test_config(
            f'{self.config_dir}/template.config', output_path=test_config,
            add_params={
                'quant_spectra_dir': {project: f'{self.local_ms_files}/{project}' for project in projects},
                'quant_spectra_glob': '*400to1000*.raw',
                'chromatogram_library_spectra_dir': [f'{self.local_ms_files}/{project}' for project in projects],
                'chromatogram_library_spectra_regex': r'-Lib\.raw$',
                'chromatogram_library_spectra_glob': None,
                'fasta': f'{self.local_db}/uniprot_human_jan2021_yeastENO1.fasta',
                'spectral_library': f'{self.local_db}/uniprot_human_jan2021_yeastENO1.dlib',
                'replicate_metadata': f'{self.data_dir}/metadata/Sp3_Strap_combined_metadata.tsv',
            } | meta_vars
        )
        result = setup_functions.run_main(
            vpp._main, args, self.work_dir, prog=self.prog
        )
        self.assertEqual(result.returncode, 0, result.stderr)

        with self.subTest('Check metadata report'):
            self.check_metadata_report(
                f'{self.work_dir}/{test_prefix}_metadata_validation_report.json',
                **{METADATE_CONFIG_TO_ARG_PARAMS[k]: v for k, v in meta_vars.items()}
            )

        with self.subTest('Check replicate report'):
            self.check_replicate_report(
                f'{self.work_dir}/{test_prefix}_replicate_validation_report.json', projects
            )


    def test_all_remote_batched(self):
        projects = ['Strap', 'Sp3']
        meta_params = {
            'qc_report.color_vars': ['experiment', 'cellLine', 'NCI7std']
        }
        test_prefix = 'test_all_remote_batched'
        test_config = f'{self.work_dir}/test_all_local_multi_batch.config'
        args = self.common_test_args + [test_config, '--report-prefix', f'{test_prefix}_']

        setup_test_config(
            f'{self.config_dir}/template.config', output_path=test_config,
            add_params={
                'quant_spectra_dir': {'Strap': STRAP_URL, 'Sp3': SP3_URL},
                'quant_spectra_glob': None,
                'quant_spectra_regex': '400to1000.+?.mzML$',
                'chromatogram_library_spectra_dir': [STRAP_URL, SP3_URL],
                'chromatogram_library_spectra_regex': '-Lib.raw$',
                'chromatogram_library_spectra_glob': None,
                'fasta': f'{self.local_db}/uniprot_human_jan2021_yeastENO1.fasta',
                'spectral_library': f'{self.local_db}/uniprot_human_jan2021_yeastENO1.dlib',
                'replicate_metadata': f'{self.data_dir}/metadata/Sp3_Strap_combined_metadata.tsv',
            } | meta_params
        )

        if MOCK_PANORAMA:
            with mock.patch('DIA_QC_report.validate_pipeline_params.list_panorama_files',
                            side_effect=mock_list_list_panorama_files) as mock_list_files:
                result = setup_functions.run_main(
                    vpp._main, args, self.work_dir, prog=self.prog
                )
            self.assertEqual(result.returncode, 0, result.stderr)
            mock_list_files.assert_has_calls([mock.call(STRAP_URL, api_key=PANORAMA_PUBLIC_KEY),
                                              mock.call(SP3_URL, api_key=PANORAMA_PUBLIC_KEY)])

        else:
            if not have_internet():
                self.skipTest('No internet connection available for remote files.')

            result = setup_functions.run_main(
                vpp._main, args, self.work_dir, prog=self.prog
            )
            self.assertEqual(result.returncode, 0, result.stderr)


        with self.subTest('Check metadata report'):
            self.check_metadata_report(
                f'{self.work_dir}/{test_prefix}_metadata_validation_report.json',
                **{METADATE_CONFIG_TO_ARG_PARAMS[k]: v for k, v in meta_params.items()}
            )

        with self.subTest('Check replicate report'):
            self.check_replicate_report(
                f'{self.work_dir}/{test_prefix}_replicate_validation_report.json', projects
            )


    def test_all_remote_flat(self):
        project = 'Strap'
        meta_params = {
            'qc_report.control_key': 'cellLine',
            'qc_report.control_values': ['HeLa', 'H23']
        }
        test_prefix = 'test_all_remote_flat'
        test_config = f'{self.work_dir}/{test_prefix}.config'
        args = self.common_test_args + [test_config, '--report-prefix', f'{test_prefix}_']

        setup_test_config(
            f'{self.config_dir}/template.config', output_path=test_config,
            add_params={
                'quant_spectra_dir': STRAP_URL,
                'quant_spectra_glob': None,
                'quant_spectra_regex': '400to1000.+?.mzML$',
                'chromatogram_library_spectra_dir': STRAP_URL,
                'chromatogram_library_spectra_glob': '*-Lib.raw',
                'fasta': f'{self.local_db}/uniprot_human_jan2021_yeastENO1.fasta',
                'spectral_library': f'{self.local_db}/uniprot_human_jan2021_yeastENO1.dlib',
                'replicate_metadata': f'{self.data_dir}/metadata/Strap_metadata.tsv',
            } | meta_params
        )

        if MOCK_PANORAMA:
            with mock.patch('DIA_QC_report.validate_pipeline_params.list_panorama_files',
                            side_effect=mock_list_list_panorama_files) as mock_list_files:
                result = setup_functions.run_main(
                    vpp._main, args, self.work_dir, prog=self.prog
                )
            self.assertEqual(result.returncode, 0, result.stderr)
            mock_list_files.assert_has_calls([mock.call(STRAP_URL, api_key=PANORAMA_PUBLIC_KEY)])

        else:
            if not have_internet():
                self.skipTest('No internet connection available for remote files.')

            result = setup_functions.run_main(
                vpp._main, args, self.work_dir, prog=self.prog
            )
            self.assertEqual(result.returncode, 0, result.stderr)

        with self.subTest('Check metadata report'):
            self.check_metadata_report(
                f'{self.work_dir}/{test_prefix}_metadata_validation_report.json',
                **{METADATE_CONFIG_TO_ARG_PARAMS[k]: v for k, v in meta_params.items()}
            )

        with self.subTest('Check replicate report'):
            self.check_replicate_report(
                f'{self.work_dir}/{test_prefix}_replicate_validation_report.json', project
            )


    def test_no_metadata_batched(self):
        projects = ['Strap', 'Sp3']
        meta_params = {}
        test_prefix = 'test_all_remote_batched'
        test_config = f'{self.work_dir}/test_all_local_multi_batch.config'
        args = self.common_test_args + [test_config, '--report-prefix', f'{test_prefix}_']

        setup_test_config(
            f'{self.config_dir}/template.config', output_path=test_config,
            add_params={
                'quant_spectra_dir': {'Strap': STRAP_URL, 'Sp3': SP3_URL},
                'quant_spectra_glob': None,
                'quant_spectra_regex': '400to1000.+?.mzML$',
                'chromatogram_library_spectra_dir': [STRAP_URL, SP3_URL],
                'chromatogram_library_spectra_regex': '-Lib.raw$',
                'chromatogram_library_spectra_glob': None,
                'fasta': f'{self.local_db}/uniprot_human_jan2021_yeastENO1.fasta',
                'spectral_library': f'{self.local_db}/uniprot_human_jan2021_yeastENO1.dlib',
                'replicate_metadata': f'{self.data_dir}/metadata/Sp3_Strap_combined_metadata.tsv',
            } | meta_params
        )

        if MOCK_PANORAMA:
            with mock.patch('DIA_QC_report.validate_pipeline_params.list_panorama_files',
                            side_effect=mock_list_list_panorama_files) as mock_list_files:
                result = setup_functions.run_main(
                    vpp._main, args, self.work_dir, prog=self.prog
                )
            self.assertEqual(result.returncode, 0, result.stderr)
            mock_list_files.assert_has_calls([mock.call(STRAP_URL, api_key=PANORAMA_PUBLIC_KEY),
                                              mock.call(SP3_URL, api_key=PANORAMA_PUBLIC_KEY)])

        else:
            if not have_internet():
                self.skipTest('No internet connection available for remote files.')

            result = setup_functions.run_main(
                vpp._main, args, self.work_dir, prog=self.prog
            )
            self.assertEqual(result.returncode, 0, result.stderr)


        with self.subTest('Check metadata report'):
            self.check_metadata_report(
                f'{self.work_dir}/{test_prefix}_metadata_validation_report.json',
                **{METADATE_CONFIG_TO_ARG_PARAMS[k]: v for k, v in meta_params.items()}
            )

        with self.subTest('Check replicate report'):
            self.check_replicate_report(
                f'{self.work_dir}/{test_prefix}_replicate_validation_report.json', projects
            )


    def test_no_metadata_flat(self):
        project = 'Strap'
        meta_vars = {}
        test_config = f'{self.work_dir}/test_all_local_flat.config'
        test_prefix = 'test_all_local_flat'
        args = self.common_test_args + [test_config, '--report-prefix', f'{test_prefix}_']

        setup_test_config(
            f'{self.config_dir}/template.config', output_path=test_config,
            add_params={
                'quant_spectra_dir': f'{self.local_ms_files}/{project}',
                'quant_spectra_glob': '*400to1000*.raw',
                'chromatogram_library_spectra_dir': f'{self.local_ms_files}/{project}',
                'chromatogram_library_spectra_glob': '*-Lib.raw',
                'fasta': f'{self.local_db}/uniprot_human_jan2021_yeastENO1.fasta',
                'spectral_library': f'{self.local_db}/uniprot_human_jan2021_yeastENO1.dlib',
                'replicate_metadata': f'{self.data_dir}/metadata/Strap_metadata.tsv'
            } | meta_vars
        )
        result = setup_functions.run_main(
            vpp._main, args, self.work_dir, prog=self.prog
        )
        self.assertEqual(result.returncode, 0, result.stderr)

        with self.subTest('Check metadata report'):
            self.check_metadata_report(
                f'{self.work_dir}/{test_prefix}_metadata_validation_report.json',
                **{METADATE_CONFIG_TO_ARG_PARAMS[k]: v for k, v in meta_vars.items()}
            )

        with self.subTest('Check replicate report'):
            self.check_replicate_report(
                f'{self.work_dir}/{test_prefix}_replicate_validation_report.json', project
            )


class TestConfigInvalid(unittest.TestCase):
    def setUp(self):
        pass


    def test_missing_meta_var(self):
        pass


    def test_invalid_config_var_type(self):
        pass


class TestParams(unittest.TestCase):

    def test_batched(self):
        pass


    def test_flat(self):
        pass