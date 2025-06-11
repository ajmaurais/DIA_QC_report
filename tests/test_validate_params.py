
import unittest
from unittest import mock
from types import SimpleNamespace
import os
import re
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

METADATA_CONFIG_TO_ARG_PARAMS = {
    'qc_report.color_vars': 'color_vars', 'batch_report.covariate_vars': 'covariate_vars',
    'batch_report.batch1': 'batch1', 'batch_report.batch2': 'batch2',
    'batch_report.control_key': 'control_key', 'batch_report.control_values': 'control_values'
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


def _is_null(value):
    return value is None or pd.isna(value) or value == ''


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
        self, report_path, metadata_types=None, metadata_df=None, total_reps=None,
        color_vars=None, control_key=None, control_values=None,
        batch1=None, batch2=None, covariate_vars=None
    ):
        if metadata_types is None:
            metadata_types = self.metadata_types
        if metadata_df is None:
            metadata_df = self.combined_metadata

        metadata_counts = {}
        for row in metadata_df.itertuples():
            if row.annotationKey not in metadata_counts:
                metadata_counts[row.annotationKey] = 0
            if not(_is_null(row.annotationValue)):
                metadata_counts[row.annotationKey] += 1

        meta_params = []
        if color_vars:
            if isinstance(color_vars, str):
                color_vars = [color_vars]
            meta_params = [(var, 'color_vars') for var in color_vars]
        if covariate_vars:
            if isinstance(covariate_vars, str):
                covariate_vars = [covariate_vars]
            meta_params += [(var, 'covariate_vars') for var in covariate_vars]
        if control_key:
            meta_params.append((control_key, 'control_key'))
        if batch1:
            meta_params.append((batch1, 'batch1'))
        if batch2:
            meta_params.append((batch2, 'batch2'))

        self.assertTrue(os.path.exists(report_path), f'Metadata report does not exist: {report_path}')
        with open(report_path, 'r') as f:
            report = json.load(f)
        total_reps = len(metadata_df['Replicate'].drop_duplicates())

        self.assertEqual(len(report), len(meta_params))
        for row in report:
            for key in ['variable', 'parameter', 'type', 'missing_in_n_replicates', 'found_in_n_replicates']:
                self.assertIn(key, row, f'Missing key in metadata report row: {key}')

            param_key = (row['variable'], row['parameter'])
            self.assertIn(param_key, meta_params, f'Unexpected metadata parameter: {param_key}')
            self.assertEqual(row['type'], str(metadata_types[row['variable']]),
                             f'Mismatched type for variable {row["variable"]} in metadata report')

            self.assertEqual(
                row['missing_in_n_replicates'], total_reps - metadata_counts[row['variable']],
                f'Mismatched missing count for variable {row["variable"]} in metadata report'
            )


    def check_replicate_report(self, report_path, projects, metadata_df=None, multi_batch=False):
        if isinstance(projects, str):
            projects = [projects]

        if metadata_df is None:
            metadata_df = self.combined_metadata
        metadata_dict = {}
        for row in metadata_df.itertuples():
            if row.Replicate not in metadata_dict:
                metadata_dict[row.Replicate] = {}
            metadata_dict[row.Replicate][row.annotationKey] = row.annotationValue

        self.assertTrue(os.path.exists(report_path), f'Replicate report does not exist: {report_path}')
        with open(report_path, 'r') as f:
            report = json.load(f)

        report_grouped = {}
        for row in report:
            batch = row.pop('ParameterBatch')
            if batch not in report_grouped:
                report_grouped[batch] = []
            report_grouped[batch].append(row)

        target_rept_grouped = {}
        for project in projects:
            if multi_batch:
                target_rept_grouped[project] = self.project_replicates[project]
            else:
                if None not in target_rept_grouped:
                    target_rept_grouped[None] = []
                target_rept_grouped[None].extend(self.project_replicates[project])

        self.assertEqual(len(report_grouped), len(target_rept_grouped), 'Unexpected number of batches in report')
        for project in target_rept_grouped:
            if multi_batch:
                self.assertIn(project, report_grouped, f'Missing batch: {project}')
            else:
                self.assertIn(None, report_grouped, f'Report should have only 1 batch')

            self.assertEqual(
                len(report_grouped[project]), len(target_rept_grouped[project]),
                f'Unexpected number of replicates for project {project} in report'
            )
            report_project_reps = {row['Replicate'] for row in report_grouped[project]}
            self.assertEqual(
                report_project_reps, set(target_rept_grouped[project]),
                f'Mismatched replicates for project {project} in report'
            )

            # Check replicate metadata
            for row in report_grouped[project]:
                replicate = row['Replicate']
                self.assertIn(replicate, metadata_dict, f'Missing metadata for replicate {replicate}')
                meta_keys = set(row.keys()) - {'ParameterBatch', 'Replicate', 'FileName'}
                for key in meta_keys:
                    self.assertIn(key, row, f'Missing key {key} in replicate report row')
                    if _is_null(row[key]):
                        self.assertTrue(_is_null(metadata_dict[replicate][key]))
                    else:
                        self.assertEqual(
                            row[key], metadata_dict[replicate][key],
                            f'Mismatched value for {key} in replicate report row for {replicate}'
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
        cls.local_db_files = {
            'fasta': f'{cls.local_db}/uniprot_human_jan2021_yeastENO1.fasta',
            'spectral_library': f'{cls.local_db}/uniprot_human_jan2021_yeastENO1.dlib',
            'replicate_metadata': f'{cls.data_dir}/metadata/Sp3_Strap_combined_metadata.tsv'
        }
        setup_functions.make_work_dir(cls.work_dir, clear_dir=True)


    def test_all_local_flat(self):
        project = 'Strap'
        meta_params = {
            'qc_report.color_vars': ['experiment', 'cellLine'],
            'batch_report.control_key': 'cellLine',
            'batch_report.control_values': ['A549', 'CCRF-CEM', 'COLO-205', 'H226', 'H23', 'HeLa',
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
                'chromatogram_library_spectra_glob': '*-Lib.raw'
            } | meta_params | self.local_db_files
        )
        result = setup_functions.run_main(
            vpp._main, args, self.work_dir, prog=self.prog
        )
        self.assertEqual(result.returncode, 0, result.stderr)

        with self.subTest('Check metadata report'):
            self.check_metadata_report(
                f'{self.work_dir}/{test_prefix}_metadata_validation_report.json',
                **{METADATA_CONFIG_TO_ARG_PARAMS[k]: v for k, v in meta_params.items()}
            )

        with self.subTest('Check replicate report'):
            self.check_replicate_report(
                f'{self.work_dir}/{test_prefix}_replicate_validation_report.json',
                project, multi_batch=False
            )


    def test_all_local_batched(self):
        projects = ['Strap', 'Sp3']
        meta_params = {
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
                'chromatogram_library_spectra_glob': None
            } | meta_params | self.local_db_files
        )
        result = setup_functions.run_main(
            vpp._main, args, self.work_dir, prog=self.prog
        )
        self.assertEqual(result.returncode, 0, result.stderr)

        with self.subTest('Check metadata report'):
            self.check_metadata_report(
                f'{self.work_dir}/{test_prefix}_metadata_validation_report.json',
                **{METADATA_CONFIG_TO_ARG_PARAMS[k]: v for k, v in meta_params.items()}
            )

        with self.subTest('Check replicate report'):
            self.check_replicate_report(
                f'{self.work_dir}/{test_prefix}_replicate_validation_report.json',
                projects, multi_batch=True
            )


    def test_all_remote_batched(self):
        projects = ['Strap', 'Sp3']
        meta_params = {
            'qc_report.color_vars': 'cellLine',
            'batch_report.covariate_vars': ['experiment', 'cellLine']
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
                'chromatogram_library_spectra_glob': None
            } | meta_params | self.local_db_files
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
                **{METADATA_CONFIG_TO_ARG_PARAMS[k]: v for k, v in meta_params.items()}
            )

        with self.subTest('Check replicate report'):
            self.check_replicate_report(
                f'{self.work_dir}/{test_prefix}_replicate_validation_report.json',
                projects, multi_batch=True
            )


    def test_mixed_batched(self):
        projects = ['Strap', 'Sp3']
        meta_params = {
            'qc_report.color_vars': ['experiment', 'cellLine', 'NCI7std'],
            'batch_report.covariate_vars': 'cellLine'
        }
        test_prefix = 'test_all_remote_batched'
        test_config = f'{self.work_dir}/test_all_local_multi_batch.config'
        args = self.common_test_args + [test_config, '--report-prefix', f'{test_prefix}_']

        setup_test_config(
            f'{self.config_dir}/template.config', output_path=test_config,
            add_params={
                'quant_spectra_dir': {'Strap': f'{self.local_ms_files}/Strap', 'Sp3': SP3_URL},
                'quant_spectra_glob': None,
                'quant_spectra_regex': '400to1000.+?.raw$',
                'chromatogram_library_spectra_dir': [f'{self.local_ms_files}/Strap', SP3_URL],
                'chromatogram_library_spectra_regex': '-Lib.mzML$',
                'chromatogram_library_spectra_glob': None
            } | meta_params | self.local_db_files
        )

        if MOCK_PANORAMA:
            with mock.patch('DIA_QC_report.validate_pipeline_params.list_panorama_files',
                            side_effect=mock_list_list_panorama_files) as mock_list_files:
                result = setup_functions.run_main(
                    vpp._main, args, self.work_dir, prog=self.prog
                )
            self.assertEqual(result.returncode, 0, result.stderr)
            mock_list_files.assert_has_calls([mock.call(SP3_URL, api_key=PANORAMA_PUBLIC_KEY)])

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
                **{METADATA_CONFIG_TO_ARG_PARAMS[k]: v for k, v in meta_params.items()}
            )

        with self.subTest('Check replicate report'):
            self.check_replicate_report(
                f'{self.work_dir}/{test_prefix}_replicate_validation_report.json',
                projects, multi_batch=True
            )


    def test_all_remote_flat(self):
        project = ['Strap', 'Sp3']
        meta_params = {
            'batch_report.control_key': 'cellLine',
            'batch_report.control_values': ['HeLa', 'H23', 'A549', 'H226'],
            'batch_report.batch1': 'experiment'
        }
        test_prefix = 'test_all_remote_flat'
        test_config = f'{self.work_dir}/{test_prefix}.config'
        args = self.common_test_args + [test_config, '--report-prefix', f'{test_prefix}_']

        setup_test_config(
            f'{self.config_dir}/template.config', output_path=test_config,
            add_params={
                'quant_spectra_dir': [STRAP_URL, SP3_URL],
                'quant_spectra_glob': None,
                'quant_spectra_regex': '400to1000.+?.mzML$'
            } | meta_params | self.local_db_files
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
                **{METADATA_CONFIG_TO_ARG_PARAMS[k]: v for k, v in meta_params.items()}
            )

        with self.subTest('Check replicate report'):
            self.check_replicate_report(
                f'{self.work_dir}/{test_prefix}_replicate_validation_report.json',
                project, multi_batch=False
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
                'chromatogram_library_spectra_glob': None
            } | meta_params | self.local_db_files
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
                **{METADATA_CONFIG_TO_ARG_PARAMS[k]: v for k, v in meta_params.items()}
            )

        with self.subTest('Check replicate report'):
            self.check_replicate_report(
                f'{self.work_dir}/{test_prefix}_replicate_validation_report.json',
                projects, multi_batch=True
            )


    def test_no_metadata_flat(self):
        project = 'Strap'
        meta_params = {}
        test_config = f'{self.work_dir}/test_all_local_flat.config'
        test_prefix = 'test_all_local_flat'
        args = self.common_test_args + ['--report-prefix', f'{test_prefix}_', test_config]

        setup_test_config(
            f'{self.config_dir}/template.config', output_path=test_config,
            add_params={
                'quant_spectra_dir': f'{self.local_ms_files}/{project}',
                'quant_spectra_glob': '*400to1000*.raw',
                'chromatogram_library_spectra_dir': f'{self.local_ms_files}/{project}',
                'chromatogram_library_spectra_glob': '*-Lib.raw'
            } | meta_params | self.local_db_files
        )
        result = setup_functions.run_main(
            vpp._main, args, self.work_dir, prog=self.prog
        )
        self.assertEqual(result.returncode, 0, result.stderr)

        with self.subTest('Check metadata report'):
            self.check_metadata_report(
                f'{self.work_dir}/{test_prefix}_metadata_validation_report.json',
                **{METADATA_CONFIG_TO_ARG_PARAMS[k]: v for k, v in meta_params.items()}
            )

        with self.subTest('Check replicate report'):
            self.check_replicate_report(
                f'{self.work_dir}/{test_prefix}_replicate_validation_report.json', project
            )


    def test_remote_schema(self):
        project = 'Strap'
        meta_params = {}
        test_config = f'{self.work_dir}/test_all_local_flat.config'
        test_prefix = 'test_all_local_flat'
        args = [
            'config', test_config, '--report-prefix', f'{test_prefix}_',
            '--pipeline', 'ajmaurais/nf-skyline-dia-ms', '--revision', 'nf-schema'
        ]

        setup_test_config(
            f'{self.config_dir}/template.config', output_path=test_config,
            add_params={
                'quant_spectra_dir': f'{self.local_ms_files}/{project}',
                'quant_spectra_glob': '*400to1000*.raw',
                'chromatogram_library_spectra_dir': f'{self.local_ms_files}/{project}',
                'chromatogram_library_spectra_glob': '*-Lib.raw'
            } | meta_params | self.local_db_files
        )
        result = setup_functions.run_main(
            vpp._main, args, self.work_dir, prog=self.prog
        )
        self.assertEqual(result.returncode, 0, result.stderr)


    def test_missing_meta_values(self):
        project = 'Strap'
        meta_params = {
            'qc_report.color_vars': ['string_var', 'bool_var', 'int_var', 'float_var', 'na_var']
        }
        test_config = f'{self.work_dir}/test_missing_meta_values.config'
        test_prefix = 'test_missing_meta_values'
        args = self.common_test_args + [test_config, '--report-prefix', f'{test_prefix}_']
        replicate_metadata = f'{self.data_dir}/metadata/Strap_missing_multi_var_metadata.tsv'

        setup_test_config(
            f'{self.config_dir}/template.config', output_path=test_config,
            add_params={
                'quant_spectra_dir': f'{self.local_ms_files}/{project}',
                'quant_spectra_glob': '*400to1000*.raw',
                'fasta': f'{self.local_db}/uniprot_human_jan2021_yeastENO1.fasta',
                'spectral_library': f'{self.local_db}/uniprot_human_jan2021_yeastENO1.dlib',
                'replicate_metadata': replicate_metadata
            } | meta_params
        )
        result = setup_functions.run_main(
            vpp._main, args, self.work_dir, prog=self.prog
        )
        self.assertEqual(result.returncode, 0, result.stderr)

        meta_reader = Metadata()
        if not meta_reader.read(replicate_metadata):
            raise FileNotFoundError(f"Metadata file '{replicate_metadata}' not found or could not be read.")

        for var in meta_params['qc_report.color_vars']:
            self.assertRegex(
                result.stdout,
                rf"(?m)^\[WARNING\]\s+{str(meta_reader.types[var])} variable '{var}' is missing in [0-9]+ of [0-9]+ replicates\.$",
            )

        with self.subTest('Check metadata report'):
            self.check_metadata_report(
                f'{self.work_dir}/{test_prefix}_metadata_validation_report.json',
                **{METADATA_CONFIG_TO_ARG_PARAMS[k]: v for k, v in meta_params.items()},
                metadata_df=meta_reader.df, metadata_types=meta_reader.types
            )

        with self.subTest('Check replicate report'):
            self.check_replicate_report(
                f'{self.work_dir}/{test_prefix}_replicate_validation_report.json', project,
                metadata_df=meta_reader.df, multi_batch=False
            )


class TestConfigInvalid(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data_dir = f'{setup_functions.TEST_DIR}/data'
        cls.config_dir = f'{cls.data_dir}/validate_pipeline_params/config_files'
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_validate_pipeline_config_invalid'
        setup_functions.make_work_dir(cls.work_dir, clear_dir=True)


    def test_invalid_pipeline_args(self):
        test_config = f'{self.config_dir}/template.config'

        with self.subTest('Both --pipeline and --schema'):
            args = [
                'config', '--schema', f'{self.config_dir}/nextflow_schema.json',
                '--pipeline', vpp.DEFAULT_PIPELINE, test_config
            ]
            result = setup_functions.run_main(vpp._main, args, self.work_dir, prog='dia_qc validate')
            self.assertEqual(result.returncode, 2)
            self.assertIn('--pipeline and --schema options conflict', result.stderr)


        with self.subTest('--revision on local pipeline'):
            args = [
                'config', '--revision', 'main',
                '--pipeline', f'{self.config_dir}/template.config',
                test_config
            ]
            result = setup_functions.run_main(vpp._main, args, self.work_dir, prog='dia_qc validate')
            self.assertEqual(result.returncode, 2)
            self.assertIn('--revision cannot be used for a local pipeline script.', result.stderr)


    def test_missing_meta_var(self):
        pass


    def test_invalid_config_var_type(self):
        pass


class TestParams(TestValidateSetup):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_validate_params'
        cls.common_test_args = [
            '--report-format=json'
        ]
        setup_functions.make_work_dir(cls.work_dir, clear_dir=True)


    def write_ms_file_json(self, projects, output_path,
                           multi_batch=False, file_regex=None):
        if isinstance(projects, str):
            projects = [projects]

        if multi_batch:
            files = {}
            for project in projects:
                files[project] = [
                    f for f in os.listdir(f'{self.local_ms_files}/{project}')
                    if re.search(file_regex, f)
                ]
        else:
            files = []
            for project in projects:
                files.extend([
                    f for f in os.listdir(f'{self.local_ms_files}/{project}')
                    if re.search(file_regex, f)
                ])

        with open(output_path, 'w') as f:
            json.dump(files, f, indent=4)


    def test_flat(self):
        project = 'Strap'
        color_vars = ['experiment', 'cellLine']
        test_prefix = 'test_flat'
        quant_file_args = [
            f'-q={f}' for f in os.listdir(f'{self.local_ms_files}/{project}')
            if re.search(r'400to1000.+?.raw$', f)
        ]
        args = [
            'params', '--report-prefix', f'{test_prefix}_',
            '--metadata', f'{setup_functions.TEST_DIR}/data/metadata/Strap_metadata.json',
        ] + self.common_test_args + [f'--addColorVar={var}' for var in color_vars] + quant_file_args

        restult = setup_functions.run_main(
            vpp._main, args, self.work_dir, prog='dia_qc validate'
        )
        self.assertEqual(restult.returncode, 0, restult.stderr)

        with self.subTest('Check metadata report'):
            self.check_metadata_report(
                f'{self.work_dir}/{test_prefix}_metadata_validation_report.json',
                color_vars=color_vars
            )

        with self.subTest('Check replicate report'):
            self.check_replicate_report(
                f'{self.work_dir}/{test_prefix}_replicate_validation_report.json', project
            )


    def test_flat_json(self):
        project = 'Strap'
        color_vars = ['experiment', 'cellLine']
        test_prefix = 'test_flat'
        quant_files_path = f'{self.work_dir}/{test_prefix}_quant_files.json'
        chrom_files_path = f'{self.work_dir}/{test_prefix}_lib_files.json'
        self.write_ms_file_json(project, quant_files_path, file_regex=r'400to1000.+?.raw$')
        self.write_ms_file_json(project, chrom_files_path, file_regex=r'-Lib\.raw$')
        args = [
            'params', '--report-prefix', f'{test_prefix}_',
            '--quant-spectra-json', quant_files_path, '--chrom-lib-spectra-json', chrom_files_path,
            '--metadata', f'{setup_functions.TEST_DIR}/data/metadata/Strap_metadata.json',
        ] + self.common_test_args + [f'--addColorVar={var}' for var in color_vars]

        restult = setup_functions.run_main(
            vpp._main, args, self.work_dir, prog='dia_qc validate'
        )
        self.assertEqual(restult.returncode, 0, restult.stderr)

        with self.subTest('Check metadata report'):
            self.check_metadata_report(
                f'{self.work_dir}/{test_prefix}_metadata_validation_report.json',
                color_vars=color_vars
            )

        with self.subTest('Check replicate report'):
            self.check_replicate_report(
                f'{self.work_dir}/{test_prefix}_replicate_validation_report.json', project
            )


    def test_batched(self):
        projects = ['Strap', 'Sp3']
        color_vars = ['experiment', 'cellLine']

        test_prefix = 'test_batched'
        quant_files_path = f'{self.work_dir}/{test_prefix}_quant_files.json'
        chrom_files_path = f'{self.work_dir}/{test_prefix}_lib_files.json'
        self.write_ms_file_json(
            projects, quant_files_path, multi_batch=True, file_regex=r'400to1000.+?.raw$'
        )
        self.write_ms_file_json(projects, chrom_files_path, file_regex=r'-Lib\.raw$')
        args = [
            'params', '--report-prefix', f'{test_prefix}_',
            '--quant-spectra-json', quant_files_path, '--chrom-lib-spectra-json', chrom_files_path,
            '--metadata', f'{setup_functions.TEST_DIR}/data/metadata/Sp3_Strap_combined_metadata.tsv'
        ] + self.common_test_args + [f'--addColorVar={var}' for var in color_vars]

        restult = setup_functions.run_main(
            vpp._main, args, self.work_dir, prog='dia_qc validate'
        )
        self.assertEqual(restult.returncode, 0, restult.stderr)

        with self.subTest('Check metadata report'):
            self.check_metadata_report(
                f'{self.work_dir}/{test_prefix}_metadata_validation_report.json',
                color_vars=color_vars
            )

        with self.subTest('Check replicate report'):
            self.check_replicate_report(
                f'{self.work_dir}/{test_prefix}_replicate_validation_report.json',
                projects, multi_batch=True
            )