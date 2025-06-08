
import unittest
from unittest import mock
from types import SimpleNamespace
import os

from DIA_QC_report import validate_pipeline_params as vpp
from DIA_QC_report.submodules.pipeline_config import PipelineConfig
from DIA_QC_report.submodules.panorama import have_internet
from DIA_QC_report.submodules.panorama import PANORAMA_PUBLIC_KEY

import setup_functions

STRAP_URL = 'https://panoramaweb.org/_webdav/ICPC/NCI-7%20Joint%20Project/NCI-7%20Data%20Harmonization/LFQ-Analyses/USA-UW/%40files/RawFiles/S-Trap/'
SP3_URL = 'https://panoramaweb.org/_webdav/ICPC/NCI-7%20Joint%20Project/NCI-7%20Data%20Harmonization/LFQ-Analyses/USA-UW/%40files/RawFiles/SP3/'
PUBLIC_URL = 'https://panoramaweb.org/_webdav/Panorama%20Public/2024/Thermo%20Fisher%20Research%20and%20Development%20-%202024_Stellar_Instrument_Platform/Label%20Free%20-%20E.%20coli/%40files/RawFiles/ReplicatesSmall/'

MOCK_PANORAMA = True

def mock_list_list_panorama_files(url, **kwargs):
    if url == STRAP_URL:
        return os.listdir(setup_functions.TEST_DIR + '/data/validate_pipeline_params/mock_local_ms_files/Strap')
    elif url == SP3_URL:
        return os.listdir(setup_functions.TEST_DIR + '/data/validate_pipeline_params/mock_local_ms_files/SP3')
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


class TestConfig(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data_dir = f'{setup_functions.TEST_DIR}/data'
        cls.config_dir = f'{cls.data_dir}/validate_pipeline_params/config_files'
        cls.local_ms_files = f'{cls.data_dir}/validate_pipeline_params/mock_local_ms_files'
        cls.local_db = f'{cls.data_dir}/validate_pipeline_params/mock_local_db'
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
        test_config = f'{self.work_dir}/test_all_local_flat.config'
        setup_test_config(
            f'{self.config_dir}/flat.config', output_path=test_config,
            add_params={
                'quant_spectra_dir': f'{self.local_ms_files}/{project}',
                'quant_spectra_glob': '*400to1000*.raw',
                'chromatogram_library_spectra_dir': f'{self.local_ms_files}/{project}',
                'chromatogram_library_spectra_glob': '*-Lib.raw',
                'fasta': f'{self.local_db}/uniprot_human_jan2021_yeastENO1.fasta',
                'spectral_library': f'{self.local_db}/uniprot_human_jan2021_yeastENO1.dlib',
                'replicate_metadata': f'{self.data_dir}/metadata/Strap_metadata.tsv'
            }
        )
        result = setup_functions.run_main(
            vpp._main, ['config', '--schema', f'{self.config_dir}/nextflow_schema.json', test_config],
            self.work_dir, prog=self.prog
        )
        self.assertEqual(result.returncode, 0, result.stderr)


    def test_all_local_batched(self):
        projects = ['Strap', 'SP3']
        test_config = f'{self.work_dir}/test_all_local_multi_batch.config'
        setup_test_config(
            f'{self.config_dir}/multi_batch.config', output_path=test_config,
            add_params={
                'quant_spectra_dir': {project: f'{self.local_ms_files}/{project}' for project in projects},
                'quant_spectra_glob': '*400to1000*.raw',
                'chromatogram_library_spectra_dir': [f'{self.local_ms_files}/{project}' for project in projects],
                'chromatogram_library_spectra_regex': r'-Lib\.raw$',
                'chromatogram_library_spectra_glob': None,
                'qc_report.color_vars': ['experiment', 'cellLine'],
                'fasta': f'{self.local_db}/uniprot_human_jan2021_yeastENO1.fasta',
                'spectral_library': f'{self.local_db}/uniprot_human_jan2021_yeastENO1.dlib',
                'replicate_metadata': f'{self.data_dir}/metadata/Sp3_Strap_combined_metadata.tsv'
            }
        )
        result = setup_functions.run_main(
            vpp._main, self.common_test_args + [test_config],
            self.work_dir, prog=self.prog
        )
        self.assertEqual(result.returncode, 0, result.stderr)


    def test_all_remote_batched(self):
        projects = ['Strap', 'SP3']
        test_config = f'{self.work_dir}/test_all_local_multi_batch.config'
        setup_test_config(
            f'{self.config_dir}/multi_batch.config', output_path=test_config,
            add_params={
                'quant_spectra_dir': {'Strap': STRAP_URL, 'SP3': SP3_URL},
                'quant_spectra_glob': None,
                'quant_spectra_regex': '400to1000.+?.mzML$',
                'chromatogram_library_spectra_dir': [STRAP_URL, SP3_URL],
                'chromatogram_library_spectra_regex': '-Lib.raw$',
                'chromatogram_library_spectra_glob': None,
                'qc_report.color_vars': ['experiment', 'cellLine', 'NCI7std'],
                'fasta': f'{self.local_db}/uniprot_human_jan2021_yeastENO1.fasta',
                'spectral_library': f'{self.local_db}/uniprot_human_jan2021_yeastENO1.dlib',
                'replicate_metadata': f'{self.data_dir}/metadata/Sp3_Strap_combined_metadata.tsv'
            }
        )

        if MOCK_PANORAMA:
            with mock.patch('DIA_QC_report.validate_pipeline_params.list_panorama_files',
                            side_effect=mock_list_list_panorama_files) as mock_list_files:
                result = setup_functions.run_main(
                    vpp._main, self.common_test_args + [test_config],
                    self.work_dir, prog=self.prog
                )
            self.assertEqual(result.returncode, 0, result.stderr)
            mock_list_files.assert_has_calls([mock.call(STRAP_URL, api_key=PANORAMA_PUBLIC_KEY),
                                              mock.call(SP3_URL, api_key=PANORAMA_PUBLIC_KEY)])

        else:
            if not have_internet():
                self.skipTest('No internet connection available for remote files.')

            result = setup_functions.run_main(
                vpp._main, self.common_test_args + [test_config],
                self.work_dir, prog=self.prog
            )
            self.assertEqual(result.returncode, 0, result.stderr)


    def test_all_remote_flat(self):
        project = 'Strap'
        test_config = f'{self.work_dir}/test_all_local_flat.config'
        setup_test_config(
            f'{self.config_dir}/flat.config', output_path=test_config,
            add_params={
                'quant_spectra_dir': STRAP_URL,
                'quant_spectra_glob': None,
                'quant_spectra_regex': '400to1000.+?.mzML$',
                'chromatogram_library_spectra_dir': STRAP_URL,
                'chromatogram_library_spectra_glob': '*-Lib.raw',
                'qc_report.control_key': 'cellLine',
                'qc_report.control_values': ['HeLa', 'H23'],
                'fasta': f'{self.local_db}/uniprot_human_jan2021_yeastENO1.fasta',
                'spectral_library': f'{self.local_db}/uniprot_human_jan2021_yeastENO1.dlib',
                'replicate_metadata': f'{self.data_dir}/metadata/Strap_metadata.tsv'
            }
        )

        if MOCK_PANORAMA:
            with mock.patch('DIA_QC_report.validate_pipeline_params.list_panorama_files',
                            side_effect=mock_list_list_panorama_files) as mock_list_files:
                result = setup_functions.run_main(
                    vpp._main, self.common_test_args + [test_config],
                    self.work_dir, prog=self.prog
                )
            self.assertEqual(result.returncode, 0, result.stderr)
            mock_list_files.assert_has_calls([mock.call(STRAP_URL, api_key=PANORAMA_PUBLIC_KEY)])

        else:
            if not have_internet():
                self.skipTest('No internet connection available for remote files.')

            result = setup_functions.run_main(
                vpp._main, self.common_test_args + [test_config],
                self.work_dir, prog=self.prog
            )
            self.assertEqual(result.returncode, 0, result.stderr)



    def test_no_metadata_batched(self):
        pass


    def test_no_metadata_flat(self):
        pass


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