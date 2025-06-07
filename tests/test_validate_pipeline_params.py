
import unittest
from unittest import mock
from types import SimpleNamespace

from DIA_QC_report import validate_pipeline_params as vpp
from DIA_QC_report.submodules.pipeline_config import PipelineConfig
from DIA_QC_report.submodules.panorama import have_internet
from DIA_QC_report.submodules.panorama import PANORA_PUBLIC_KEY

import setup_functions

STRAP_URL = 'https://panoramaweb.org/_webdav/ICPC/NCI-7%20Joint%20Project/NCI-7%20Data%20Harmonization/LFQ-Analyses/USA-UW/%40files/RawFiles/S-Trap/'
SP3_URL = 'https://panoramaweb.org/_webdav/ICPC/NCI-7%20Joint%20Project/NCI-7%20Data%20Harmonization/LFQ-Analyses/USA-UW/%40files/RawFiles/SP3/'
PUBLIC_URL = 'https://panoramaweb.org/_webdav/Panorama%20Public/2024/Thermo%20Fisher%20Research%20and%20Development%20-%202024_Stellar_Instrument_Platform/Label%20Free%20-%20E.%20coli/%40files/RawFiles/ReplicatesSmall/'


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
            vpp._main, ['config', '--debug', '--schema', f'{self.config_dir}/nextflow_schema.json', test_config],
            self.work_dir, prog=self.prog
        )
        self.assertEqual(result.returncode, 0, msg=result.stderr)


    def test_all_local_batched(self):
        pass


    def test_all_remote_batched(self):
        pass


    def test_all_remote_flat(self):
        pass


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