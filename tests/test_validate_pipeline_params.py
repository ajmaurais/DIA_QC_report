
import unittest
from unittest import mock
from types import SimpleNamespace

from DIA_QC_report import validate_pipeline_params as vpp
from DIA_QC_report.submodules import pipeline_config as npc
from DIA_QC_report.submodules.panorama import have_internet
from DIA_QC_report.submodules.panorama import PANORA_PUBLIC_KEY

import setup_functions

STRAP_URL = 'https://panoramaweb.org/_webdav/ICPC/NCI-7%20Joint%20Project/NCI-7%20Data%20Harmonization/LFQ-Analyses/USA-UW/%40files/RawFiles/S-Trap/'
SP3_URL = 'https://panoramaweb.org/_webdav/ICPC/NCI-7%20Joint%20Project/NCI-7%20Data%20Harmonization/LFQ-Analyses/USA-UW/%40files/RawFiles/SP3/'
PUBLIC_URL = 'https://panoramaweb.org/_webdav/Panorama%20Public/2024/Thermo%20Fisher%20Research%20and%20Development%20-%202024_Stellar_Instrument_Platform/Label%20Free%20-%20E.%20coli/%40files/RawFiles/ReplicatesSmall/'


def setup_test_config(base_config, output_path, add_params=None):
    config = npc.parse_config(base_config)

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
        add_attribute(config, path, value)

    npc.write_params_namespace(config, output_path)


class TestMainConfig(unittest.TestCase):

    def test_main(self):
        pass


    def test_all_local_batched(self):
        pass


    def test_all_remote_batched(self):
        pass


    def test_all_local_flat(self):
        pass


    def test_all_remote_flat(self):
        pass


    def test_no_metadata_batched(self):
        pass


    def test_no_metadata_flat(self):
        pass


class TestMainConfigInvalid(unittest.TestCase):
    def setUp(self):
        pass


    def test_missing_meta_var(self):
        pass


class TestMainParams(unittest.TestCase):

    def test_batched(self):
        pass


    def test_flat(self):
        pass