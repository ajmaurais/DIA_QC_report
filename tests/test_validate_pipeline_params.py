
import unittest
from unittest import mock

from DIA_QC_report import validate_pipeline_params as vpp
from DIA_QC_report.submodules.panorama import have_internet
from DIA_QC_report.submodules.panorama import PANORA_PUBLIC_KEY

import setup_functions

STRAP_URL = 'https://panoramaweb.org/_webdav/ICPC/NCI-7%20Joint%20Project/NCI-7%20Data%20Harmonization/LFQ-Analyses/USA-UW/%40files/RawFiles/S-Trap/'
SP3_URL = 'https://panoramaweb.org/_webdav/ICPC/NCI-7%20Joint%20Project/NCI-7%20Data%20Harmonization/LFQ-Analyses/USA-UW/%40files/RawFiles/SP3/'
PUBLIC_URL = 'https://panoramaweb.org/_webdav/Panorama%20Public/2024/Thermo%20Fisher%20Research%20and%20Development%20-%202024_Stellar_Instrument_Platform/Label%20Free%20-%20E.%20coli/%40files/RawFiles/ReplicatesSmall/'


class TestMainConfig(unittest.TestCase):
    def setup_test_config(base_config, output_path, add_params=None):
        pass

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