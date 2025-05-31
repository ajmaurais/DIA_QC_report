
import unittest
from unittest import mock

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