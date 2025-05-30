
import unittest

from DIA_QC_report.submodules import panorama
import setup_functions

class TestListPanoramaFiles(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'
        cls.target_files = {
            '240419_P1_Neo_60SPD_EColiHeLa_Small_Assay_Rep1.raw',
            '240419_P1_Neo_60SPD_EColiHeLa_Small_Assay_Rep2.raw',
            '240419_P1_Neo_60SPD_EColiHeLa_Small_Assay_Rep3.raw',
            '240419_P1_Neo_60SPD_EColiHeLa_Small_Assay_Rep4.raw',
            '240419_P1_Neo_60SPD_EColiHeLa_Small_Assay_Rep5.raw',
            '240419_P1_Neo_60SPD_EColiHeLa_Small_Assay_Rep6.raw',
            '240419_P1_Neo_60SPD_EColiHeLa_Small_Assay_Rep7.raw',
            '240419_P1_Neo_60SPD_EColiHeLa_Small_Assay_Rep8.raw',
            'injection_times_ecoli_replicates.zip'
        }

        with open(f'{cls.data_dir}/panorama/flat_response.xml', 'r') as f:
            cls.flat_response = f.read()


    @staticmethod
    def _make_mock_response(text):
        mock = unittest.mock.MagicMock()
        mock.ok = True
        mock.status_code = 207 # WebDAV multistatus
        mock.text = text
        return mock


    def test_list_panorama_files_anonymous(self):
        with unittest.mock.patch('DIA_QC_report.submodules.panorama.requests.request',
                                 return_value=self._make_mock_response(self.flat_response)):
            files = panorama.list_panorama_files(
                url='https://panoramaweb.org/_webdav/dummy/path'
            )
            self.assertIsInstance(files, list)
            self.assertEqual(set(files), self.target_files)


    def test_list_panorama_files_api_key(self):
        with unittest.mock.patch('DIA_QC_report.submodules.panorama.requests.Session.request',
                                 return_value=self._make_mock_response(self.flat_response)):
            files = panorama.list_panorama_files(
                url='https://panoramaweb.org/_webdav/dummy/path',
                api_key=panorama.PANORA_PUBLIC_KEY
            )
            self.assertIsInstance(files, list)
            self.assertEqual(set(files), self.target_files)


    def test_list_panorama_files_full_url(self):
        with unittest.mock.patch('DIA_QC_report.submodules.panorama.requests.request',
                                 return_value=self._make_mock_response(self.flat_response)):
            files = panorama.list_panorama_files(
                url='https://panoramaweb.org/_webdav/dummy/path'
                full_url=True
            )
            self.assertIsInstance(files, list)
            self.assertEqual(set(files), self.target_files)