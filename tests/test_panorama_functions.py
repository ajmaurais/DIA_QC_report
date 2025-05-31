
import unittest
from unittest import mock
import os

from DIA_QC_report.submodules import panorama
import setup_functions

class TestListPanoramaFiles(unittest.TestCase):
    def setUp(self):
        self.data_dir = f'{setup_functions.TEST_DIR}/data/'
        self.target_files = {
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

        with open(f'{self.data_dir}/validate_pipeline_params/panorama/flat_response.xml', 'r') as f:
            self.flat_response = f.read()


    @staticmethod
    def _make_mock_response(text):
        mock = unittest.mock.MagicMock()
        mock.ok = True
        mock.status_code = 207 # WebDAV multistatus
        mock.text = text
        return mock


    def test_list_panorama_files_anonymous(self):
        with mock.patch('DIA_QC_report.submodules.panorama.requests.request',
                                 return_value=self._make_mock_response(self.flat_response)):
            files = panorama.list_panorama_files(
                url='https://panoramaweb.org/_webdav/dummy/path'
            )
            self.assertIsInstance(files, list)
            self.assertEqual(set(files), self.target_files)


    def test_list_panorama_files_api_key(self):
        with mock.patch('DIA_QC_report.submodules.panorama.requests.Session.request',
                                 return_value=self._make_mock_response(self.flat_response)):
            files = panorama.list_panorama_files(
                url='https://panoramaweb.org/_webdav/dummy/path',
                api_key=panorama.PANORA_PUBLIC_KEY
            )
            self.assertIsInstance(files, list)
            self.assertEqual(set(files), self.target_files)


    def test_list_panorama_files_full_url(self):
        test_data_url = 'https://panoramaweb.org/_webdav/Panorama Public/2024/Thermo Fisher Research and Development - 2024_Stellar_Instrument_Platform/Label Free - E. coli/@files/RawFiles/ReplicatesSmall/'
        target_files = {f'{test_data_url}{file_name}' for file_name in self.target_files}

        with mock.patch('DIA_QC_report.submodules.panorama.requests.request',
                                 return_value=self._make_mock_response(self.flat_response)):
            files = panorama.list_panorama_files(
                url='https://panoramaweb.org/_webdav/dummy/path',
                full_url=True
            )
            self.assertIsInstance(files, list)
            self.assertEqual(set(files), target_files)


class TestDownloadWebDAVFile(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_download_panorama_file/'

        setup_functions.make_work_dir(cls.work_dir, clear_dir=True)

    def _fake_response(self, status_code=200, data=b''):
        resp = unittest.mock.MagicMock()
        resp.status_code = status_code
        resp.ok = status_code == 200
        resp.iter_content = lambda chunk_size: [data[i:i + chunk_size]
                                                for i in range(0, len(data), chunk_size)]
        # When status_code != 200 we need .text for error printing
        resp.text = data.decode() if isinstance(data, bytes) else data
        return resp


    def test_download_with_api_key(self):
        test_file = 'Sp3_metadata.json'
        with open(f'{self.data_dir}/metadata/{test_file}', 'rb') as inF:
            payload = inF.read()

        with mock.patch('DIA_QC_report.submodules.panorama.requests.Session.get',
                   return_value=self._fake_response(data=payload)) as mock_get:
            out_path = panorama.download_webdav_file(
                f'https://server/_webdav/lab/%40files/{test_file}',
                dest_path=f'{self.work_dir}/{test_file}',
                api_key='DUMMY_KEY'
            )

            # ensure file is written
            self.assertTrue(os.path.isfile(out_path))
            with open(out_path, 'rb') as outF:
                self.assertEqual(outF.read(), payload)

            # verify the URL passed to get()
            mock_get.assert_called_once_with(
                f'https://server/_webdav/lab/%40files/{test_file}',
                headers={'Authorization': 'Basic YXBpa2V5OkRVTU1ZX0tFWQ=='},
                stream=True,
                verify=True
            )


    def test_download_anonymous(self):
        test_file = 'Sp3_metadata.tsv'
        with open(f'{self.data_dir}/metadata/{test_file}', 'rb') as inF:
            payload = inF.read()

        with mock.patch('DIA_QC_report.submodules.panorama.requests.get',
                   return_value=self._fake_response(data=payload)) as mock_get:
            out_path = panorama.download_webdav_file(
                f'https://server/_webdav/lab/%40files/{test_file}',
                dest_path=f'{self.work_dir}/{test_file}'
            )

            # ensure file is written
            self.assertTrue(os.path.isfile(out_path))
            with open(out_path, 'rb') as outF:
                self.assertEqual(outF.read(), payload)

            # verify the URL passed to get()
            mock_get.assert_called_once_with(
                f'https://server/_webdav/lab/%40files/{test_file}',
                headers={},
                stream=True,
                verify=True
            )


    def test_download_file_to_string(self):
        test_file = 'Sp3_metadata.tsv'
        with open(f'{self.data_dir}/metadata/{test_file}', 'rb') as inF:
            payload = inF.read()

        with mock.patch('DIA_QC_report.submodules.panorama.requests.get',
                   return_value=self._fake_response(data=payload)) as mock_get:
            file_text = panorama.download_webdav_file(
                f'https://server/_webdav/lab/%40files/{test_file}',
                return_text=True
            )

            target_text = payload.decode('utf-8')
            self.assertEqual(file_text, target_text)

            # verify the URL passed to get()
            mock_get.assert_called_once_with(
                f'https://server/_webdav/lab/%40files/{test_file}',
                headers={},
                stream=True,
                verify=True
            )