
import unittest
from unittest import mock
import os

from DIA_QC_report.submodules import panorama
import setup_functions

PUBLIC_URL = 'https://panoramaweb.org/_webdav/Panorama%20Public/2024/Thermo%20Fisher%20Research%20and%20Development%20-%202024_Stellar_Instrument_Platform/Label%20Free%20-%20E.%20coli/%40files/RawFiles/ReplicatesSmall/'
PUBLIC_FILE = 'https://panoramaweb.org/_webdav/Panorama%20Public/2022/MacCoss%20-%20Human%20AD%20Clean%20Diagnosis%20DIA%20Data/SMTG/%40files/MetaData/Clean-SMTG-B1-1410-Oct2022-3qtrans-Meta.csv'

MOCK_PANORAMA = True

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


    def test_list_panorama_files(self):
        if MOCK_PANORAMA:
            with mock.patch('DIA_QC_report.submodules.panorama.requests.request',
                            return_value=self._make_mock_response(self.flat_response)) as mock_request:
                files = panorama.list_panorama_files(PUBLIC_URL)

            headers={'Depth': '1', 'Content-Type': 'text/xml; charset=utf-8'}
            call_args = ('PROPFIND', PUBLIC_URL)
            call_kwargs = {
                'headers': headers, 'data': panorama._LIST_BODY,
                'verify': True, 'auth': None
            }
            mock_request.assert_called_once_with(*call_args, **call_kwargs)

        else:
            if not panorama.have_internet():
                self.skipTest("Internet connection is required for this test.")

            files = panorama.list_panorama_files(url=PUBLIC_URL)

        self.assertIsInstance(files, list)
        self.assertEqual(set(files), self.target_files)


    def test_list_panorama_files_api_key(self):
        if MOCK_PANORAMA:
            with mock.patch('DIA_QC_report.submodules.panorama.requests.Session.request',
                            return_value=self._make_mock_response(self.flat_response)) as mock_request:
                files = panorama.list_panorama_files(PUBLIC_URL, api_key=panorama.PANORA_PUBLIC_KEY)

            headers={'Depth': '1', 'Content-Type': 'text/xml; charset=utf-8'}
            token = panorama.b64encode(f'apikey:{panorama.PANORA_PUBLIC_KEY}'.encode()).decode()
            headers['Authorization'] = f'Basic {token}'

            call_args = ('PROPFIND', PUBLIC_URL)
            call_kwargs = {
                'headers': headers, 'data': panorama._LIST_BODY,
                'verify': True, 'auth': ()
            }
            mock_request.assert_called_once_with(*call_args, **call_kwargs)

        else:
            if not panorama.have_internet():
                self.skipTest("Internet connection is required for this test.")

            files = panorama.list_panorama_files(
                url=PUBLIC_URL, api_key=panorama.PANORA_PUBLIC_KEY
            )

        self.assertIsInstance(files, list)
        self.assertEqual(set(files), self.target_files)


    def test_list_panorama_files_full_url(self):
        target_files = {f'{PUBLIC_URL}{file_name}' for file_name in self.target_files}

        with mock.patch('DIA_QC_report.submodules.panorama.requests.request',
                                 return_value=self._make_mock_response(self.flat_response)):
            files = panorama.list_panorama_files(
                url=PUBLIC_URL,
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
        test_file = 'Clean-SMTG-B1-1410-Oct2022-3qtrans-Meta.csv'
        with open(f'{self.data_dir}/validate_pipeline_params/panorama/{test_file}', 'rb') as inF:
            payload = inF.read()

        if MOCK_PANORAMA:
            with mock.patch('DIA_QC_report.submodules.panorama.requests.Session.get',
                    return_value=self._fake_response(data=payload)) as mock_get:
                out_path = panorama.download_webdav_file(
                    PUBLIC_FILE,
                    dest_path=f'{self.work_dir}/{test_file}',
                    api_key=panorama.PANORA_PUBLIC_KEY
                )

            # verify the URL passed to get()
            mock_get.assert_called_once_with(
                PUBLIC_FILE,
                headers={'Authorization': f'Basic {panorama._encode_api_key(panorama.PANORA_PUBLIC_KEY)}'},
                stream=True, verify=True
            )

        else:
            if not panorama.have_internet():
                self.skipTest("Internet connection is required for this test.")

            out_path = panorama.download_webdav_file(
                PUBLIC_FILE,
                dest_path=f'{self.work_dir}/{test_file}',
                api_key=panorama.PANORA_PUBLIC_KEY
            )

        # ensure file is written
        self.assertTrue(os.path.isfile(out_path))
        with open(out_path, 'rb') as outF:
            self.assertEqual(outF.read(), payload)


    def test_download_anonymous(self):
        test_file = 'Clean-SMTG-B1-1410-Oct2022-3qtrans-Meta.csv'
        with open(f'{self.data_dir}/validate_pipeline_params/panorama/{test_file}', 'rb') as inF:
            payload = inF.read()

        if MOCK_PANORAMA:
            with mock.patch('DIA_QC_report.submodules.panorama.requests.get',
                    return_value=self._fake_response(data=payload)) as mock_get:
                out_path = panorama.download_webdav_file(
                    PUBLIC_FILE, dest_path=f'{self.work_dir}/{test_file}'
                )

            # verify the URL passed to get()
            mock_get.assert_called_once_with(
                PUBLIC_FILE, headers={}, stream=True, verify=True
            )
        else:
            if not panorama.have_internet():
                self.skipTest("Internet connection is required for this test.")

            out_path = panorama.download_webdav_file(
                PUBLIC_FILE, dest_path=f'{self.work_dir}/{test_file}'
            )

        # ensure file is written
        self.assertTrue(os.path.isfile(out_path))
        with open(out_path, 'rb') as outF:
            self.assertEqual(outF.read(), payload)


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

        # verify the URL passed to get()
        mock_get.assert_called_once_with(
            f'https://server/_webdav/lab/%40files/{test_file}',
            headers={}, stream=True, verify=True
        )

        target_text = payload.decode('utf-8')
        self.assertEqual(file_text, target_text)