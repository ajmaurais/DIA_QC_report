
import unittest
from unittest import mock
from types import SimpleNamespace
import io

from DIA_QC_report.submodules import nextflow_pipeline_config as npc
import setup_functions


class TestGToPy(unittest.TestCase):
    def test_plain_string(self):
        self.assertEqual(
            npc._g_to_py("/home/user/data/file.txt"),
            "/home/user/data/file.txt",
        )
        self.assertEqual(
            npc._g_to_py("https://example.com/path"),
            "https://example.com/path",
        )


    def test_quoted_string(self):
        self.assertEqual(npc._g_to_py("'hello world'"), "hello world")
        self.assertEqual(npc._g_to_py('"hi"'), "hi")


    def test_groovy_list(self):
        self.assertEqual(npc._g_to_py("['x', 'y']"), ["x", "y"])
        self.assertEqual(npc._g_to_py("['https://aurl.com', 'y']"), ["https://aurl.com", "y"])
        self.assertEqual(npc._g_to_py("['https://example.com/raw/', '/local/dir']"),
                                                           ['https://example.com/raw/', '/local/dir'])



    def test_groovy_map(self):
        self.assertEqual(
            npc._g_to_py("[foo:'bar', baz: 'qux']"),
            {"foo": "bar", "baz": "qux"},
        )
        self.assertEqual(
            npc._g_to_py("['k 1':'v1', 'k2':'v2']"),
            {"k 1": "v1", "k2": "v2"},
        )
        self.assertEqual(
            npc._g_to_py("['k 1': 1, 'k2':false]"),
            {"k 1": 1, "k2": False},
        )
        self.assertEqual(
            npc._g_to_py("['k1': 'https://aurl.com', 'k2':'/local/path']"),
            {"k1": "https://aurl.com", "k2": "/local/path"}
        )


    def test_booleans(self):
        self.assertIs(npc._g_to_py("true"), True)
        self.assertIs(npc._g_to_py("false"), False)


    def test_null_literal(self):
        self.assertIsNone(npc._g_to_py("null"))


    def test_numbers(self):
        self.assertEqual(npc._g_to_py("42"), 42)
        self.assertAlmostEqual(npc._g_to_py("3.14"), 3.14)


class TestFindParamsBlock(unittest.TestCase):
    def test_empty(self):
        tree = npc.parse_and_digest_groovy_content("")
        with self.assertRaises(ValueError) as e:
            npc._find_params_block(tree)
        self.assertIn('No top-level params { ... } block found in the config', str(e.exception))


    def test_no_params(self):
        tree = npc.parse_and_digest_groovy_content('not_a_params_block { hello = "world" }')
        with self.assertRaises(ValueError) as e:
            npc._find_params_block(tree)
        self.assertIn('No top-level params { ... } block found in the config', str(e.exception))


    def test_params_block(self):
        config = '''
            params {
                input = "data.txt"
                output = "results.txt"
            }

            docker {
                enabled = true
            }
            includeConfig '/path/to/another/pipeline.config'
        '''
        tree = npc.parse_and_digest_groovy_content(config)
        params_tree = npc._find_params_block(tree)
        self.assertIsInstance(params_tree, dict)

        data = npc.parse_params(text=config)
        self.assertIsInstance(data, SimpleNamespace)
        self.assertTrue(hasattr(data, 'input'))
        self.assertTrue(hasattr(data, 'output'))
        self.assertEqual(data.input, 'data.txt')
        self.assertEqual(data.output, 'results.txt')
        self.assertEqual(len(vars(data)), 2)


class TestNodeHasIdentifier(unittest.TestCase):
    def test_node_has_identifier(self):
        node = {
            'leaf': 'ASSIGN',
            'children': [
                {'leaf': 'IDENTIFIER', 'value': 'foo'},
                {'leaf': 'STRING', 'value': 'bar'}
            ]
        }
        self.assertTrue(npc._node_has_ident(node, 'foo'))
        self.assertFalse(npc._node_has_ident(node, 'baz'))


    def test_node_does_not_have_identifier(self):
        node = {
            'leaf': 'ASSIGN',
            'children': [
                {'leaf': 'IDENTIFIER', 'value': 'foo'},
                {'leaf': 'STRING', 'value': 'bar'}
            ]
        }
        self.assertFalse(npc._node_has_ident(node, 'baz'))


class TestParseParams(unittest.TestCase):
    def test_parse_params(self):
        config = '''
            // The params section
            params {
                input = "data.txt"
                output = "results.txt"
                threshold = 0.05
                a_boolean = true
                a_null_value = null
                n_files = 10
            }
        '''
        with mock.patch('pathlib.Path.read_text', return_value=config):
            params = npc.parse_params('pipeline.config')
            self.assertEqual(params.input, 'data.txt')
            self.assertEqual(params.output, 'results.txt')
            self.assertEqual(params.threshold, 0.05)
            self.assertEqual(params.a_boolean, True)
            self.assertIsNone(params.a_null_value)
            self.assertEqual(params.n_files, 10)


    def test_parse_complex_types(self):
        url_str = 'https://panoramaweb.org/_webdav/path/to/dir'
        config = '''
            /* The params
               section */
            params {
                string_param = "%s"
                list = ["/local/dir", "%s"]
                map = [key1: "/local/dir", key2: "%s"]
            } ''' % (url_str, url_str, url_str)

        with mock.patch('pathlib.Path.read_text', return_value=config):
            params = npc.parse_params('pipeline.config')
            self.assertEqual(params.string_param, url_str)
            self.assertEqual(params.list, ["/local/dir", url_str])
            self.assertEqual(params.map, {'key1': '/local/dir', 'key2': url_str})


    def test_multi_line_string(self):
        url_str = 'https://panoramaweb.org/_webdav/path/to/dir'
        config = """
            params {
                multi_line = '''
                    /dir/one
                    %s
                    '''
            } """ % url_str

        params = npc.parse_params(text=config)
        self.assertEqual(
            npc.param_to_list(params.multi_line),
            ['/dir/one', url_str]
        )


    def test_tripple_quoted_string(self):
        target = SimpleNamespace(tripple_quoted='text')

        config_single = "params { tripple_quoted = '''text''' }"
        config_data = npc.parse_params(text=config_single)
        self.assertEqual(config_data, target)

        config_double = 'params { tripple_quoted = """text""" }'
        config_data = npc.parse_params(text=config_double)
        self.assertEqual(config_data, target)


    def test_parse_string_list(self):
        config = '''
            params {
                string_list = ['https://example.com/raw/', '/local/dir']
            } '''
        target = SimpleNamespace(string_list=['https://example.com/raw/', '/local/dir'])

        config_data = npc.parse_params(text=config)
        self.assertEqual(config_data, target)


    def test_parse_map(self):
        config = '''
            params {
                map = [ key1: 'https://example.com/raw/', key2: '/local/dir' ]
            } '''
        target = SimpleNamespace(map={'key1': 'https://example.com/raw/', 'key2': '/local/dir'})

        config_data = npc.parse_params(text=config)
        self.assertEqual(config_data, target)


    def test_parse_nested_identifiers(self):
        config = '''
            params {
                map = [ key1: 'https://example.com/raw/', key2: '/local/dir' ]
                list = ['https://example.com/raw/', '/local/dir']
                carafe {
                    spectra_file         = '/path/to/spectra.mzML'
                    peptide_results_file = 'results.tsv'
                }
            } '''
        target = SimpleNamespace(
            list=['https://example.com/raw/', '/local/dir'],
            map={'key1': 'https://example.com/raw/', 'key2': '/local/dir'},
            carafe=SimpleNamespace(spectra_file='/path/to/spectra.mzML',
                                   peptide_results_file='results.tsv')
        )
        config_data = npc.parse_params(text=config)
        self.assertEqual(config_data, target)


    def test_dot_and_brace_nesting_identical(self):
        brace_config = '''
            params {
                carafe {
                    spectra_file         = '/path/to/spectra.mzML'
                    peptide_results_file = 'results.tsv'
                }
            } '''

        dot_config = '''
            params {
                carafe.peptide_results_file = 'results.tsv'
                carafe.spectra_file         = '/path/to/spectra.mzML'
            } '''

        target = SimpleNamespace(
            carafe=SimpleNamespace(peptide_results_file='results.tsv', spectra_file='/path/to/spectra.mzML')
        )
        brace_data = npc.parse_params(text=brace_config)
        dot_data = npc.parse_params(text=dot_config)
        self.assertEqual(brace_data, dot_data)
        self.assertEqual(brace_data, target)


    def test_complex_config(self):
        config_file = f'{setup_functions.TEST_DIR}/data/validate_pipeline_params/config_files/panorama_multi_batch_diann.config'

        target = SimpleNamespace(
            carafe=SimpleNamespace(
                spectra_file = "<path>",
                peptide_results_file = None,
                carafe_fasta = None,
                diann_fasta = None
            ),
            quant_spectra_dir = {'SP3': '<path>', 'Strap': '<path>'},
            quant_spectra_glob = '<glob>',
            quant_spectra_regex = None,
            chromatogram_library_spectra_dir = ['<path>', '<path>'],
            chromatogram_library_spectra_glob = '<glob>',
            chromatogram_library_spectra_regex = None,
            images = SimpleNamespace(
                diann = "quay.io/mauraisa/diann:1.8.1",
                proteowizard = 'quay.io/mauraisa/pwiz-skyline-i-agree-to-the-vendor-licenses:imputation_special_build_5.25'
            ),
            spectral_library = '<path>',
            fasta = '<path>',
            search_engine = 'diann',
            msconvert = SimpleNamespace(do_demultiplex = True, do_simasspectra = True),
            encyclopedia = SimpleNamespace(
                chromatogram = SimpleNamespace(params=None),
                quant = SimpleNamespace(params=None),
                save_output = False
            ),
            replicate_metadata = '<path>',
            qc_report = SimpleNamespace(
                skip = False,
                color_vars = ['cellLine'],
                export_tables = False
            ),
            batch_report = SimpleNamespace(skip=False),
            skyline = SimpleNamespace(
                minimize = False,
                group_by_gene = True,
                protein_parsimony = True,
                document_name = 'US_UW',
                use_hardlinks = True,
                skyr_file = '<path>',
                skip = False
            ),
            panorama = SimpleNamespace(upload = False, upload_url = '<path>', import_skyline = False)
        )

        config_data = npc.parse_params(config_file)
        self.assertEqual(config_data, target)


    def test_write_params(self):
        config_file = f'{setup_functions.TEST_DIR}/data/validate_pipeline_params/config_files/panorama_multi_batch_diann.config'
        config_data = npc.parse_params(config_file)

        out = io.StringIO()
        npc.write_params_namespace(SimpleNamespace(params=config_data), file=out)
        out.seek(0)
        written_config = out.read()
        written_data = npc.parse_params(text=written_config)
        self.assertEqual(config_data, written_data)


    def test_namespace_to_dict(self):
        config = '''
            params {
                map = [ key1: 'https://example.com/raw/', key2: '/local/dir' ]
                list = ['https://example.com/raw/', '/local/dir']
                carafe {
                    spectra_file         = '/path/to/spectra.mzML'
                    peptide_results_file = 'results.tsv'
                }
            } '''
        target = {
            'map': { 'key1': 'https://example.com/raw/', 'key2': '/local/dir' },
            'list': ['https://example.com/raw/', '/local/dir'],
            'carafe': {
                'spectra_file'         : '/path/to/spectra.mzML',
                'peptide_results_file' : 'results.tsv'
            }
        }
        config_data = npc.parse_params(text=config)
        config_dict = npc.namespace_to_dict(config_data)
        self.assertEqual(config_dict, target)