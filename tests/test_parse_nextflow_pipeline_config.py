
import unittest
from unittest import mock

from DIA_QC_report.submodules import nextflow_pipeline_config


class TestGToPy(unittest.TestCase):
    def test_plain_string(self):
        self.assertEqual(
            nextflow_pipeline_config._g_to_py("/home/user/data/file.txt"),
            "/home/user/data/file.txt",
        )
        self.assertEqual(
            nextflow_pipeline_config._g_to_py("https://example.com/path"),
            "https://example.com/path",
        )


    def test_quoted_string(self):
        self.assertEqual(nextflow_pipeline_config._g_to_py("'hello world'"), "hello world")
        self.assertEqual(nextflow_pipeline_config._g_to_py('"hi"'), "hi")


    def test_brace_list(self):
        self.assertEqual(nextflow_pipeline_config._g_to_py("{a,b,c}"), ["a", "b", "c"])


    def test_groovy_list(self):
        self.assertEqual(nextflow_pipeline_config._g_to_py("['x', 'y']"), ["x", "y"])
        self.assertEqual(nextflow_pipeline_config._g_to_py("['https://aurl.com', 'y']"), ["https://aurl.com", "y"])


    def test_groovy_map(self):
        self.assertEqual(
            nextflow_pipeline_config._g_to_py("[foo:'bar', baz: 'qux']"),
            {"foo": "bar", "baz": "qux"},
        )
        self.assertEqual(
            nextflow_pipeline_config._g_to_py("['k 1':'v1', 'k2':'v2']"),
            {"k 1": "v1", "k2": "v2"},
        )
        self.assertEqual(
            nextflow_pipeline_config._g_to_py("['k 1': 1, 'k2':false]"),
            {"k 1": 1, "k2": False},
        )
        self.assertEqual(
            nextflow_pipeline_config._g_to_py("['k1': 'https://aurl.com', 'k2':'/local/path']"),
            {"k1": "https://aurl.com", "k2": "/local/path"}
        )


    def test_booleans(self):
        self.assertIs(nextflow_pipeline_config._g_to_py("true"), True)
        self.assertIs(nextflow_pipeline_config._g_to_py("false"), False)


    def test_null_literal(self):
        self.assertIsNone(nextflow_pipeline_config._g_to_py("null"))


    def test_numbers(self):
        self.assertEqual(nextflow_pipeline_config._g_to_py("42"), 42)
        self.assertAlmostEqual(nextflow_pipeline_config._g_to_py("3.14"), 3.14)


class TestFindParamsBlock(unittest.TestCase):
    def test_empty(self):
        tree = nextflow_pipeline_config.parse_and_digest_groovy_content("")
        with self.assertRaises(ValueError) as e:
            nextflow_pipeline_config._find_params_block(tree)
        self.assertIn('No `params { ... }` block found in the config.', str(e.exception))


    def test_no_params(self):
        tree = nextflow_pipeline_config.parse_and_digest_groovy_content('not_a_params_block { hello = "world" }')
        with self.assertRaises(ValueError) as e:
            nextflow_pipeline_config._find_params_block(tree)
        self.assertIn('No `params { ... }` block found in the config.', str(e.exception))


    def test_params_block(self):
        config = '''
            params {
                input = "data.txt"
                output = "results.txt"
            }
             includeConfig '/path/to/another/pipeline.config'
        '''
        tree = nextflow_pipeline_config.parse_and_digest_groovy_content(config)
        params_tree = nextflow_pipeline_config._find_params_block(tree)
        self.assertIsInstance(params_tree, dict)


class TestNodeHasIdentifier(unittest.TestCase):
    def test_node_has_identifier(self):
        node = {
            'leaf': 'ASSIGN',
            'children': [
                {'leaf': 'IDENTIFIER', 'value': 'foo'},
                {'leaf': 'STRING', 'value': 'bar'}
            ]
        }
        self.assertTrue(nextflow_pipeline_config._node_has_ident(node, 'foo'))
        self.assertFalse(nextflow_pipeline_config._node_has_ident(node, 'baz'))


    def test_node_does_not_have_identifier(self):
        node = {
            'leaf': 'ASSIGN',
            'children': [
                {'leaf': 'IDENTIFIER', 'value': 'foo'},
                {'leaf': 'STRING', 'value': 'bar'}
            ]
        }
        self.assertFalse(nextflow_pipeline_config._node_has_ident(node, 'baz'))


class TestContainsLeaf(unittest.TestCase):
    def test_contains_leaf(self):
        node = {
            'leaf': 'ASSIGN',
            'children': [
                {'leaf': 'IDENTIFIER', 'value': 'foo'},
                {'leaf': 'STRING', 'value': 'bar'}
            ]
        }
        self.assertTrue(nextflow_pipeline_config._contains_leaf(node, 'IDENTIFIER'))
        self.assertFalse(nextflow_pipeline_config._contains_leaf(node, 'NON_EXISTENT_LEAF'))


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
            params = nextflow_pipeline_config.parse_params('pipeline.config')
            self.assertEqual(params['input'], 'data.txt')
            self.assertEqual(params['output'], 'results.txt')
            self.assertEqual(params['threshold'], 0.05)
            self.assertEqual(params['a_boolean'], True)
            self.assertIsNone(params['a_null_value'])
            self.assertEqual(params['n_files'], 10)


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
            params = nextflow_pipeline_config.parse_params('pipeline.config')
            self.assertEqual(params['string_param'], url_str)
            self.assertEqual(params['list'], ["/local/dir", url_str])
            self.assertEqual(params['map'], {'key1': '/local/dir', 'key2': url_str})


    def test_multi_line_string(self):
        url_str = 'https://panoramaweb.org/_webdav/path/to/dir'
        config = """
            params {
                multi_line = '''
                    /dir/one
                    %s
                    '''
            } """ % url_str

        with mock.patch('pathlib.Path.read_text', return_value=config):
            params = nextflow_pipeline_config.parse_params('pipeline.config')
            self.assertEqual(
                nextflow_pipeline_config.param_to_list(params['multi_line']),
                ['/dir/one', url_str]
            )