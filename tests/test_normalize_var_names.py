
import unittest

from DIA_QC_report.generate_qc_qmd import normalize_var_names


class TestNormalizeVarNames(unittest.TestCase):
    def test_invalid_chars(self):
        meta_keys = {'1var1': '_var1', '123var2': '_var2',
                     '1var123': '_var123', '123va1r': '_va1r',
                     'var123': 'var123', 'v1ar': 'v1ar',
                     'v$a': 'v_a', 'v#.k': 'v_k', 'abc.123': 'abc_123'}

        test_keys = normalize_var_names(meta_keys.keys())

        self.assertEqual(len(test_keys), len(meta_keys))
        for var, expected_fixed_var in meta_keys.items():
            self.assertEqual(expected_fixed_var, test_keys[var])


    def test_duplicate_vars(self):
        meta_keyss = [{'1va1r': '_va1r_1', '12va1r': '_va1r_2'},
                      {'1var': '_var_1', '123var': '_var_2', 'var123': 'var123', 'v1ar': 'v1ar',
                       '1var123': '_var123_1', '1va1r': '_va1r_1',
                       '12var123': '_var123_2', '12va1r': '_va1r_2'},
                      {'abc@#_()*': 'abc____1', 'abc___': 'abc____2'}]

        for key_set in meta_keyss:
            test_keys = normalize_var_names(key_set)

            self.assertEqual(len(key_set), len(test_keys))
            for var, expected_fixed_var in key_set.items():
                self.assertEqual(expected_fixed_var, test_keys[var])
