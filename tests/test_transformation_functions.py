
import unittest

import setup_functions

import DIA_QC_report.submodules.transformation as transformation


VALID_OPTIONS = {'weights':      {'default': 'uniform',
                                  'choices': ('uniform', 'weighted'),
                                  'dtype': str},
                 'n_neighbors':  {'default': 5,
                                  'dtype': int,
                                  'min_value': 0,
                                  'min_inclusive': False},
                 'max_missing': {'dtype': float,
                                 'default': 0.5,
                                 'min_value': 0,
                                 'max_value': 1,
                                 'max_inclusive': True},
                 'min_i':       {'dtype': int,
                                 'min_value': 0,
                                 'min_inclusive': True},
                 'max_i':       {'dtype': int,
                                 'max_value': 10,
                                 'max_inclusive': False},
                 'bool_test':   {'dtype': bool,
                                 'default': True}}


def get_err_message(opt, err_type, value):
    if err_type == 'value':
        return f"Cannot convert value <{value}> to {str(opt.dtype)}"
    elif err_type == 'choice':
        return f"Value '{value}' must be in {str(opt.choices)}"
    elif err_type == 'range':
        return opt._format_range_err()
    else:
        raise RuntimeError('Unknown err_type!')


class TestCammelCase(unittest.TestCase):
    def test(self):
        self.assertEqual(transformation.cammel_case('log2', 'area'), 'log2Area')
        norm_value = transformation.cammel_case('normalized', 'area')
        self.assertEqual(norm_value, 'normalizedArea')
        self.assertEqual(transformation.cammel_case('log2', norm_value), 'log2NormalizedArea')


class TestOption(unittest.TestCase):
    def test_invalid_init(self):
        init_args = [['invalid_type', {'choices': ('a', 'b'), 'dtype': int}],
                     ['invalid_min', {'choices': (1, 2), 'min_value': 1}],
                     ['invalid_max', {'choices': (1, 2), 'max_value': 1}],
                     ['invalid_max', {'dtype': str, 'max_value': 1}],
                     ['invalid_default', {'default': 'a'}],
                     ['invalid_type', {'choices': (1, 2), 'dtype': str}]]

        error_msgs = {'invalid_type': r"Choice '[a-zA-Z\-_0-9]+' is not of type <class '[a-zA-Z]+'>",
                      'invalid_min': 'Cannot set min/max value for option `dtype` or `choices`!',
                      'invalid_max': 'Cannot set min/max value for option `dtype` or `choices`!',
                      'invalid_default': r"Invalid default argument: '[a-zA-Z0-9\_-]+'!"}

        for name, kwargs in init_args:
            self.assertRaisesRegex(ValueError, error_msgs[name],
                                   transformation.Option, name, **kwargs)


    def test_valid_choice(self):
        tests = {'max_missing': [(0.5, True), (9, True), (0, True), (1, True)],
                 'n_neighbors': [(5, True), (-1, True), (0, True), (400, True)],
                 'min_i': [(0, True), (-1, True), (1000, True)],
                 'bool_test': [('dummy', True), (True, True), (True, True)],
                 'weights': [('uniform', True), ('weighted', True), ('dummy', False)]}

        for name, args in tests.items():
            o = transformation.Option(name, **VALID_OPTIONS[name])
            for value, target in args:
                self.assertEqual(o.valid_choice(value), target)


    def test_can_convert(self):
        tests = {'max_missing': [(0.3, True), ('5', True), ('0.5', True), (1, True)],
                 'n_neighbors': [('5', True), (-1, True), ('abc', False), ('1a', False)],
                 'min_i': [(0, True), ('-1', True), ('0.1', False), ('1.1', False), ('h7', False)],
                 'bool_test': [('dummy', False), ('true', True), ('on', True), ('FALSE', True)],
                 'weights': [('uniform', True), (7, True), ('dummy', True)]}

        for name, args in tests.items():
            o = transformation.Option(name, **VALID_OPTIONS[name])
            for value, target in args:
                self.assertEqual(o.can_convert(value), target)


    def test_in_range(self):
        tests = {'max_missing': [(0.5, True), (9, False), (0, False), (1, True)],
                 'n_neighbors': [(5, True), (-1, False), (0, False), (400, True)],
                 'min_i': [(0, True), (-1, False), (1000, True)],
                 'bool_test': [('dummy', True), (True, True), (False, True)],
                 'weights': [('uniform', True), ('weighted', True), ('dummy', True)]}

        for name, args in tests.items():
            o = transformation.Option(name, **VALID_OPTIONS[name])
            for value, target in args:
                self.assertEqual(o.in_range(value), target)


    def test_format_range_err(self):
        messages = {'max_missing': 'Value must be > 0 and <= 1!',
                    'n_neighbors': 'Value must be > 0!',
                    'min_i': 'Value must be >= 0!',
                    'max_i': 'Value must be < 10!',
                    'bool_test': '',
                    'weights': ''}

        for name, msg in messages.items():
            o = transformation.Option(name, **VALID_OPTIONS[name])
            self.assertEqual(o._format_range_err(), msg)


    def test_parse(self):
        good_tests = {'max_missing': [(0.5, 0.5), ('0.9', 0.9), ('0.25', 0.25), ('.5', 0.5)],
                      'n_neighbors': [('5', 5), (5, 5), ('100', 100), ('010', 10)],
                      'min_i': [(0, 0), ('1', 1), ('05', 5), ('7', 7)],
                      'bool_test': [('true', True), ('on', True), ('FALSE', False), (False, False)],
                      'weights': [('uniform', 'uniform'), ('weighted', 'weighted')]}

        for name, args in good_tests.items():
            o = transformation.Option(name, **VALID_OPTIONS[name])
            for value, target in args:
                self.assertEqual(o.parse(value), target)

        bad_tests = {'max_missing': [('2', 'range'), (2, 'range'), ('e4', 'value')],
                     'n_neighbors': [(-1, 'range'), ('0', 'range'), ('e5', 'value'), ('0.1', 'value')],
                     'min_i': [('-1', 'range'), ('ke1', 'value')],
                     'max_i': [(11, 'range'), ('12', 'range'), ('ke7', 'value')],
                     'bool_test': [('dummy', 'value')],
                     'weights': [('dummy', 'choice')]}

        for name, args in bad_tests.items():
            o = transformation.Option(name, **VALID_OPTIONS[name])
            for arg, err_type in args:
                with self.assertNoLogs(transformation.LOGGER):
                    ret = o.parse(arg, quiet=True)
                self.assertIsNone(ret)

                with self.assertLogs(transformation.LOGGER, level='ERROR') as cm:
                    ret = o.parse(arg)
                self.assertIsNone(ret)
                target_msg = get_err_message(o, err_type, arg)
                self.assertEqual(len(cm.output), 1)
                self.assertTrue(target_msg in l for l in cm.output[0])


class TestMethodOptionsRegex(unittest.TestCase):
    def test_option_regex_valid(self):
        test_strs = {'a=b': ('a', 'b'),
                     'a = b': ('a', 'b'),
                     'a= b': ('a', 'b'),
                     'a =b': ('a', 'b'),
                     'a b': ('a', 'b'),
                     'abc 123-456': ('abc', '123-456'),
                     "abc '123 456'": ('abc', '123 456'),
                     "abc 0.56": ('abc', '0.56'),
                     "abc = 0.56": ('abc', '0.56'),
                     "abc=0.56": ('abc', '0.56'),
                     ' a b': ('a', 'b'),
                     ' a  b ': ('a', 'b'),
                     '\t a  b ': ('a', 'b')}

        for s, (k, v) in test_strs.items():
            match = transformation.MethodOptions.OPTION_RE.search(s)
            self.assertIsNotNone(match)

            name, value = transformation.MethodOptions._get_re_key_value_str(match)
            self.assertEqual(name, k)
            self.assertEqual(value, v)


class TestMethodOptionsParseOption(unittest.TestCase):
    def setUp(self):
        self.options = transformation.MethodOptions()

        for name, kwargs in VALID_OPTIONS.items():
            self.options.add_option(name, **kwargs)


    def test_parse_options_valid(self):
        test_args = ['weights=uniform', 'weights=weighted',
                     'weights = uniform', "weights 'weighted'", "weights='weighted'",
                     'n_neighbors=6', 'n_neighbors "5"', ' n_neighbors 7',
                     'max_missing=0.1', 'max_missing 0.5']

        with self.assertNoLogs(transformation.LOGGER, level='ERROR'):
            self.assertTrue(self.options.parse_strings(test_args))

        final_values = {'weights': 'weighted', 'n_neighbors': 7, 'max_missing': 0.5}
        for name, value in final_values.items():
            self.assertEqual(self.options.values[name], value)


    def test_parse_options_invalid(self):
        test_args = [('dummy', "Could not parse option argument: 'dummy'"),
                     ('max_missing =', "Could not parse option argument: 'max_missing ='"),
                     ('dummy on', "Unknown option: 'dummy'"),
                     ('max_missing=0', get_err_message(self.options.options['max_missing'], 'range', 0)),
                     ('weights dummy', get_err_message(self.options.options['weights'], 'choice', 'dummy'))]

        with self.assertLogs(transformation.LOGGER, level='ERROR') as cm:
            self.assertFalse(self.options.parse_strings([a[0] for a in test_args]))

        for i, (_, msg) in enumerate(test_args):
            self.assertTrue(msg in cm.output[i])
