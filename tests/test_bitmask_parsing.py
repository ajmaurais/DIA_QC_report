
import unittest
from unittest import mock
from itertools import product

from DIA_QC_report.submodules.dia_db_utils import validate_bit_mask
from DIA_QC_report.submodules.dia_db_utils import parse_bitmask_options
from DIA_QC_report.submodules.dia_db_utils import LOGGER as db_utils_logger

class TestBitMaskParsing(unittest.TestCase):

    def test_invalid_mask(self):
        masks = [('89', 3, 2),
                 ('45', 1, 2),
                 ('45', 2, 2)]
        max_values = (0, 1, 3, 7)

        for mask, n_options, n_digits in masks:
            with self.assertLogs(db_utils_logger, level='ERROR') as cm:
                result = validate_bit_mask(mask, n_options=n_options, n_digits=n_digits)

            self.assertFalse(result)

            self.assertEqual(len(cm.output), 1)
            expected_message = f'Bit mask digits must be between 0 and {max_values[n_options]}!'
            self.assertTrue(expected_message in cm.output[0])


    def test_invalid_n_digits(self):
        masks = [('12', 3, 1),
                 ('45', 3, 3)]

        for mask, n_options, n_digits in masks:
            with self.assertLogs(db_utils_logger, level='ERROR') as cm:
                result = validate_bit_mask(mask, n_options=n_options, n_digits=n_digits)

            self.assertFalse(result)
            self.assertEqual(len(cm.output), 1)
            self.assertTrue(f'Bit mask must be {n_digits} digits!', cm.output[0])


    def test_invalid_n_options(self):
        masks = [('00', 4, 2),
                 ('00', 0, 2)]

        for mask, n_options, n_digits in masks:
            with self.assertRaises(ValueError):
                validate_bit_mask(mask, n_options=n_options, n_digits=n_digits)


    @staticmethod
    def permute_options(digit_names, option_map, method_names):
        std = {}
        for digits in product(*[range(len(option_map)) for _ in range(len(digit_names))]):
            key = ''.join(str(x) for x in digits)
            std[key] = {}
            for i in range(len(digit_names)):
                std[key][digit_names[i]] = dict(zip(method_names, option_map[digits[i]]))

        return std


    def test_valid_1_options(self):
        # setup ground truth
        option_map = [(False,), # 0
                      (True,)]  # 1

        digit_names = ('wide', 'long', 'extra_long')
        method_names = ('write',)

        for i in range(len(digit_names)):
            digit_slice = digit_names[0:i+1]
            std = self.permute_options(digit_slice, option_map, method_names)
            for option, result in std.items():
                test_result = parse_bitmask_options(option, digit_slice, method_names)
                for digit in digit_slice:
                    self.assertDictEqual(result[digit], test_result[digit])


    def test_valid_2_options(self):
        # setup ground truth
        option_map = [(False, False), # 0
                      (True, False),  # 1
                      (False, True),  # 2
                      (True, True)]   # 3

        method_names = ('unnormalized', 'normalized')
        digit_names = ('wide', 'long', 'extra_long')

        for i in range(len(digit_names)):
            digit_slice = digit_names[0:i+1]
            std = self.permute_options(digit_slice, option_map, method_names)
            for option, result in std.items():
                test_result = parse_bitmask_options(option, digit_slice, method_names)
                for digit in digit_slice:
                    self.assertDictEqual(result[digit], test_result[digit])


    def test_valid_3_options(self):
        # setup ground truth
        option_map = [(False, False, False), # 0
                      (True, False, False),  # 1
                      (False, True, False),  # 2
                      (True, True, False),   # 3
                      (False, False, True),  # 4
                      (True, False, True),   # 5
                      (False, True, True),   # 6
                      (True, True, True)]    # 7

        digit_names = ('wide', 'long', 'extra_long')
        method_names = ('unnormalized', 'normalized', 'batch_corrected')

        for i in range(len(digit_names)):
            digit_slice = digit_names[0:i+1]
            std = self.permute_options(digit_slice, option_map, method_names)
            for option, result in std.items():
                test_result = parse_bitmask_options(option, digit_slice, method_names)
                for digit in digit_slice:
                    self.assertDictEqual(result[digit], test_result[digit])


if __name__ == '__main__':
    unittest.main()
