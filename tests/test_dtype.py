
import unittest

from DIA_QC_report.submodules.dtype import Dtype

class TestDtype(unittest.TestCase):
    @classmethod
    def setUpClass(TestDtype):
        TestDtype.na_strings = ['NA', 'NULL', '#N/A', 'NaN']

        TestDtype.bool_strings = ['true', 'True', 'TRUE', 'False', 'false', 'FALSE']
        TestDtype.bool_bools = [True, True, True, False, False, False]

        TestDtype.int_strings = ['1', '01', '+8', '-8', '10029', '-69', '-01747', '+98133']
        TestDtype.int_ints = [1, 1, 8, -8, 10029, -69, -1747, 98133]

        TestDtype.float_strings = ['0.1', '+0.102', '-0.15', '1.', '1.4e5', '3.14E0', '7.5e-2']
        TestDtype.float_floats = [0.1, 0.102, -0.15, 1.0, 1.4e5, 3.14e0, 7.5e-2]

        TestDtype.string_strings = ['0.5r', 'op90']


    def test_infer_type_NA(self):
        for s in self.na_strings:
            self.assertIs(Dtype.infer_type(s), Dtype.NULL)


    def test_infer_type_bool(self):
        for s in self.bool_strings:
            self.assertIs(Dtype.infer_type(s), Dtype.BOOL)


    def test_infer_type_int(self):
        for s in self.int_strings:
            self.assertIs(Dtype.infer_type(s), Dtype.INT)


    def test_infer_type_float(self):
        for s in self.float_strings:
            self.assertIs(Dtype.infer_type(s), Dtype.FLOAT)


    def test_infer_type_string(self):
        for s in self.string_strings:
            self.assertIs(Dtype.infer_type(s), Dtype.STRING)


    def test_invalid_bool_converstions_fail(self):
        for s in self.float_strings + self.string_strings + self.na_strings:
            self.assertRaises(ValueError, Dtype.BOOL.convert, s)


    def test_invalid_int_converstions_fail(self):
        for s in self.string_strings:
            self.assertRaises(ValueError, Dtype.INT.convert, s)


    def test_invalid_float_converstions_fail(self):
        for s in self.string_strings:
            self.assertRaises(ValueError, Dtype.FLOAT.convert, s)


    def test_na_conversions(self):
        for n in self.na_strings:
            self.assertIs(Dtype.convert(Dtype.NULL, n), None)


    def test_int_conversions(self):
        for s, i in zip(self.int_strings, self.int_ints):
            self.assertEqual(Dtype.convert(Dtype.INT, s), i)


    def test_float_conversions(self):
        for s, d in zip(self.float_strings, self.float_floats):
            self.assertEqual(Dtype.convert(Dtype.FLOAT, s), d)


    def test_bool_conversions(self):
        for s, b in zip(self.bool_strings, self.bool_bools):
            self.assertEqual(Dtype.convert(Dtype.BOOL, s), b)

        self.assertEqual(Dtype.convert(Dtype.BOOL, '1'), True)
        self.assertEqual(Dtype.convert(Dtype.BOOL, '0'), False)


    def test_comparison(self):
        self.assertTrue(Dtype.INT >= Dtype.BOOL)
        self.assertTrue(Dtype.INT >= Dtype.NULL)
        self.assertTrue(Dtype.INT >= Dtype.INT)
        self.assertFalse(Dtype.INT >= Dtype.FLOAT)
        self.assertFalse(Dtype.INT >= Dtype.STRING)

        self.assertFalse(Dtype.INT < Dtype.BOOL)
        self.assertFalse(Dtype.INT < Dtype.NULL)
        self.assertFalse(Dtype.INT < Dtype.INT)
        self.assertTrue(Dtype.INT < Dtype.FLOAT)
        self.assertTrue(Dtype.INT < Dtype.STRING)

        self.assertRaises(ValueError, Dtype.FLOAT.__lt__, (Dtype.FLOAT, 0.1))
        self.assertRaises(ValueError, Dtype.FLOAT.__ge__, (Dtype.FLOAT, 0.1))


    def test_var_to_type(self):
        self.assertRaises(ValueError, Dtype.var_to_type, ['a', 'list'])
        self.assertRaises(ValueError, Dtype.var_to_type, {'a': 'dict'})
        self.assertRaises(ValueError, Dtype.var_to_type, {'a', 'set'})

        self.assertIs(Dtype.var_to_type(None), Dtype.NULL)
        self.assertIs(Dtype.var_to_type(True), Dtype.BOOL)
        self.assertIs(Dtype.var_to_type(False), Dtype.BOOL)
        self.assertIs(Dtype.var_to_type(69), Dtype.INT)
        self.assertIs(Dtype.var_to_type(3.14), Dtype.FLOAT)
        self.assertIs(Dtype.var_to_type('Hello there General Kenobi'), Dtype.STRING)


if __name__ == '__main__':
    unittest.main()
