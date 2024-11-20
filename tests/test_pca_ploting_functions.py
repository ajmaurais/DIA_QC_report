
import unittest

from pandas import NA, isna

from DIA_QC_report.submodules.pca_plot import map_discrete_colors
from DIA_QC_report.submodules.pca_plot import GREY_RGB
from DIA_QC_report.submodules.pca_plot import NA_CATEGORY_NAME


class TestMapDiscreteColors(unittest.TestCase):

    def do_no_nas_test(self, categories):
        maped_categories, colors = map_discrete_colors(categories)

        self.assertEqual(len(maped_categories), len(categories))
        self.assertEqual(len(colors), len(set(maped_categories.values())))
        self.assertFalse(GREY_RGB in colors.values())


    def do_na_test(self, categories):
        maped_categories, colors = map_discrete_colors(categories)

        self.assertEqual(len(maped_categories), len(categories))
        self.assertEqual(len(colors), len(set(maped_categories.values())))

        for rep, category in maped_categories.items():
            if categories[rep] is None or isna(categories[rep]):
                self.assertEqual(category, NA_CATEGORY_NAME)
                self.assertEqual(colors[category], GREY_RGB)
            else:
                self.assertNotEqual(colors[category], GREY_RGB)


    def test_no_nas_int(self):
        categories = dict(enumerate([0, 1, 2, 3, 3, 4, 5, 5, 5, 7, 7]))
        self.do_no_nas_test(categories)


    def test_no_nas_bool(self):
        categories = dict(enumerate([True, False, True, False, True, False]))
        self.do_no_nas_test(categories)


    def test_no_nas_string(self):
        categories = dict(enumerate(['a', 'b', 'c', 'd', 'd', 'a', 'b', 'a']))
        self.do_no_nas_test(categories)


    def test_no_nas_float(self):
        categories = dict(enumerate([0.1, 1.1, 2.2, 3.3, 3.3, 4.4, 5.5, 5.6, 5.6, 7.7, 7.7]))
        self.do_no_nas_test(categories)


    def test_nas_int(self):
        categories = dict(enumerate([0, 1, 2, None, None, 4, 5, 5, 5, 7, 7]))
        self.do_na_test(categories)


    def test_nas_bool(self):
        categories = dict(enumerate([True, False, None, None, True, False]))
        self.do_na_test(categories)


    def test_nas_string(self):
        categories = dict(enumerate(['a', 'b', 'c', 'd', NA, NA, 'a', 'b', 'a']))
        self.do_na_test(categories)


    def test_nas_float(self):
        categories = dict(enumerate([None, None, 2.2, 3.3, 3.3, 4.4, 5.5, 5.6, 5.6, 7.7, 7.7]))
        self.do_na_test(categories)


    def test_invalid_format_fails(self):
        invalid_inputs = [{'sample.name': {0: 'RPMI-8226', 1: 'T47D', 2: 'H226', 3: '7-pool'}},
                          'hello', {'a': [], 'b': []}]

        for i in invalid_inputs:
            with self.assertRaises(AssertionError) as cm:
                _, _ = map_discrete_colors(i)
