import unittest

class TestFuzzyDictLookup(unittest.TestCase):
    def setUp(self):
        self.test_dict = {
            'glucose': 'C6H12O6',
            'fructose': 'C6H12O6',
            'sucrose': 'C12H22O11'
        }

    def test_exact_match(self):
        result = fuzzy_dict_lookup('glucose', self.test_dict)
        self.assertEqual(result, 'C6H12O6')

    def test_close_match_above_threshold(self):
        result = fuzzy_dict_lookup('glocose', self.test_dict)  # minor typo
        self.assertIsNone(result)


    def test_no_match_below_threshold(self):
        result = fuzzy_dict_lookup('xyz', self.test_dict, threshold=85)
        self.assertIsNone(result)

    def test_custom_threshold(self):
        result = fuzzy_dict_lookup('glocose', self.test_dict, threshold=80)
        self.assertEqual(result, 'C6H12O6')


if __name__ == '__main__':
    unittest.main()
