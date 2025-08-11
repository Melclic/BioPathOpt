import unittest
import os
import json

class TestExactPubChemSearch(unittest.TestCase):
    def setUp(self):
        self.data_obj = biopathopt.Data()
        self.base_dir = os.path.dirname(os.path.realpath(__file__))
        self.glucose_pubchem_str = json.load(open(os.path.join(self.base_dir, 'test_files/glucose_pubchem_str.json')))

    def test_single_result(self):
        glucose_dict = self.data_obj.exact_pubchem_search('glucose')
        self.assertDictEqual(self.glucose_pubchem_str, glucose_dict)

    def test_invalid_result(self):
        glucose_dict = self.data_obj.exact_pubchem_search('thisisnotacorrectcidname')
        self.assertDictEqual({}, glucose_dict)


if __name__ == '__main__':
    unittest.main()
