import unittest
import os
import json
import cobra
import gzip
import biopathopt

class TestSubunit(unittest.TestCase):
    def setUp(self):
        self.base_dir = os.path.dirname(os.path.realpath(__file__))
        self.uniprot_subunit_test = json.load(open(os.path.join(self.base_dir, 'test_files/e_coli_core_uniprot_subunit.json')))
        with gzip.open(os.path.join(self.base_dir, 'test_files/e_coli_core.json.gz'), 'rt', encoding='utf-8') as gz_file:
            model = cobra.io.load_json_model(gz_file)
        self.ecm = biopathopt.EnzymeConstrainedModel(model)
        
    def test_get_subunit_data(self):
        self.ecm._get_subunit_data()
        subunit_data = {g.id: g.annotation.get('number_of_subunits') for g in ecm.model.genes}
        self.assertDictEqual(subunit_data, self.uniprot_subunit_test)

    def test_get_subunit_num_1(self):
        #res = '1'
        protein_names = ['Acetate kinase', 'Acetokinase']
        gene_names = ['ackA']
        txt_subunit = 'Homodimer'
        res = ecm._get_subunit_num(
            protein_names, gene_names, txt_subunit
        )
        self.assertEqual(res, '1')
    
    def test_get_subunit_num_2(self):
        #res = '2'
        protein_names = ['Phosphoenolpyruvate-protein phosphotransferase', 'Phosphotransferase system, enzyme I']
        gene_names = ['ptsI']
        txt_subunit = 'Homodimer (PubMed:12705838, PubMed:17053069). Interacts with the pole-localizer protein TmaR (PubMed:33376208). Binding to TmaR is reversible as long as TmaR can get phosphorylated, whereas binding to non-phosphorylated TmaR is very strong and shifts the equilibrium toward binding (PubMed:33376208)'
        res = ecm._get_subunit_num(
            protein_names, gene_names, txt_subunit
        )
        self.assertEqual(res, '2')

    def test_get_subunit_num_4(self):
        #res = '4'
        protein_names = ['Fructose-1,6-bisphosphatase class 1', 'FBPase class 1', 'D-fructose-1,6-bisphosphate 1-phosphohydrolase class 1']
        gene_names = ['fbp']
        txt_subunit = 'Homotetramer. Phosphoenolpyruvate stabilizes the homotetramer'
        res = ecm._get_subunit_num(
            protein_names, gene_names, txt_subunit
        )
        self.assertEqual(res, '4')
