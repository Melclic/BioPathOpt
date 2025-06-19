import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)
import json
import logging
import difflib
import os
import tempfile
import time
import urllib.request
from tqdm import tqdm
import re
import gzip
import copy
import pickle

import compress_json
import networkx as nx
import numpy as np
import pandas as pd
import requests
import pubchempy as pcp

from typing import Tuple, Dict, List, Union, Optional, Any

from biopathopt.utils import fuzzy_dict_lookup

"""Collection of functions that fetch data from MetaNetX

TODO: use a data lake

Note:
mnxm --> MetaNetX molecule id
mnxr --> MetaNetX reaction id
biggm --> BIGG molecule id
biggr --> BIGG reaction id
"""

class Data:
    """Class to hold all the parsers and the methods to fetch public available
    data for reactions and metabolites etc...
    """

    def __init__(self):
        self.base_dir = os.path.dirname(os.path.realpath(__file__))
        self._chem_xref = None
        self._reac_xref = None
        self._chem_prop = None
        self._reac_prop = None
        self._biggm_mnxm = None
        self._biggr_mnxr = None
        self._inchikey_mnxm = None
        self._mnxm_inchikey = None
        self._inchikey2_mnxm = None
        self._g_depr_mnxm = None
        self._g_depr_mnxr = None
        self._keggr_mnxr = None
        self._keggm_mnxm = None
        self._chebim_mnxm = None
        self._molname_mnxm = None
        self._json_brenda = None
        self._brenda_ec_inchikey_kcat = None
        self._brenda_ec_inchikey_sa = None
        #  these are the user overwrites
        self.inchikey_overwrite = {
            "GPRLSGONYQIRFK-FTGQXOHASA-N": "GPRLSGONYQIRFK-UHFFFAOYSA-N",
            "OMHUCGDTACNQEX-YXHQYMSVSA-M": "OMHUCGDTACNQEX-OSHKXICASA-M",
            "QSIDJGUAAUSPMG-JYUJVNTISA-N": "QSIDJGUAAUSPMG-CULFPKEHSA-M",
            "ONVABDHFQKWOSV-ZJXASSBSSA-N": "ONVABDHFQKWOSV-HPUSYDDDSA-N",
        }
        self.mnxm_overwrite = {
                "MNXM38": "MNXM1105762",
                "MNXM01": "MNXM1",
        }
        self.biggr_overwrite = {"GLUDyi": "GLUDy"}
        #  other function specfic
        self.pubchem_search_cache = {}
        self.pubchem_min_start = 0.0
        self.pubchem_min_count = 0
        self.pubchem_sec_start = 0.0
        self.pubchem_sec_count = 0
        #self.ncbi = None

        """
        Includes the following plus others defined in the ECMpy package
                # get all the hardcoded cofactors that exist
                df_cofactors = pd.read_csv(os.path.join(self.base_dir, 'flatfiles', 'cofactor_inchi_201811.tsv'), sep='\t', comment='#', header=None)
                df_cofactors.columns = ['InChI', 'name', 'MetaNetX']
                self.metanetx_cofactors = [self.single_depr_mnxm(i) for y in df_cofactors.MetaNetX for i in str(y).split(',') if not pd.isna(y)]
        """
        self.mnxm_cofactors = ['MNXM1102419', 'MNXM3', 'MNXM1105762', 'MNXM731165', 'MNXM13', 'MNXM1103285', 'MNXM796', 'MNXM729302', 'MNXM90960', 'MNXM734750', 'MNXM1102191', 'MNXM257', 'MNXM02', 'MNXM1', 'MNXM728294', 'MNXM728062', 'MNXM51918', 'MNXM738702', 'MNXM1107902', 'MNXM360', 'MNXM1104559', 'MNXM739756', 'MNXM1107906', 'MNXM1105927', 'MNXM740692', 'MNXM10', 'MNXM230', 'MNXM232', 'MNXM732398', 'MNXM735437', 'MNXM95', 'MNXM8975', 'MNXM11', 'MNXM255', 'MNXM436', 'MNXM9', 'MNXM731949', 'MNXM178', 'MNXM1102152', 'MNXM124865', 'MNXM2229', 'MNXM128', 'MNXM1101474', 'MNXM2255', 'MNXM39', 'MNXM3654', 'MNXM726339', 'MNXM1108018', 'MNXM1103428', 'MNXM1102167', 'MNXM8978', 'MNXM35', 'MNXM729214', 'MNXM36', 'MNXM729215', 'MNXM924', 'MNXM286', 'MNXM162231', 'MNXM4133', 'MNXM490899', 'MNXM572', 'MNXM1104555', 'MNXM58', 'MNXM731166', 'MNXM741485', 'MNXM1231', 'MNXM732620', 'MNXM1103718', 'MNXM332', 'MNXM40333', 'MNXM727276', 'MNXM1104823', 'MNXM1101285', 'MNXM1103302', 'MNXM726712', 'MNXM730586', 'MNXM266', 'MNXM8', 'MNXM53428', 'MNXM137', 'MNXM5', 'MNXM653', 'MNXM537', 'MNXM452', 'MNXM740736', 'MNXM152', 'MNXM191', 'MNXM1104385', 'MNXM727888', 'MNXM2174', 'MNXM432', 'MNXM411', 'MNXM726711', 'MNXM107', 'MNXM726710', 'MNXM652', 'MNXM736415', 'MNXM1101868', 'MNXM169', 'MNXM27', 'MNXM735438', 'MNXM394', 'MNXM738430', 'MNXM1922', 'MNXM40414', 'MNXM4835', 'MNXM1102128', 'MNXM736654', 'MNXM1105936', 'MNXM1562', 'MNXM234', 'WATER', 'MNXM1107622', 'MNXM735978', 'MNXM4041', 'MNXM1102072', 'MNXM1103458', 'MNXM735047', 'MNXM727224', 'MNXM733186', 'MNXM1103553', 'MNXM344', 'MNXM1092965']

    # ############ 

    def construct_reaction_string(
            self,
            reactants: List[Tuple[str, str]], 
            products: List[Tuple[str, str]], 
            inchikey_levels: int = 3
    ) -> str:
        """Constructs a reaction string from reactants and products with truncated InChIKeys.

        This function generates a reaction string from lists of reactants and products. 
        Each reactant and product is represented by a coefficient and an InChIKey. 
        The InChIKey is truncated to the specified number of levels.

        Args:
            reactants (List[Tuple[str, str]]): A list of tuples where each tuple contains 
                the coefficient (str) and the InChIKey (str) for a reactant.
            products (List[Tuple[str, str]]): A list of tuples where each tuple contains 
                the coefficient (str) and the InChIKey (str) for a product.
            inchikey_levels (int, optional): The number of levels of the InChIKey to include 
                in the output string. Defaults to 3.

        Returns:
            str: The constructed reaction string in the format 'reactants = products'.

        Example:
            >>> reactants = [("2", "ABCDEF-ABCDEF-N"), ("1", "GHIJKL-GHIJKL-F")]
            >>> products = [("1", "MNOPQR-MNOPQR-N")]
            >>> construct_reaction_string(reactants, products, inchikey_levels=2)
            '2 ABCDEF-ABCDEF + 1 GHIJKL-GHIJKL = 1 MNOPQR-MNOPQR'
        """
        
        # Construct the reaction string for reactants
        reacts_str = ''
        for coefficient, inchikey in reactants:
            # Truncate the InChIKey to the specified number of levels
            truncated_inchikey = '-'.join(inchikey.split('-')[:inchikey_levels])
            reacts_str += f"{coefficient} {truncated_inchikey} + "
        reacts_str = reacts_str[:-3]  # Remove the trailing ' + '

        # Construct the reaction string for products
        prods_str = ''
        for coefficient, inchikey in products:
            # Truncate the InChIKey to the specified number of levels
            truncated_inchikey = '-'.join(inchikey.split('-')[:inchikey_levels])
            prods_str += f"{coefficient} {truncated_inchikey} + "
        prods_str = prods_str[:-3]  # Remove the trailing ' + '

        # Combine reactants and products into the final equation string
        return f"{reacts_str} = {prods_str}"

    def convert_mnxr_equation(
        self,
        reac_str: str,
        inchikey_levels: int = 3
    ) -> str:
        """Converts a MetaNetX reaction equation to a simplified InChIKey-based equation.

        This function parses a MetaNetX reaction equation, retrieves the InChIKeys
        for the reactants and products, and converts the equation to a simplified
        format based on the specified number of InChIKey levels.

        Args:
            reac_str (str): The MetaNetX reaction equation string.
            inchikey_levels (int, optional): The number of InChIKey levels to include
                in the output string. Defaults to 3.

        Returns:
            str: The simplified reaction equation string using truncated InChIKeys.
        """
        # Parse the reaction string into reactants and products
        reactants, products = self.parse_mnxr_equation(reac_str)

        # Convert reactants to their InChIKeys with the specified level of detail
        tmp_reactants = []
        for coefficient, mnxm in reactants:
            inchikey = self.mnxm_inchikey[mnxm]
            tmp_reactants.append([coefficient, inchikey])

        # Convert products to their InChIKeys with the specified level of detail
        tmp_products = []
        for coefficient, mnxm in products:
            inchikey = self.mnxm_inchikey[mnxm]
            tmp_products.append([coefficient, inchikey])

        return self.construct_reaction_string(tmp_reactants, tmp_products, inchikey_levels=inchikey_levels)

    def parse_mnxr_equation(
            self, 
            reac_str: str
    ) -> Tuple[List[List[str]], List[List[str]]]:
        """Parses a metabolic network equation string from MetaNetX into reactants and products.

        This function takes a metabolic network equation string and parses it into
        lists of reactants and products. The input string should be in the form
        'reactant1 + reactant2 = product1 + product2', where each component may
        contain quantities and compartment information separated by spaces and '@' symbols.

        Args:
            reac_str (str): The metabolic network equation string.

        Returns:
            Tuple[List[List[str]], List[List[str]]]: A tuple containing two lists:
                - The first list contains lists of parsed reactants.
                - The second list contains lists of parsed products.
        """
        # Split the equation string into reactants and products
        reactants_str, products_str = reac_str.split('=')

        # Split the reactants and products by the '+' symbol and strip any surrounding whitespace
        reactants = [r.strip() for r in reactants_str.split('+')]
        products = [p.strip() for p in products_str.split('+')]

        # Further split each reactant and product by spaces
        reactants = [r.split(' ') for r in reactants]
        products = [p.split(' ') for p in products]

        # Remove compartment information and other annotations following the '@' symbol
        reactants = [[item.split('@')[0] for item in r] for r in reactants]
        products = [[item.split('@')[0] for item in p] for p in products]

        return reactants, products

    # ############ MetaNetX specific #################

    def single_depr_mnxr(self, mnxr: str) -> str:
        """Check that the MetaNetX reaction ID is not deprecated and return the correct one.
        Because there may be multiple mnxr for a depr mnx and recursion
        because we can have a linked deprecated values

        Args:
            mnxr (str): MetaNetX reaction ID
        Returns:
            str: The non-deprecated MetaNetX reaction id
        """
        logging.debug("---- single_depr_mnxr ------")
        if mnxr:
            try:
                tmp = []
                for i in list(nx.dfs_tree(self.g_depr_mnxr, source=mnxr)):
                    if not list(self.g_depr_mnxr.successors(i)):
                        tmp.append(i)
                if len(tmp) == 1:
                    if tmp[0] != "EMPTY":
                        return tmp[0]
                elif len(tmp) > 1:
                    logging.debug(
                        "Cannot determine with certainty the ID "
                        + str(mnxr)
                        + ": "
                        + str(tmp)
                        + ". Returning the smallest value"
                    )
                    return sorted(tmp)[0]
            except (KeyError, nx.exception.NetworkXError) as e:
                pass
        return mnxr

    def single_depr_mnxm(self, mnxm, strict=False):
        """Because there may be multiple mnxm for a depr mnx. Recursion
        because we can have a linked deprecated values

        Args:
            mnxm (str): MetaNetX molecule ID
        Returns:
            str: The non-deprecated MetaNetX molecule id
        """
        #  logging.debug('---- single_depr_mnxm ------')
        if mnxm:
            try:
                return self.mnxm_overwrite[mnxm]
            except KeyError:
                pass
            try:
                tmp = []
                for i in list(nx.dfs_tree(self.g_depr_mnxm, source=mnxm)):
                    if not list(self.g_depr_mnxm.successors(i)):
                        tmp.append(i)
                if len(tmp) == 1:
                    return tmp[0]
                elif len(tmp) > 1:
                    if not strict:
                        logging.debug(
                            "Cannot determine with certainty the ID "
                            + str(mnxm)
                            + ": "
                            + str(tmp)
                            + ". Returning the smallest value"
                        )
                        return sorted(tmp)[0]
            except (KeyError, nx.exception.NetworkXError) as e:
                pass
        return mnxm

    def single_chem_prop(self, mnxm):
        """Unify the search of MetaNetX molecule id's. First assume that
        the id is correct, then check that its not an old ID and if not,
        check the deprecated maps and return the closest. Because of
        protonation you may have duplicates and thus the deprecated map.

        Args:
            mnxm (str): MetaNetX molecule ID
        Returns:
            str: The non-deprecated MetaNetX molecule id
        """
        #  logging.debug('------ single_chem_prop -------')
        try:
            return self.chem_prop[mnxm]
        except KeyError:
            try:
                return self.chem_prop[self.single_depr_mnxm(mnxm)]
            except KeyError:
                #  if all else fails, find the closest member of the deprecated map
                for i in list(
                    nx.bfs_tree(
                        self.g_depr_mnxm, self.single_depr_mnxm(mnxm), reverse=True
                    )
                ):
                    try:
                        return self.chem_prop[i]
                    except KeyError:
                        pass
        logging.warning("Cannot find chem properties for: " + str(mnxm))
        return {}

    def single_reac_prop(self, mnxr):
        """Unify the search of MetaNetX reaction id's. First assume that
        the id is correct, then check that its not an old ID. Because of
        protonation you may have duplicates and thus the deprecated map.

        Args:
            mnxm (str): MetaNetX reaction ID
        Returns:
            str: The non-deprecated MetaNetX reaction id
        """
        logging.debug("------ single_reac_prop ------")
        try:
            return self.reac_prop[mnxr]
        except KeyError:
            try:
                return self.reac_prop[self.single_depr_mnxr(mnxr)]
            except KeyError:
                pass
                """
                G = self.g_depr_mnxr()
                c_mnxr = single_depr_mnxr(mnxr)
                if c_mnxr in G.nodes:
                    for i in list(nx.bfs_tree(G, c_mnxr, reverse=True)):
                        try:
                            return rp[i]
                        except KeyError:
                            pass
                """
        logging.warning("Cannot find chem properties for: " + str(mnxr))
        return {}

    # ############# BRENDA ######################
    # Here we parse the brenda file and return the substrate kcat values
    
    @property
    def json_brenda(self):
        if not self._json_brenda:
            logging.debug("------ json_brenda -----")
            logging.debug("\t-> Populating...")
            path_file = os.path.join(
                self.base_dir, "flatfiles/json_brenda.json.gz"
            )
            if not os.path.exists(path_file):
                raise TypeError('Need to download the file at: https://brenda-enzymes.org/download.php')
            with gzip.open(path_file, 'rt', encoding='utf-8') as f:
                self._json_brenda = json.load(f)
        return self._json_brenda


    @property
    def brenda_ec_inchikey_kcat(self):
        if not self._brenda_ec_inchikey_kcat:
            logging.debug("------ brenda_ec_inchikey_kcat -----")
            logging.debug("\t-> Populating...")
            path_file = os.path.join(
                #self.base_dir, "flatfiles/brenda_ec_inchikey_kcat.pkl"
                self.base_dir, "flatfiles/brenda_kcat.pkl"
            )
            if os.path.exists(path_file):
                with open(path_file, 'rb') as f:
                    self._brenda_ec_inchikey_kcat = pickle.load(f)
            else:
                self._brenda_ec_inchikey_kcat = self._generate_brenda_kinetics(parse_type='kcat')
                with open(path_file, 'wb') as f:
                    pickle.dump(self._brenda_ec_inchikey_kcat, f, pickle.HIGHEST_PROTOCOL)
        return self._brenda_ec_inchikey_kcat


    @property
    def brenda_ec_inchikey_sa(self):
        if not self._brenda_ec_inchikey_sa:
            logging.debug("------ brenda_ec_inchikey_sa -----")
            logging.debug("\t-> Populating...")
            path_file = os.path.join(
                #self.base_dir, "flatfiles/brenda_ec_inchikey_sa.pkl"
                self.base_dir, "flatfiles/brenda_sa.pkl"
            )
            if os.path.exists(path_file):
                with open(path_file, 'rb') as f:
                    self._brenda_ec_inchikey_sa = pickle.load(f)
            else:
                self._brenda_ec_inchikey_sa = self._generate_brenda_kinetics(parse_type='sa')
                with open(path_file, 'wb') as f:
                    pickle.dump(self._brenda_ec_inchikey_sa, f, pickle.HIGHEST_PROTOCOL)
        return self._brenda_ec_inchikey_sa


    def _get_protein_reactants_inchikey(self, ec_number, protein_ids=[]):
        #because the comments do not always have the right substrates, get the original reaction from the protein

        def extract_reactants_inchikey(reactions_list, cofactor_list=[], filter_proteins_ids=[]):
            reactant_dict = {}
            for reac in reactions_list:
                reactants = [y.strip() for y in reac.get('value').split('=')[0].split('+')]
                for cof in cofactor_list:
                    if set(cof.get('proteins')) & set(reac.get('proteins')):
                        reactants = [i for i in reactants if i!=cof.get('value')]
                inchikeys = []
                for reactant in reactants:
                    inchikey = self.molecule_name_search_inchikey(reactant.lower())
                    if inchikey:
                        inchikeys.append(inchikey)
                if inchikeys:
                    for p in reac.get('proteins'):
                        if filter_proteins_ids:
                            if p not in filter_proteins_ids:
                                continue
                        if p not in reactant_dict:
                            reactant_dict[p] = []
                        if tuple(sorted(inchikeys)) not in reactant_dict[p]:
                            reactant_dict[p].append(tuple(sorted(inchikeys)))
                        reactant_dict[p] = list(set(reactant_dict[p]))
            return reactant_dict

        reactions_list = []
        if not reactions_list:
            if 'substrates_products' in self.json_brenda['data'][ec_number]:
                if protein_ids:
                    try:
                        reactions_list = [y for y in self.json_brenda['data'][ec_number]['substrates_products'] if set(protein_ids) & set(y['proteins'])]
                    except KeyError:
                        reactions_list = self.json_brenda['data'][ec_number]['substrates_products']
                else:
                    reactions_list = self.json_brenda['data'][ec_number]['substrates_products']
        if not reactions_list:
            if 'reaction' in self.json_brenda['data'][ec_number]:
                if protein_ids:
                    try:
                        reactions_list = [y for y in self.json_brenda['data'][ec_number]['reaction'] if set(protein_ids) & set(y['proteins'])]
                    except KeyError:
                        reactions_list = self.json_brenda['data'][ec_number]['reaction']
                else:
                    reactions_list = self.json_brenda['data'][ec_number]['reaction']
        if 'cofactor' in self.json_brenda['data'][ec_number]:
            if protein_ids:
                try:
                    cofactor_list = [y for y in self.json_brenda['data'][ec_number]['cofactor'] if set(protein_ids) & set(y['proteins'])]
                except KeyError:
                    cofactor_list = self.json_brenda['data'][ec_number]['cofactor']
            else:
                cofactor_list = self.json_brenda['data'][ec_number]['cofactor']
        else:
            cofactor_list = []
        reactant_dict = extract_reactants_inchikey(
            reactions_list=reactions_list, 
            cofactor_list=cofactor_list,
            filter_proteins_ids=protein_ids,
        )
        return reactant_dict


    def _extract_brenda_substrates_sa(self, comment):
        def is_single_word(text: str) -> bool:
            return bool(re.fullmatch(r'\S+', text.strip()))
        
        def extract_substrate(text: str) -> str:
            if is_single_word(text):
                return text
            match = re.search(r'substrate:\s*([A-Za-z0-9\-]+)', text, re.IGNORECASE)
            if match:
                return match.group(1)
            match = re.search(r'substrate\s*([A-Za-z0-9\-]+)', text, re.IGNORECASE)
            if match:
                return match.group(1)
            match = re.search(r'([A-Za-z0-9\-]+) as substrate', text, re.IGNORECASE)
            if match:
                return match.group(1)
            return None
        if 'substrate' in comment.lower():
            ### Try to extract the substrate
            result = extract_substrate(comment)
            if result:
                result = result.replace('.', '-').lower()
            return result
        return None

    def _extract_brenda_substrates_kcat(self, entry):
        match = re.search(r'([+-]?\d*\.?\d+)\s*\{([^}]+)\}', entry)
        if match:
            value = float(match.group(1))
            substrate = match.group(2).replace('.', '-').lower()
            return value, substrate
        return None, None

    def _generate_brenda_kinetics(self, parse_type='kcat'):
        if parse_type not in ['kcat', 'sa']:
            raise TypeError(f'Input Must be kcat or sa: {parse_type}')
        res = {}
        for ec_number in tqdm(self.json_brenda['data'], desc=f"Generate Brenda {parse_type} File"):
            if parse_type=='sa':
                try:
                    kinetics_list = copy.deepcopy(self.json_brenda['data'][ec_number]['specific_activity'])
                except KeyError:
                    continue
                try:
                    protein_dict = copy.deepcopy(self.json_brenda['data'][ec_number]['protein'])
                except KeyError:
                    continue
            elif parse_type=='kcat':
                try:
                    kinetics_list = copy.deepcopy(self.json_brenda['data'][ec_number]['kcat_km_value'])
                except KeyError:
                    continue
                try:
                    protein_dict = copy.deepcopy(self.json_brenda['data'][ec_number]['protein'])
                except KeyError:
                    continue
            ### get the reactants from the reaction description
            try:
                default_reactants = self._get_protein_reactants_inchikey(
                    ec_number=ec_number 
                )
            except KeyError:
                default_reactants = {}
            res[ec_number] = {}
            for kinetic_entry in kinetics_list:
                kinetic_entry['proteins'] = [protein_dict[y] for y in kinetic_entry.get('proteins')]
                if 'mutant' in kinetic_entry.get('comment').lower() or 'mutated' in kinetic_entry.get('comment').lower():
                    continue
                to_add_org = [y.get('organism') for y in kinetic_entry.get('proteins')]
                to_add_org = [y for y in to_add_org if y]       
                sub_id = [('no_identifiable_substrate',)]
                #extract from the comment
                if parse_type=='sa':
                    extracted_sub = self._extract_brenda_substrates_sa(kinetic_entry.get('comment').lower())
                    value = kinetic_entry.get('value')
                    if value and value>0.0:
                        value = float(value)
                    else:
                        continue
                elif parse_type=='kcat':
                    value, extracted_sub = self._extract_brenda_substrates_kcat(kinetic_entry.get('value').lower())
                    if not value and value>0.0:
                        continue
                if extracted_sub:
                    search_inchikey = self.molecule_name_search_inchikey(extracted_sub.lower())
                    if search_inchikey:
                        sub_id = [search_inchikey]
                #because the comments do not always have the right substrates, get the original reaction from the protein
                if sub_id==[('no_identifiable_substrate',)]:
                    if default_reactants:
                        tmp_s = []
                        prots = kinetic_entry.get('proteins')
                        if prots:
                            try:
                                for y in [default_reactants[i.get('id')] for i in prots]:
                                    #make sure there are no None
                                    tmp_s += [x for x in y if x]
                            except KeyError:
                                pass
                        if tmp_s:
                            sub_id = tmp_s
                #### if uniprit is input check
                try:
                    to_add_uniprot = [y.get('accessions') for y in kinetic_entry.get('proteins') if y.get('source')=='uniprot']
                    tmp_u = []
                    for i in to_add_uniprot:
                        if i:
                            if isinstance(i, list):
                                for y in i:
                                    if y:
                                        tmp_u.append(y)
                            elif isinstance(i, str):
                                tmp_u.append(i)
                    to_add_uniprot = tmp_u
                except (TypeError, KeyError) as e:
                    to_add_uniprot = []
                sub_id = tuple(sub_id)
                if not sub_id in res[ec_number]:
                    res[ec_number][sub_id] = {}
                for org in to_add_org:
                    if org not in res[ec_number][sub_id]:
                        res[ec_number][sub_id][org] = {'uniprot': [], 'values': [], 'comments': []}
                    res[ec_number][sub_id][org]['values'].append(value)
                    res[ec_number][sub_id][org]['uniprot'] += to_add_uniprot
                    res[ec_number][sub_id][org]['comments'].append(kinetic_entry.get('comment').lower())
        return res


    ##### taxonomy ###

    '''
    def _get_species_name(self, taxid: int) -> str:
        """Return the species-level name from a taxonomy ID."""
        if not self.ncbi:
            self.ncbi = NCBITaxa()
        lineage = self.ncbi.get_lineage(taxid)
        names = self.ncbi.get_taxid_translator(lineage)
        ranks = self.ncbi.get_rank(lineage)

        for tid in lineage:
            if ranks[tid] == "species":
                return names[tid]
        return names.get(taxid, None)


    def _get_taxid_from_species(self, species_name: str) -> Optional[int]:
        """
        Retrieve the taxonomy ID (taxid) for a given species name using ete3.

        Args:
            species_name (str): The scientific name of the species (e.g., 'Escherichia coli').

        Returns:
            Optional[int]: The taxonomy ID if found, otherwise None.

        Example:
            >>> get_taxid_from_species("Escherichia coli")
            562
        """
        if not self.ncbi:
            self.ncbi = NCBITaxa()
        try:
            name2taxid = self.ncbi.get_name_translator([species_name])
            return name2taxid[species_name][0] if species_name in name2taxid else None
        except Exception as e:
            logging.warning(f"Failed to retrieve taxid for '{species_name}': {e}")
            return None
    '''


    # ######### data properties ###################################
    #  These are designed to behave as parameters. The first time its called it
    #  will load, and the next time around it will pass the saved parameter

    #  ######### MetaNetX ###################
    #  taken from metnetx website. All the deprecated version of the ID's to the current id's
    #  important because we use an older version of metanetx id's

    @property
    def g_depr_mnxr(self):
        """The metanetx and the deprecated metanetx id's

        The function will download https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_depr.tsv
        MetaNetX file and reorganize it into dictionnaries

        Args:
        Returns:
            nx.Graph: Get the reaction deprecated network maps
        """
        if not self._g_depr_mnxr:
            logging.debug("------ g_depr_mnxr -----")
            logging.debug("\t-> Populating...")
            path_g_depr_mnxr = os.path.join(
                self.base_dir, "flatfiles/g_depr_mnxr.json.gz"
            )
            if not os.path.exists(path_g_depr_mnxr):
                with tempfile.TemporaryDirectory() as tmpdirname:
                    ori_reac_depr_path = os.path.join(tmpdirname, "reac_depr.tsv")
                    urllib.request.urlretrieve(
                        "https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_depr.tsv",
                        ori_reac_depr_path,
                    )
                    reac_depr = pd.read_csv(
                        ori_reac_depr_path,
                        comment="#",
                        sep="\t",
                        header=None,
                    )
                    reac_depr.columns = ["deprecated_ID", "ID", "version"]
                    self._g_depr_mnxr = nx.DiGraph()
                    for n in np.unique(
                        reac_depr["deprecated_ID"].to_list() + reac_depr["ID"].to_list()
                    ):
                        self._g_depr_mnxr.add_node(n)
                    for index, row in reac_depr.iterrows():
                        self._g_depr_mnxr.add_edge(row["deprecated_ID"], row["ID"])
                    compress_json.dump(
                        nx.cytoscape_data(self._g_depr_mnxr), path_g_depr_mnxr
                    )
            self._g_depr_mnxr = nx.cytoscape_graph(compress_json.load(path_g_depr_mnxr))
        return self._g_depr_mnxr

    @property
    def chem_xref(self):
        """The chem xref describing the different database of the cross reference of a chemical
        species

        Downloading the follwing file https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_xref.tsv from
        MetaNetX describing the cross reference of chemical species. The function will download the file
        if it does not exists. Then it will make the data frame

        Args:
        Returns:
            pd.DataFrame: the cross reference for MetaNetX
        """
        #  taken from the metanetx website. all chemical cross links
        if not isinstance(self._chem_xref, pd.DataFrame):
            logging.debug("------ chem_xref -----")
            logging.debug("\t-> Populating...")
            chem_xref_path = os.path.join(self.base_dir, "flatfiles/chem_xref.tsv.xz")
            if not os.path.exists(chem_xref_path):
                logging.info(
                    "The MetaNetX file chem_xref.tsv does not exist... downloading it"
                )
                with tempfile.TemporaryDirectory() as tmpdirname:
                    ori_chem_xref_path = os.path.join(tmpdirname, "chem_xref.tsv")
                    urllib.request.urlretrieve(
                        "https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_xref.tsv",
                        ori_chem_xref_path,
                    )
                    chem_xref = pd.read_csv(
                        ori_chem_xref_path,
                        comment="#",
                        sep="\t",
                        header=None,
                    )
                    chem_xref.to_csv(chem_xref_path)
            chem_xref = pd.read_csv(
                chem_xref_path,
                index_col=0,
            )
            chem_xref.columns = ["source", "ID", "description"]
            tmp_db = []
            tmp_entry = []
            for i in chem_xref["source"]:
                try:
                    tmp_db.append(i.split(":")[0])
                except (KeyError, IndexError):
                    tmp_db.append(np.nan)
                try:
                    tmp_entry.append(i.split(":")[1])
                except (KeyError, IndexError):
                    tmp_entry.append(i)
            chem_xref.index = chem_xref["source"]
            chem_xref["source_db"] = tmp_db
            chem_xref["source_entry"] = tmp_entry
            self._chem_xref = chem_xref
        return self._chem_xref

    @property
    def reac_xref(self):
        """The chem xref describing the different database of the cross reference of a chemical
        species

        Downloading the follwing file https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_xref.tsv from
        MetaNetX describing the cross reference of chemical species. The function will download the file
        if it does not exists. Then it will make the data frame

 
        Args:
        Returns:
            pd.DataFrame: The reaction cross-reference dataframe
        """
        #  taken from the metanetx website. all chemical cross links
        if not isinstance(self._reac_xref, pd.DataFrame):
            logging.debug("------ reac_xref -----")
            logging.debug("\t-> Populating...")
            reac_xref_path = os.path.join(self.base_dir, "flatfiles/reac_xref.tsv.xz")
            if not os.path.exists(reac_xref_path):
                logging.info(
                    "The MetaNetX file reac_xref.tsv.xz does not exist... downloading it"
                )
                with tempfile.TemporaryDirectory() as tmpdirname:
                    ori_reac_xref_path = os.path.join(tmpdirname, "reac_xref.tsv")
                    urllib.request.urlretrieve(
                        "https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_xref.tsv",
                        ori_reac_xref_path,
                    )
                    reac_xref = pd.read_csv(
                        ori_reac_xref_path,
                        comment="#",
                        sep="\t",
                        header=None,
                    )
                    reac_xref.to_csv(reac_xref_path)
            reac_xref = pd.read_csv(
                reac_xref_path,
                index_col=0,
            )
            reac_xref.columns = ["source", "ID", "description"]
            tmp_db = []
            tmp_entry = []
            for i in reac_xref["source"]:
                try:
                    tmp_db.append(i.split(":")[0])
                except (KeyError, IndexError):
                    tmp_db.append(np.nan)
                try:
                    tmp_entry.append(i.split(":")[1])
                except (KeyError, IndexError):
                    tmp_entry.append(i)
            reac_xref.index = reac_xref["source"]
            reac_xref["source_db"] = tmp_db
            reac_xref["source_entry"] = tmp_entry
            self._reac_xref = reac_xref
        return self._reac_xref

    #  taken from metnetx website. All the deprecated version of the ID's to the current id's
    #  important because we use an older version of metanetx id's
    @property
    def g_depr_mnxm(self):
        """The metanetx and the deprecated metanetx id's

        The function will download https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_depr.tsv
        MetaNetX file and reorganize it into dictionnaries

        Args:
        Returns:
            dict: Return dictionnary of Metanetx deprecated id's
        """
        if not self._g_depr_mnxm:
            logging.debug('------ g_depr_mnxm -----')
            logging.debug("\t-> Populating...")
            path_g_depr_mnxm = os.path.join(
                self.base_dir, "flatfiles/g_depr_mnxm.json.gz"
            )
            if not os.path.exists(path_g_depr_mnxm):
                with tempfile.TemporaryDirectory() as tmpdirname:
                    logging.info(
                        "The MetaNetX file mnxm_deprmnxm.tsv does not exist... downloading it"
                    )
                    ori_chem_depr_path = os.path.join(tmpdirname, "chem_depr.tsv")
                    urllib.request.urlretrieve(
                        "https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_depr.tsv",
                        ori_chem_depr_path,
                    )
                    chem_depr = pd.read_csv(
                        ori_chem_depr_path,
                        comment="#",
                        sep="\t",
                        header=None,
                    )
                    chem_depr.columns = ["deprecated_ID", "ID", "version"]
                    G = nx.DiGraph()
                    for n in np.unique(
                        chem_depr["deprecated_ID"].to_list() + chem_depr["ID"].to_list()
                    ):
                        G.add_node(n)
                    for index, row in chem_depr.iterrows():
                        G.add_edge(row["deprecated_ID"], row["ID"])
                    compress_json.dump(nx.cytoscape_data(G), path_g_depr_mnxm)
            self._g_depr_mnxm = nx.cytoscape_graph(compress_json.load(path_g_depr_mnxm))
        return self._g_depr_mnxm

    @property
    def chem_prop(self):
        """Return the chemical properties of molecules

        This function will download the following file https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_prop.tsv that
        describes the chemical structure, etc...

        Args:
        Returns:
            dict: Chemical properties from MetaNetX
        """
        if not self._chem_prop:
            logging.debug("------ chem_prop ----")
            logging.debug("\t-> Populating...")
            chem_prop_path = os.path.join(self.base_dir, "flatfiles/chem_prop.json.gz")
            if not os.path.exists(chem_prop_path):
                logging.info("Chem prop does not exist... generating it")
                with tempfile.TemporaryDirectory() as tmpdirname:
                    chem_prop_ori_path = os.path.join(tmpdirname, "chem_prop.tsv")
                    urllib.request.urlretrieve(
                        "https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_prop.tsv",
                        chem_prop_ori_path,
                    )
                    chem_prop = pd.read_csv(
                        chem_prop_ori_path, comment="#", sep="\t", header=None
                    )
                    #  Need to do this since id's are now mutiple and contain conflicts
                    #  chem_prop["ID"] = [single_depr_mnxm(i) for i in chem_prop["ID"]]
                    chem_prop.columns = [
                        "ID",
                        "name",
                        "reference",
                        "formula",
                        "charge",
                        "mass",
                        "InChI",
                        "InChIKey",
                        "SMILES",
                    ]
                    chem_prop["InChIKey"] = [
                        i.replace("InChIKey=", "") if not pd.isna(i) else i
                        for i in chem_prop["InChIKey"]
                    ]
                    chem_prop = chem_prop.set_index("ID")
                    chem_prop = chem_prop.transpose().to_dict()
                    compress_json.dump(chem_prop, chem_prop_path)
                    chem_prop = None
            self._chem_prop = compress_json.load(chem_prop_path)
        return self._chem_prop

    @property
    def reac_prop(self):
        """Return the chemical properties of molecules

        This function will download the following file https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_prop.tsv that
        describes the chemical structure, etc...

        Args:
        Returns:
            dict: Reaction properties from MetaNetX
        """
        if not self._reac_prop:
            logging.debug("------ reac_prop -----")
            logging.debug("\t-> Populating...")
            reac_prop_path = os.path.join(self.base_dir, "flatfiles/reac_prop.json.gz")
            if not os.path.exists(reac_prop_path):
                with tempfile.TemporaryDirectory() as tmpdirname:
                    reac_prop_ori_path = os.path.join(tmpdirname, "reac_prop.tsv")
                    if not os.path.exists(reac_prop_ori_path):
                        logging.info(
                            "The MetaNetX file reac_prop.tsv does not exist... downloading it"
                        )
                        urllib.request.urlretrieve(
                            "https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_prop.tsv",
                            reac_prop_ori_path,
                        )
                    reac_prop = pd.read_csv(
                        reac_prop_ori_path,
                        comment="#",
                        sep="\t",
                        header=None,
                    )
                    reac_prop.columns = [
                        "ID",
                        "mnx_equation",
                        "reference",
                        "classifs",
                        "is_balanced",
                        "is_transport",
                    ]
                    reac_prop = reac_prop.set_index("ID")
                    reac_prop = reac_prop.transpose().to_dict()
                    for i in reac_prop:
                        try:
                            reac_prop[i]["classifs"] = reac_prop[i]["classifs"].split(
                                ";"
                            )
                        except AttributeError:
                            reac_prop[i]["classifs"] = ""
                    ### generate the reaction based in inchikeys
                    for i in reac_prop:
                        try:
                            reac_prop[i]['inchikey2_equation'] = \
                                self.convert_mnxr_equation(reac_prop[i]['mnx_equation'], inchikey_levels=2)
                        except (ValueError, KeyError) as e:
                            reac_prop[i]['inchikey2_equation'] = ''
                        try:
                            reac_prop[i]['inchikey_equation'] = \
                                self.convert_mnxr_equation(reac_prop[i]['mnx_equation'], inchikey_levels=3)
                        except (ValueError, KeyError) as e:
                            reac_prop[i]['inchikey_equation'] = ''
                    #save it
                    compress_json.dump(reac_prop, reac_prop_path)
            self._reac_prop = compress_json.load(reac_prop_path)
        return self._reac_prop

    @property
    def biggm_mnxm(self):
        """Return the data.frame with bigg id

        The two have a 1:1 relationship

        Args:
        Returns:
            dict: BIGG molecule to MetNetX molecule id
        """
        if not self._biggm_mnxm:
            logging.debug("------ biggm_mnxm -----")
            logging.debug("\t-> Populating...")
            tmp = self.chem_xref[self.chem_xref["source_db"] == "bigg.metabolite"]
            tmp = tmp[["ID", "source_entry"]]
            # tmp['source_entry'] = [i.lower() for i in tmp['source_entry']]
            tmp = tmp.groupby("source_entry")["ID"].apply(list).to_dict()
            for i in tmp:
                if not len(tmp[i]) == 1:
                    raise KeyError(
                        "The 1:1 assumption for biggm_mnxm is not respected for "
                        + str(i)
                    )
            self._biggm_mnxm = {i: tmp[i][0] for i in tmp}
        return self._biggm_mnxm

    @property
    def biggr_mnxr(self):
        """Return the data.frame with bigg id

        Assumption that there is a 1:1 relationship between the two

        Args:
        Returns:
            dict: BIGG reaction ID to MetNetX reaction ID
        """
        if not self._biggr_mnxr:
            logging.debug('------ biggr_mnxr ----')
            logging.debug("\t-> Populating...")
            tmp = self.reac_xref[self.reac_xref["source_db"] == "bigg.reaction"]
            tmp = tmp[["ID", "source_entry"]]
            # tmp['source_entry'] = [i.lower() for i in tmp['source_entry']]
            tmp = tmp.groupby("source_entry")["ID"].apply(list).to_dict()
            for i in tmp:
                if not len(tmp[i]) == 1:
                    raise KeyError(
                        "The 1:1 assumption for biggr_mnxr not respected for "
                        + str(i)
                    )
            self._biggr_mnxr = {i: tmp[i][0] for i in tmp}
        return self._biggr_mnxr

    @property
    def inchikey2_mnxm(self):
        """Based on the inchikey return the mnxm id

        Args:
        Returns:
            dict: InChIkey ID to MetNetX molecule id
        """
        if not self._inchikey2_mnxm:
            logging.debug("------ inchikey_mnxm -----")
            logging.debug("\t-> Populating...")
            self._inchikey2_mnxm = {}
            for i in self.chem_prop:
                try:
                    if not pd.isna(self.chem_prop[i]["InChIKey"]) or not self.chem_prop[i]["InChIKey"]=='nan':
                        self._inchikey2_mnxm[
                            '-'.join(str( self.chem_prop[i]['InChIKey'] ).split('-')[:2])
                        ] = self.single_depr_mnxm(i)
                except KeyError:
                    pass
        return self._inchikey2_mnxm


    @property
    def inchikey_mnxm(self):
        """Based on the inchikey return the mnxm id

        Args:
        Returns:
            dict: InChIkey ID to MetNetX molecule id
        """
        if not self._inchikey_mnxm:
            logging.debug("------ inchikey_mnxm -----")
            logging.debug("\t-> Populating...")
            self._inchikey_mnxm = {}
            for i in self.chem_prop:
                try:
                    if not pd.isna(self.chem_prop[i]["InChIKey"]):
                        self._inchikey_mnxm[
                            self.chem_prop[i]["InChIKey"]
                        ] = self.single_depr_mnxm(i)
                except KeyError:
                    pass
        return self._inchikey_mnxm

    @property
    def mnxm_inchikey(self):
        """Based on the mxnm return the inchikey

        Args:
        Returns:
            dict: the mnxm to inchikey
        """
        if not self._mnxm_inchikey:
            logging.debug("------ mnxm_inchikey -----")
            logging.debug("\t-> Populating...")
            self._mnxm_inchikey = {}
            for i in self.chem_prop:
                try:
                    if not pd.isna(self.chem_prop[i]["InChIKey"]):
                        self._mnxm_inchikey[
                                self.single_depr_mnxm(i)
                        ] = self.chem_prop[i]["InChIKey"]
                except KeyError:
                    pass
        return self._mnxm_inchikey

    @property
    def keggr_mnxr(self):
        """KEGG reaction ID to MetaNetX reaction ID

        Assuming a 1:1 relationship

        Args:
        Returns:
            dict: KEGG reaction ID to MetaNetX reaction ID
        """
        if not self._keggr_mnxr:
            logging.debug("------ keggr_mnxr ------")
            logging.debug("\t-> Populating...")
            tmp = self.reac_xref[self.reac_xref["source_db"] == "kegg.reaction"]
            tmp = tmp[["source_entry", "ID"]]
            tmp = tmp[[not pd.isna(i) for i in tmp["ID"]]]
            tmp = tmp.groupby("source_entry")["ID"].apply(list).to_dict()
            for i in tmp:
                if not len(tmp[i]) == 1:
                    raise KeyError(
                        "The 1:1 assumption for keggr_mnxr is not respected for "
                        + str(i)
                    )
            self._keggr_mnxr = {i: tmp[i][0] for i in tmp}
        return self._keggr_mnxr

    @property
    def keggm_mnxm(self):
        """KEGG molecule ID to MetaNetX molecule ID

        Assuming a 1:1 relationship

        Args:
        Returns:
            dict: KEGG molecule ID to MetaNetX molecule ID
        """
        if not self._keggm_mnxm:
            logging.debug("------ keggm_mnxm ------")
            logging.debug("\t-> Populating...")
            tmp = self.chem_xref[self.chem_xref["source_db"] == "kegg.compound"]
            tmp = tmp[["source_entry", "ID"]]
            tmp = tmp[[not pd.isna(i) for i in tmp["ID"]]]
            tmp = tmp.groupby("source_entry")["ID"].apply(list).to_dict()
            for i in tmp:
                if not len(tmp[i]) == 1:
                    raise KeyError(
                        "The 1:1 assumption for keggm_mnxm is not respected for "
                        + str(i)
                    )
            self._keggm_mnxm = {i: tmp[i][0] for i in tmp}
        return self._keggm_mnxm

    @property
    def chebim_mnxm(self):
        """KEGG molecule ID to MetaNetX molecule ID

        Assuming a 1:1 relationship

        Args:
        Returns:
            dict: KEGG molecule ID to MetaNetX molecule ID
        """
        if not self._chebim_mnxm:
            logging.debug("------ keggm_mnxm ------")
            logging.debug("\t-> Populating...")
            tmp1 = self.chem_xref[self.chem_xref["source_db"]=="CHEBI"]
            tmp2 = self.chem_xref[self.chem_xref["source_db"]=="chebi"]
            tmp = pd.concat([tmp1, tmp2])
            tmp = tmp[["source_entry", "ID"]]
            tmp = tmp[[not pd.isna(i) for i in tmp["ID"]]]
            tmp = tmp.groupby("source_entry")["ID"].apply(list).to_dict()
            tmp = {i: list(np.unique(tmp[i])) for i in tmp}
            for i in tmp:
                if not len(tmp[i]) == 1:
                    raise KeyError(
                        "the 1:1 assumption for chebim_mnxm is not respected for "
                        + str(i)
                    )
            self._chebim_mnxm = {i: tmp[i][0] for i in tmp}
        return self._chebim_mnxm

    @property
    def molname_mnxm(self):
        """molecule name to mnxm
        """
        if not self._molname_mnxm:
            logging.debug("------ molname_mnxm ------")
            logging.debug("\t-> Populating...")
            self._molname_mnxm = {}
            for mnxm in self.chem_prop:
                if 'name' in self.chem_prop[mnxm]:
                    if not pd.isna(self.chem_prop[mnxm]['name']):
                        if self.chem_prop[mnxm]['name'] not in self._molname_mnxm:
                            self._molname_mnxm[self.chem_prop[mnxm]['name']] = mnxm
                        #else:
                        #    logging.warning(f"Duplicate name entries for {self.chem_prop[mnxm]['name']}")
        return self._molname_mnxm

    # ################# search functions ###############

    # #### Pubchem #####

    def exact_pubchem_search(self, query: str, itype: str = 'name') -> Dict[str, Any]:
        """
        Perform an exact search on PubChem using the given identifier.

        Args:
            query (str): The compound name or identifier to search.
            itype (str): The type of identifier (e.g., 'name', 'smiles', 'inchi'). Defaults to 'name'.

        Returns:
            Dict[str, Any]: A dictionary containing the compound data if found. Empty dict if not found or ambiguous.

        Raises:
            KeyError: If multiple compounds are returned for the query.
        """
        logging.debug(f'Query: {query} - itype: {itype}')
        
        if not query:
            return {}

        # Check if result is already cached
        try:
            return self.pubchem_search_cache[query.lower()]
        except KeyError:
            pass

        try:
            # Perform the PubChem search
            cids = pcp.get_compounds(identifier=query, namespace=itype)
        except ValueError:
            # Cache the empty result and return
            self.pubchem_search_cache[query.lower()] = {}
            return {}

        if len(cids) == 0:
            return {}
        elif len(cids) == 1:
            # Cache and return the single result
            cid = cids[0].to_dict()
            self.pubchem_search_cache[query.lower()] = cid
            return cid
        else:
            # Cache the empty result and raise an error for ambiguity
            self.pubchem_search_cache[query.lower()] = {}
            raise KeyError(f'Multiple cids {cids} for {query}')
        return {}
    

    def molecule_name_search_inchikey(self, name: str) -> Optional[str]:
        """
        Search for the InChIKey of a molecule based on its name using cached data,
        fuzzy matching with MNX identifiers, and PubChem search.

        Args:
            name (str): The name of the molecule to search.

        Returns:
            Optional[str]: The InChIKey of the molecule if found, otherwise None.
        """
        try:
            # Attempt lookup in cached PubChem search results
            try:
                cid = self.pubchem_search_cache[name.lower()]
                inchikey = cid.get('inchikey')
                if inchikey:
                    return inchikey
            except KeyError:
                pass
            k


            # Fuzzy match against internal MNX mapping
            fuzzy_mnxm = fuzzy_dict_lookup(name.lower(), self.molname_mnxm)
            if fuzzy_mnxm:
                try:
                    return self.mnxm_inchikey[fuzzy_mnxm]
                except KeyError:
                    pass

            # Fallback to exact PubChem search by name
            try:
                cid = self.exact_pubchem_search(name.lower(), 'name')
                inchikey = cid.get('inchikey')
                if inchikey:
                    return inchikey
            except (KeyError, AttributeError):
                pass

        except Exception:
            pass

        return None


    # ########### find cross-references ###################

    def keggr_xref(self, kegg_id):
        """Return the cross-reference of the KEGG ID of reactions"""
        _kegg_id = kegg_id.upper()
        if not _kegg_id[0] == "R":
            _kegg_id = "R" + str(_kegg_id)
        try:
            _mnxr = self.keggr_mnxr[_kegg_id]
        except KeyError:
            logging.error("Cannot find the following KEGG id: " + str(_kegg_id))
            return {}, None
        xref = self.mnxr_xref(_mnxr)
        if not xref:
            _mnxr = self.single_depr_mnxr(_mnxr, strict=True)
            xref = self.mnxr_xref(_mnxr)
        return self.mnxr_xref(_mnxr)

    def keggm_xref(self, kegg_id):
        """Return the cross-reference of the KEGG ID of molecule"""
        # _kegg_id = kegg_id.lower()
        _kegg_id = kegg_id.upper()
        if not _kegg_id[0] == "C":
            _kegg_id = "C" + str(_kegg_id)
        try:
            _mnxm = self.keggm_mnxm[_kegg_id]
        except KeyError:
            logging.error("Cannot find the following KEGG id: " + str(_kegg_id))
            return {}, None
        return self.mnxm_xref(_mnxm)

    def biggm_xref(self, bigg_id):
        """Find the cross references associated with a molecule BIGG ID

        Args:
            bigg_id (str): The BIGG ID of a molecule
        Returns:
            tuple: A dictionnary of cross references and the MetaNetX ID
        """
        try:
            mnxm = self.biggm_mnxm[bigg_id]
            return self.mnxm_xref(mnxm)
        except KeyError:
            pass
        logging.warning("Cannot find the following BIGG molecule: " + str(bigg_id))
        return {}, None

    def biggr_xref(self, bigg_id):
        """Find the cross references associated with a reaction BIGG ID

        Args:
            bigg_id (str): The BIGG ID of a reaction
        Returns:
            tuple: A dictionnary of cross references and the MetaNetX ID
        """
        #  TODO: this is a small fix but not great
        _bigg_id = bigg_id.replace("_copy1", "")
        _bigg_id = _bigg_id.replace("_copy2", "")
        try:
            _bigg_id = self.biggr_overwrite[_bigg_id]
        except KeyError:
            pass
        try:
            mnxr = self.biggr_mnxr[_bigg_id]
            return self.mnxr_xref(mnxr)
        except KeyError:
            pass
        logging.warning("Cannot find the following BIGG molecule: " + str(bigg_id))
        return {}, None

    def inchikey_xref(self, inchikey, ignore_stereo=True):
        """Find all the cross references associated with a given inchikey

        Args:
            inchikey (str): The inchikey structure description of a molecule.
                Must be a string containing two - and a single letter at the end
            ignore_stereo (bool): Include or ignore the last dimension of the structure
        Returns:
            tuple: A dictionnary of cross references and the MetaNetX ID
        """
        #  test to see if we should overwrite a give inchikey
        try:
            inchikey = self.inchikey_overwrite[inchikey]
        except KeyError:
            pass
        try:
            mnxm = self.inchikey_mnxm[inchikey]
            return self.mnxm_xref(mnxm)
        except KeyError:
            if ignore_stereo:
                tmp = [
                    i
                    for i in self.inchikey_mnxm
                    if "-".join(inchikey.split("-")[:2]) in i
                ]
                if len(tmp) == 1:
                    mnxm = self.inchikey_mnxm[tmp[0]]
                    return self.mnxm_xref(mnxm)
                else:
                    logging.debug(
                        "Cannot find the unique inchi when ignoring the stereo: "
                        + str(tmp)
                    )
                    pass
        logging.warning("Cannot find the follwing inchikey: " + str(inchikey))
        return {}, None

    def mnxm_xref(self, mnxm):
        """Return the cross reference for a MetaNetX molecule ID

        Args:
            mnxm (str): MetaNetX species id
        Returns:
            Dictionnary with the cross reference
        """
        _mnxm = self.single_depr_mnxm(mnxm.upper())
        tmp = self.chem_xref[self.chem_xref["ID"] == _mnxm]["source"].to_dict()
        xref = {tmp[i].split(":")[0].lower(): [] for i in tmp if "MNX" not in tmp[i]}
        for i in tmp:
            if "MNX" not in tmp[i]:
                if ":" in tmp[i]:
                    xref[tmp[i].split(":")[0].lower()].append(
                        tmp[i].split(":")[1].lower()
                    )
        #  properties
        cp = self.single_chem_prop(_mnxm)
        cp["xref"] = xref
        if "metanetx.chemical" not in cp["xref"]:
            cp["xref"]["metanetx.chemical"] = [_mnxm]
        return cp, _mnxm

    def mnxr_xref(self, mnxr):
        """Return the cross reference for a given molecule BIGG or MNXM id

        Args:
            mnxr (str): MetaNetX species id
        Returns:
            Dictionnary with the cross reference
        """
        #  logging.debug('------ mnxr_xref ------')
        #  mnxr depr test
        _mnxr = self.single_depr_mnxr(mnxr.upper())
        #  xref
        tmp = self.reac_xref[self.reac_xref["ID"] == mnxr]["source"].to_dict()
        xref = {tmp[i].split(":")[0].lower(): [] for i in tmp if "MNX" not in tmp[i]}
        for i in tmp:
            if "MNX" not in tmp[i]:
                if ":" in tmp[i]:
                    xref[tmp[i].split(":")[0].lower()].append(
                        tmp[i].split(":")[1].lower()
                    )
        #  properties
        cp = self.single_reac_prop(mnxr)
        if cp:
            try:
                xref["ec-code"] = cp["classifs"]
                del cp["classifs"]
            except KeyError:
                pass
        cp["xref"] = xref
        if "metanetx.reaction" not in cp["xref"]:
            cp["xref"]["metanetx.reaction"] = [mnxr]
        else:
            if _mnxr not in cp["xref"]["mnxr"]:
                cp["xref"]["mnxr"].append(mnxr)
        return cp, _mnxr

    #  refresh the cache #

    def refresh_cache(self, delete_old_files=False):
        """Refresh all the cache files"""
        logging.info("Refreshing the cache. This may take a while...")
        if delete_old_files:
            cache_dir = os.path.join(self.base_dir, "flatfiles")
            for f in os.listdir(cache_dir):
                logging.info(
                    "Deleting the following file: " + str(os.path.join(cache_dir, f))
                )
                os.remove(os.path.join(cache_dir, f))
        _ = self.g_depr_mnxm
        _ = self.g_depr_mnxr
        _ = self.chem_xref
        _ = self.reac_xref
        _ = self.chem_prop
        _ = self.reac_prop
        _ = self.biggm_mnxm
        _ = self.biggr_mnxr
        _ = self.mnxm_inchikey
        _ = self.inchikey_mnxm
        _ = self.inchikey2_mnxm
        _ = self.keggm_mnxm
        _ = self.chebim_mnxm
        _ = self.name_mnxm
