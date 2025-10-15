import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)
import json
import logging
import difflib
import os
import tempfile
import tarfile
import time
import urllib.request
from tqdm import tqdm
import re
import gzip
import copy
import pickle
import sys

import compress_json
import networkx as nx
import numpy as np
import pandas as pd
import requests
import pubchempy as pcp
from urllib.error import URLError
from http.client import RemoteDisconnected
from ete3 import NCBITaxa

from typing import Tuple, Dict, List, Union, Optional, Any

from biopathopt.utils import fuzzy_dict_lookup, inchikey_layer_extract, stream_json

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

    def __init__(
            self, 
            low_memory_mode: bool = False, 
            use_progressbar: bool = True
        ):
        self.use_progressbar = use_progressbar
        self.low_memory_mode = low_memory_mode
        self.base_dir = os.path.dirname(os.path.realpath(__file__))
        self._chem_xref = None
        self._reac_xref = None
        self._mnxm_prop = None
        self._mnxr_prop = None
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
        self._brenda_ec_inchikey_kcat = {}
        self._brenda_ec_inchikey_sa = {}
        self._brenda_ec_g = None
        self._rr_prop = None
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
            "MNXM739518": "MNXM1371312",
        }
        self.mnxr_overwrite = {
            "MNXR94682": "MNXR198841",
        }
        self.biggr_overwrite = {
            "GLUDyi": "GLUDy",
        }
        #  other function specfic
        self.pubchem_search_cache = {}
        self.pubchem_min_start = 0.0
        self.pubchem_min_count = 0
        self.pubchem_sec_start = 0.0
        self.pubchem_sec_count = 0
        self.ncbi = None
        self.fuzzy_search_cache = {}
        self.brenda_rest_retries = 50
        """
        Includes the following plus others defined in the ECMpy package
                # get all the hardcoded cofactors that exist
                df_cofactors = pd.read_csv(os.path.join(self.base_dir, 'flatfiles', 'cofactor_inchi_201811.tsv'), sep='\t', comment='#', header=None)
                df_cofactors.columns = ['InChI', 'name', 'MetaNetX']
                self.metanetx_cofactors = [self.single_depr_mnxm(i) for y in df_cofactors.MetaNetX for i in str(y).split(',') if not pd.isna(y)]
        """
        self.mnxm_cofactors = ['MNXM1102419', 'MNXM3', 'MNXM1105762', 'MNXM731165', 'MNXM13', 'MNXM1103285', 'MNXM796', 'MNXM729302', 'MNXM90960', 'MNXM734750', 'MNXM1102191', 'MNXM257', 'MNXM02', 'MNXM1', 'MNXM728294', 'MNXM728062', 'MNXM51918', 'MNXM738702', 'MNXM1107902', 'MNXM360', 'MNXM1104559', 'MNXM739756', 'MNXM1107906', 'MNXM1105927', 'MNXM740692', 'MNXM10', 'MNXM230', 'MNXM232', 'MNXM732398', 'MNXM735437', 'MNXM95', 'MNXM8975', 'MNXM11', 'MNXM255', 'MNXM436', 'MNXM9', 'MNXM731949', 'MNXM178', 'MNXM1102152', 'MNXM124865', 'MNXM2229', 'MNXM128', 'MNXM1101474', 'MNXM2255', 'MNXM39', 'MNXM3654', 'MNXM726339', 'MNXM1108018', 'MNXM1103428', 'MNXM1102167', 'MNXM8978', 'MNXM35', 'MNXM729214', 'MNXM36', 'MNXM729215', 'MNXM924', 'MNXM286', 'MNXM162231', 'MNXM4133', 'MNXM490899', 'MNXM572', 'MNXM1104555', 'MNXM58', 'MNXM731166', 'MNXM741485', 'MNXM1231', 'MNXM732620', 'MNXM1103718', 'MNXM332', 'MNXM40333', 'MNXM727276', 'MNXM1104823', 'MNXM1101285', 'MNXM1103302', 'MNXM726712', 'MNXM730586', 'MNXM266', 'MNXM8', 'MNXM53428', 'MNXM137', 'MNXM5', 'MNXM653', 'MNXM537', 'MNXM452', 'MNXM740736', 'MNXM152', 'MNXM191', 'MNXM1104385', 'MNXM727888', 'MNXM2174', 'MNXM432', 'MNXM411', 'MNXM726711', 'MNXM107', 'MNXM726710', 'MNXM652', 'MNXM736415', 'MNXM1101868', 'MNXM169', 'MNXM27', 'MNXM735438', 'MNXM394', 'MNXM738430', 'MNXM1922', 'MNXM40414', 'MNXM4835', 'MNXM1102128', 'MNXM736654', 'MNXM1105936', 'MNXM1562', 'MNXM234', 'WATER', 'MNXM1107622', 'MNXM735978', 'MNXM4041', 'MNXM1102072', 'MNXM1103458', 'MNXM735047', 'MNXM727224', 'MNXM733186', 'MNXM1103553', 'MNXM344', 'MNXM1092965']

    # ############ 

    def flush_parameters(self):
        self._chem_xref = None
        self._reac_xref = None
        self._mnxm_prop = None
        self._mnxr_prop = None
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
        self._brenda_ec_inchikey_kcat = {}
        self._brenda_ec_inchikey_sa = {}
        self._brenda_ec_g = None
        self._rr_prop = None

    ###########################
    ######### PARSE ###########
    ###########################

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
        """Parses a metabolic network equation string from MetaNetX or RetroRules into reactants and products.

        This function takes a metabolic network equation string and parses it into
        lists of reactants and products. The input string should be in the form
        'reactant1 + reactant2 = product1 + product2', where each component may
        contain quantities and compartment information separated by spaces and '@' symbols.

        Args:
            reac_str (str): The metabolic network equation string.

        Returns:
            Tuple[List[Tuple], List[Tuple]]: A tuple containing two lists:
                - The first list contains lists of parsed reactants.
                - The second list contains lists of parsed products.
        """
        # Rescue map for symbolic stoichiometries
        stoichio_rescue = {
            '4n': 4, '3n': 3, '2n': 2, 'n': 1, '(n)': 1, '(N)': 1, '(2n)': 2, '(x)': 1,
            'N': 1, 'm': 1, 'q': 1, '0.01': 1, '0.1': 1, '0.5': 1, '1.5': 1,
            '0.02': 1, '0.2': 1, '(n-1)': 0, '(n-2)': -1
        }
        # Example chunk: "2 C00001@MNXM" -> ("2", "C00001")
        single_reac_re = re.compile(r'(\(n-1\)|\d+|4n|3n|2n|n|\(n\)|\(N\)|\(2n\)|\(x\)|N|m|q|\(n\-2\)|\d+\.\d+) ([\w\d]+)@\w+')
        
        def parse_reaction_side(eq_side: str) -> Dict[str, int]:
            """Parse reaction side"""
            out = []
            for sto, mnxm in single_reac_re.findall(eq_side):
                mnxm = self.single_depr_mnxm(mnxm.strip())
                try:
                    out.append(
                        (
                            stoichio_rescue.get(sto.strip(), int(sto.strip())),
                            mnxm
                        )
                    )
                except ValueError:
                    logging.warning(f"Cannot convert stoichiometry {sto} from {mnxm}")
            return out
        
        # Split the equation string into reactants and products
        reactants_str, products_str = reac_str.split('=')
        reactants = parse_reaction_side(reactants_str)
        products = parse_reaction_side(products_str)

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
                return self.mnxr_overwrite[mnxr]
            except KeyError:
                pass
            try:
                tmp = []
                for i in list(nx.dfs_tree(self.g_depr_mnxr, source=mnxr)):
                    if not list(self.g_depr_mnxr.successors(i)):
                        tmp.append(i)
                if len(tmp) == 1:
                    if tmp[0] != "EMPTY":
                        return tmp[0]
                elif len(tmp) > 1:
                    logging.warning(
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
                        logging.warning(
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

    def single_mnxm_prop(self, mnxm):
        """Unify the search of MetaNetX molecule id's. First assume that
        the id is correct, then check that its not an old ID and if not,
        check the deprecated maps and return the closest. Because of
        protonation you may have duplicates and thus the deprecated map.

        Args:
            mnxm (str): MetaNetX molecule ID
        Returns:
            str: The non-deprecated MetaNetX molecule id
        """
        #  logging.debug('------ single_mnxm_prop -------')
        try:
            return self.mnxm_prop[mnxm]
        except KeyError:
            try:
                return self.mnxm_prop[self.single_depr_mnxm(mnxm)]
            except KeyError:
                try:
                    search_mnxm = self.single_depr_mnxm(mnxm)
                except KeyError:
                    search_mnxm = mnxm
                #  if all else fails, find the closest member of the deprecated map
                try:
                    depr_lst = list(nx.bfs_tree(
                            self.g_depr_mnxm, search_mnxm, reverse=True
                        ))
                    for i in depr_lst:
                        try:
                            return self.mnxm_prop[i]
                        except KeyError:
                            pass
                except nx.exception.NetworkXError as e:
                    pass
        logging.warning("Cannot find chem properties for: " + str(mnxm))
        return {}

    def single_mnxr_prop(self, mnxr):
        """Unify the search of MetaNetX reaction id's. First assume that
        the id is correct, then check that its not an old ID. Because of
        protonation you may have duplicates and thus the deprecated map.

        Args:
            mnxm (str): MetaNetX reaction ID
        Returns:
            str: The non-deprecated MetaNetX reaction id
        """
        logging.debug("------ single_mnxr_prop ------")
        try:
            return self.mnxr_prop[mnxr]
        except KeyError:
            try:
                return self.mnxr_prop[self.single_depr_mnxr(mnxr)]
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

    def _get_protein_reactants_inchikey(self, ec_entry, protein_ids=[]):
        #because the comments do not always have the right substrates, get the original reaction from the protein

        def extract_reactants_inchikey(reactions_list, cofactor_list=[], filter_proteins_ids=[]):
            reactant_dict = {}
            logging.debug(f'extract_reactants_inchikey: {len(reactions_list)}')
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
            if 'substrates_products' in ec_entry:
                if protein_ids:
                    try:
                        reactions_list = [y for y in ec_entry['substrates_products'] if set(protein_ids) & set(y['proteins'])]
                    except KeyError:
                        reactions_list = ec_entry['substrates_products']
                else:
                    reactions_list = ec_entry['substrates_products']
        if not reactions_list:
            if 'reaction' in ec_entry:
                if protein_ids:
                    try:
                        reactions_list = [y for y in ec_entry['reaction'] if set(protein_ids) & set(y['proteins'])]
                    except KeyError:
                        reactions_list = ec_entry['reaction']
                else:
                    reactions_list = ec_entry['reaction']
        if 'cofactor' in ec_entry:
            if protein_ids:
                try:
                    cofactor_list = [y for y in ec_entry['cofactor'] if set(protein_ids) & set(y['proteins'])]
                except KeyError:
                    cofactor_list = ec_entry['cofactor']
            else:
                cofactor_list = ec_entry['cofactor']
        else:
            cofactor_list = []
        logging.debug('extract_reactants_inchikey')
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

    def _generate_brenda_kinetics(
            self, 
            parse_type='kcat', 
            use_progressbar: bool = False, 
        ):
        # check that the brenda file is there
        brenda_path_file = os.path.join(
            self.base_dir, "flatfiles/brenda.json.tar.gz"
        )
        if not os.path.exists(brenda_path_file):
            raise TypeError('Need to download the file at: https://brenda-enzymes.org/download.php and save it as brenda.json.tar.gz in biopathopt/flatfiles')
        if parse_type not in ['kcat', 'sa']:
            raise TypeError(f'Input Must be kcat or sa: {parse_type}')
        pbar = None
        last_read_bytes = 0
        if use_progressbar:
            total_bytes = os.path.getsize(brenda_path_file)
            pbar = tqdm(total=total_bytes, unit="B", unit_scale=True, desc=f"Generate Brenda {parse_type} File")
        for ec_number, ec_entry, bytes_read in stream_json(brenda_path_file, pointer='data'):
            logging.debug(f'------------ {ec_number} ---------')
            res = {}
            if use_progressbar:
                pbar.set_description(f"Generate Brenda {parse_type} File: {ec_number}")
                delta = bytes_read - last_read_bytes
                if delta > 0:
                    pbar.update(delta)
                    last_read_bytes = bytes_read
            if parse_type=='kcat':
                if ec_number in self._brenda_ec_inchikey_kcat:
                    logging.debug(f'Skipping {ec_number}')
                    continue
            elif parse_type=='sa':
                if ec_number in self._brenda_ec_inchikey_sa:
                    logging.debug(f'Skipping {ec_number}')
                    continue
            logging.debug('deep copy')
            if parse_type=='sa':
                try:
                    kinetics_list = copy.deepcopy(ec_entry['specific_activity'])
                except KeyError:
                    continue
                try:
                    protein_dict = ec_entry['protein']
                except KeyError:
                    continue
            elif parse_type=='kcat':
                try:
                    kinetics_list = copy.deepcopy(ec_entry['kcat_km_value'])
                except KeyError:
                    continue
                try:
                    protein_dict = ec_entry['protein']
                except KeyError:
                    continue
            ### get the reactants from the reaction description
            logging.debug('_get_protein_reactants_inchikey')
            try:
                default_reactants = self._get_protein_reactants_inchikey(
                    ec_entry=ec_entry
                )
            except KeyError:
                default_reactants = {}
            #res[ec_number] = {}
            logging.debug('looping through the kinetics list')
            logging.debug(len(kinetics_list))
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
                    try:
                        value = float(kinetic_entry.get('value'))
                        if value<=0.0:
                            continue
                    except ValueError:
                        continue
                elif parse_type=='kcat':
                    value, extracted_sub = self._extract_brenda_substrates_kcat(kinetic_entry.get('value').lower())
                    if not value and value>0.0:
                        continue
                if extracted_sub:
                    logging.debug('molecule_name_search_inchikey')
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
                if not sub_id in res:
                    res[sub_id] = {}
                for org in to_add_org:
                    if org not in res[sub_id]:
                        res[sub_id][org] = {'uniprot': [], 'values': [], 'comments': []}
                    res[sub_id][org]['values'].append(value)
                    res[sub_id][org]['uniprot'] += to_add_uniprot
                    res[sub_id][org]['comments'].append(kinetic_entry.get('comment').lower())
            if parse_type=='kcat':
                self._brenda_ec_inchikey_kcat[ec_number] = res
            elif parse_type=='sa':
                self._brenda_ec_inchikey_sa[ec_number] = res
        if last_read_bytes < total_bytes:
            pbar.update(total_bytes - last_read_bytes)

    ##### taxonomy ###

    def _get_species_name(
        self,
        taxid: int,
    ) -> str:
        """Return the species-level name from a taxonomy ID."""
        if not self.ncbi:
            self.ncbi = NCBITaxa()
        lineage = ncbi.get_lineage(taxid)
        names = ncbi.get_taxid_translator(lineage)
        ranks = ncbi.get_rank(lineage)

        for tid in lineage:
            if ranks[tid] == "species":
                return names[tid]
        return names.get(taxid, "Unknown")

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
    '''

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

    # ###########################################
    # ###### PROPERTIES #########################
    # ###########################################
    #  These are designed to behave as parameters. The first time its called it
    #  will load, and the next time around it will pass the saved parameter

    # ############# BRENDA ######################

    @property
    def brenda_ec_g(self):
        if not self._brenda_ec_g:
            logging.debug('-------- brenda_ec_g ---------')
            logging.debug('Populating.....')
            path_file = os.path.join(
                self.base_dir, "flatfiles/brenda_ec_g.pkl"
            )
            if os.path.exists(path_file):
                with open(path_file, 'rb') as f:
                    self._brenda_ec_g = pickle.load(f)
            else:
                ### generate it ####
                def extract_history_ec(text):
                    match = re.search(r'EC\s+(\d+\.\d+\.\d+\.\d+)', text)
                    if match:
                        ec_number = match.group(1)
                        return ec_number
                    return None
                G = nx.DiGraph()
                # check that the brenda file is there
                brenda_path_file = os.path.join(
                    self.base_dir, "flatfiles/brenda.json.tar.gz"
                )
                if not os.path.exists(brenda_path_file):
                    raise TypeError('Need to download the file at: https://brenda-enzymes.org/download.php and save it as brenda.json.tar.gz in biopathopt/flatfiles')
                for ec_number, ec_entry, bytes_read in stream_json(brenda_path_file, pointer='data'):
                    G.add_node(ec_number)    
                for ec_number, ec_entry, bytes_read in stream_json(brenda_path_file, pointer='data'):
                    if 'history' in ec_entry:
                        to_ec = extract_history_ec(
                            ec_entry['history']
                        )
                        if to_ec:
                            G.add_edge(to_ec, ec_number)

                with open(path_file, 'wb') as f:
                    pickle.dump(G, f, pickle.HIGHEST_PROTOCOL)
                self._brenda_ec_g = G
        return self._brenda_ec_g

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
                for attempt in range(1, self.brenda_rest_retries + 1):
                    try:
                        self._generate_brenda_kinetics(
                                parse_type='kcat', 
                                use_progressbar=self.use_progressbar, 
                            )
                        break  # success, exit loop
                    except (pcp.PubChemHTTPError, URLError, RemoteDisconnected) as e:
                        logging.warning(f'The following eror: {e}... retrying')
                        if attempt >= self.brenda_rest_retries:
                            raise e # give up after last attempt
                        time.sleep(5.0)
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
                for attempt in range(1, self.brenda_rest_retries + 1):
                    try:
                        self._generate_brenda_kinetics(
                                parse_type='sa', 
                                use_progressbar=self.use_progressbar, 
                            )
                        break  # success, exit loop
                    except (pcp.PubChemHTTPError, URLError, RemoteDisconnected) as e:
                        logging.warning(f'The following eror: {e}... retrying')
                        if attempt >= self.brenda_rest_retries:
                            raise e# give up after last attempt
                        time.sleep(5.0)
                with open(path_file, 'wb') as f:
                    pickle.dump(self._brenda_ec_inchikey_sa, f, pickle.HIGHEST_PROTOCOL)
        return self._brenda_ec_inchikey_sa

    ###### RETRORULES #########

    @property
    def retrorules_prop(self):
        """Return the properties of retrorules2 

        This function will download the following file 
        https://zenodo.org/record/5828017/files/retrorules_rr02_rp2_hs.tar.gz
        that describe the retrorules properties of rules used byt rp2

        Args:
        Returns:
            dict: RetroRules rules properties
        """
        def _same_or_zero(lst: list[int]) -> int:
            """Return the integer if both elements in list are equal, else 0."""
            return lst[0] if len(lst) == 2 and lst[0] == lst[1] else 0

        if not self._rr_prop:
            logging.debug("------ retrorules_prop ----")
            logging.debug("\t-> Populating...")
            rr_prop_path = os.path.join(self.base_dir, "flatfiles/rr_prop.json.gz")
            if not os.path.exists(rr_prop_path):
                logging.info("RetroRules 2 prop does not exist... generating it")
                with tempfile.TemporaryDirectory() as tmpdirname:
                    rr_prop_ori_path = os.path.join(tmpdirname, "rr_prop.tsv")
                    tmp_tar_path = os.path.join(tmpdirname, 'tmp.tar.gz')
                    urllib.request.urlretrieve(
                        "https://zenodo.org/record/5828017/files/retrorules_rr02_rp2_hs.tar.gz",
                        tmp_tar_path,
                    )
                    with tarfile.open(tmp_tar_path, "r:gz") as tar:
                        tar.extractall(path=tmpdirname)
                    rr_prop = pd.read_csv(
                        os.path.join(tmpdirname, "retrorules_rr02_rp2_hs", "retrorules_rr02_rp2_flat_all.csv")
                    )
                    rr_prop['mnxr'] = [self.single_depr_mnxr(i.split('_')[0]) for i in rr_prop['Legacy ID']]
                    rr_prop['mnxm'] = [self.single_depr_mnxm(i.split('_')[1]) for i in rr_prop['Legacy ID']]
                    retrorules_prop = {}
                    for index, r in rr_prop.iterrows():
                        #merge those that are the same
                        if not r['Rule ID'] in retrorules_prop:
                            retrorules_prop[r['Rule ID']] = {}
                        if not r['mnxr'] in retrorules_prop[r['Rule ID']]:
                            retrorules_prop[r['Rule ID']][r['mnxr']] = {}
                        if r['mnxm'] in retrorules_prop[r['Rule ID']][r['mnxr']]:
                            retrorules_prop[r['Rule ID']][r['mnxr']][r['mnxm']]['ec'] = list(set(
                                    [y for y in r['EC number'].split(';') if y!='NOEC']+retrorules_prop[r['Rule ID']][r['mnxr']][r['mnxm']]['ec']
                                ))
                            retrorules_prop[r['Rule ID']][r['mnxr']][r['mnxm']]['rule_score'] = float(np.mean([
                                    r['Score normalized'],
                                    float(retrorules_prop[r['Rule ID']][r['mnxr']][r['mnxm']]['rule_score']),
                                ]))
                            retrorules_prop[r['Rule ID']][r['mnxr']][r['mnxm']]['direction'] = _same_or_zero(
                                    [
                                        int(r['Rule relative direction']), 
                                        retrorules_prop[r['Rule ID']][r['mnxr']][r['mnxm']]['direction'],
                                    ]
                                )
                        else:
                            retrorules_prop[r['Rule ID']][r['mnxr']][r['mnxm']] = {
                                'ec': [y for y in r['EC number'].split(';') if y!='NOEC'],
                                'rule_score': float(r['Score normalized']),
                                'direction': int(r['Rule relative direction']),
                            }
                    compress_json.dump(retrorules_prop, rr_prop_path)
                    rr_prop = None
            self._rr_prop = compress_json.load(rr_prop_path)
        return self._rr_prop


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
    def mnxm_prop(self):
        """Return the chemical properties of molecules

        This function will download the following file https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_prop.tsv that
        describes the chemical structure, etc...

        Args:
        Returns:
            dict: Chemical properties from MetaNetX
        """
        if not self._mnxm_prop:
            logging.debug("------ mnxm_prop ----")
            logging.debug("\t-> Populating...")
            mnxm_prop_path = os.path.join(self.base_dir, "flatfiles/mnxm_prop.json.gz")
            if not os.path.exists(mnxm_prop_path):
                logging.info("Chem prop does not exist... generating it")
                with tempfile.TemporaryDirectory() as tmpdirname:
                    mnxm_prop_ori_path = os.path.join(tmpdirname, "mnxm_prop.tsv")
                    urllib.request.urlretrieve(
                        "https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_prop.tsv",
                        mnxm_prop_ori_path,
                    )
                    mnxm_prop = pd.read_csv(
                        mnxm_prop_ori_path, comment="#", sep="\t", header=None
                    )
                    #  Need to do this since id's are now mutiple and contain conflicts
                    #  mnxm_prop["ID"] = [single_depr_mnxm(i) for i in mnxm_prop["ID"]]
                    mnxm_prop.columns = [
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
                    mnxm_prop["InChIKey"] = [
                        i.replace("InChIKey=", "") if not pd.isna(i) else i
                        for i in mnxm_prop["InChIKey"]
                    ]
                    mnxm_prop = mnxm_prop.set_index("ID")
                    mnxm_prop = mnxm_prop.transpose().to_dict()
                    compress_json.dump(mnxm_prop, mnxm_prop_path)
                    mnxm_prop = None
            self._mnxm_prop = compress_json.load(mnxm_prop_path)
        return self._mnxm_prop

    @property
    def mnxr_prop(self):
        """Return the chemical properties of molecules

        This function will download the following file https://www.metanetx.org/cgi-bin/mnxget/mnxref/mnxr_prop.tsv that
        describes the chemical structure, etc...

        Args:
        Returns:
            dict: Reaction properties from MetaNetX
        """
        if not self._mnxr_prop:
            logging.debug("------ mnxr_prop -----")
            logging.debug("\t-> Populating...")
            mnxr_prop_path = os.path.join(self.base_dir, "flatfiles/mnxr_prop.json.gz")
            if not os.path.exists(mnxr_prop_path):
                with tempfile.TemporaryDirectory() as tmpdirname:
                    mnxr_prop_ori_path = os.path.join(tmpdirname, "mnxr_prop.tsv")
                    if not os.path.exists(mnxr_prop_ori_path):
                        logging.info(
                            "The MetaNetX file mnxr_prop.tsv does not exist... downloading it"
                        )
                        urllib.request.urlretrieve(
                            "https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_prop.tsv",
                            mnxr_prop_ori_path,
                        )
                    reacs_prop = pd.read_csv(
                        mnxr_prop_ori_path,
                        comment="#",
                        sep="\t",
                        header=None,
                    )
                    reacs_prop.columns = [
                        "ID",
                        "mnx_equation",
                        "reference",
                        "classifs",
                        "is_balanced",
                        "is_transport",
                    ]
                    reacs_prop = reacs_prop.set_index("ID")
                    reacs_prop = reacs_prop.transpose().to_dict()
                    for i in reacs_prop:
                        try:
                            reacs_prop[i]["ec"] = reacs_prop[i]["classifs"].split(
                                ";"
                            )
                        except AttributeError:
                            reacs_prop[i]["ec"] = []
                    ### generate the reaction based in inchikeys
                    for i in reacs_prop:
                        try:
                            reacs_prop[i]['inchikey2_equation'] = \
                                self.convert_mnxr_equation(reacs_prop[i]['mnx_equation'], inchikey_levels=2)
                        except (ValueError, KeyError) as e:
                            reacs_prop[i]['inchikey2_equation'] = ''
                        try:
                            reacs_prop[i]['inchikey_equation'] = \
                                self.convert_mnxr_equation(reacs_prop[i]['mnx_equation'], inchikey_levels=3)
                        except (ValueError, KeyError) as e:
                            reacs_prop[i]['inchikey_equation'] = ''
                    ### generate main left and main right
                    for i in reacs_prop:
                        reactants, products = self.parse_mnxr_equation(reacs_prop[i]['mnx_equation'])
                        reacs_prop[i]['main_reactants'] = {y[1]: y[0] for y in reactants if y[1] not in self.mnxm_cofactors}
                        reacs_prop[i]['secondary_reactants'] = {y[1]: y[0] for y in reactants if y[1] in self.mnxm_cofactors}
                        reacs_prop[i]['main_products'] = {y[1]: y[0] for y in products if y[1] not in self.mnxm_cofactors}
                        reacs_prop[i]['secondary_products'] = {y[1]: y[0] for y in products if y[1]  in self.mnxm_cofactors}
                    #save it
                    compress_json.dump(reacs_prop, mnxr_prop_path)
            self._mnxr_prop = compress_json.load(mnxr_prop_path)
        return self._mnxr_prop

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
            for i in self.mnxm_prop:
                try:
                    if not pd.isna(self.mnxm_prop[i]["InChIKey"]) or not self.mnxm_prop[i]["InChIKey"]=='nan':
                        self._inchikey2_mnxm[
                            inchikey_layer_extract(self.mnxm_prop[i]['InChIKey'], 2)
                        ] = self.single_depr_mnxm(i)
                except (KeyError, AttributeError) as e:
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
            for i in self.mnxm_prop:
                try:
                    if not pd.isna(self.mnxm_prop[i]["InChIKey"]):
                        self._inchikey_mnxm[
                            self.mnxm_prop[i]["InChIKey"]
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
            for i in self.mnxm_prop:
                try:
                    if not pd.isna(self.mnxm_prop[i]["InChIKey"]):
                        self._mnxm_inchikey[
                                self.single_depr_mnxm(i)
                        ] = self.mnxm_prop[i]["InChIKey"]
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
            self._keggr_mnxr = {i: self.single_depr_mnxr(str(tmp[i][0])) for i in tmp}
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
            self._keggm_mnxm = {i: self.single_depr_mnxm(str(tmp[i][0])) for i in tmp}
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
            self._chebim_mnxm = {i: self.single_depr_mnxm(str(tmp[i][0])) for i in tmp}
        return self._chebim_mnxm

    @property
    def molname_mnxm(self):
        """molecule name to mnxm
        """
        if not self._molname_mnxm:
            logging.debug("------ molname_mnxm ------")
            logging.debug("\t-> Populating...")
            self._molname_mnxm = {}
            for mnxm in self.mnxm_prop:
                if 'name' in self.mnxm_prop[mnxm]:
                    if not pd.isna(self.mnxm_prop[mnxm]['name']):
                        if self.mnxm_prop[mnxm]['name'] not in self._molname_mnxm:
                            self._molname_mnxm[self.mnxm_prop[mnxm]['name']] = mnxm
                        #else:
                        #    logging.warning(f"Duplicate name entries for {self.mnxm_prop[mnxm]['name']}")
        return self._molname_mnxm

    # ################# search functions ###############

    # #### Pubchem #####

    def exact_pubchem_search(self, query: str, itype: str = 'name', return_lowest_cid: bool = False) -> Dict[str, Any]:
        """
        Perform an exact search on PubChem using the given identifier.

        Args:
            query (str): The compound name or identifier to search.
            itype (str): The type of identifier (e.g., 'name', 'smiles', 'inchi', 'inchikey'). Defaults to 'name'.

        Returns:
            Dict[str, Any]: A dictionary containing the compound data if found. Empty dict if not found or ambiguous.

        Raises:
            KeyError: If multiple compounds are returned for the query.
        """
        cid = {}
        cid_keys = [
            "canonical_smiles",
            "charge",
            "cid",
            "elements",
            "exact_mass",
            "inchi",
            "inchikey",
            "isomeric_smiles",
            "iupac_name",
            "molecular_formula",
            "molecular_weight",
        ]
        logging.debug(f'Query: {query} - itype: {itype}')
        
        if not query:
            return {}

        # Check if result is already cached
        try:
            if query.lower() in self.pubchem_search_cache:
                logging.debug(f'Query is in pubchem cache: {self.pubchem_search_cache[query.lower()]}')
            return self.pubchem_search_cache[query.lower()]
        except KeyError:
            pass

        try:
            # Perform the PubChem search
            logging.debug(f'Searching for {query}')
            cids = pcp.get_compounds(identifier=query, namespace=itype)
        except ValueError:
            # Cache the empty result and return
            self.pubchem_search_cache[query.lower()] = {}
            return {}

        if len(cids) == 0:
            logging.debug('There are not results')
            return {}
        elif len(cids) == 1:
            logging.debug('There is one single result')
            # Cache and return the single result
            cid = cids[0].to_dict()
            cid = {i: cid.get(i) for i in cid_keys}
        else:
            logging.warning(f'There are multiple results for query: {query} - {itype}')
            if return_lowest_cid:
                cid = min(cids, key=lambda x: int(x.cid))
                cid = cid.to_dict()
                cid = {i: cid.get(i) for i in cid_keys}
            else:
                # Cache the empty result and raise an error for ambiguity
                self.pubchem_search_cache[query.lower()] = {}
                raise KeyError(f'Multiple cids {cids} for {query}')
        #### xref ####
        if 'cid' in cid:
            r = requests.post(f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid.get("cid")}/xrefs/SBURL/JSON')
            res_list = r.json()
            xref = {}
            xref['pubchem'] = [str(cid)]
            for url in res_list['InformationList']['Information'][0]['SBURL']:
                if 'https://biocyc.org/compound?orgid=META&id=' in url:
                    if 'biocyc' not in xref:
                        xref['biocyc'] = []
                    xref['biocyc'].append(url.replace('https://biocyc.org/compound?orgid=META&id=', ''))
                if 'http://www.hmdb.ca/cidbolites/' in url:
                    if 'hmdb' not in xref:
                        xref['hmdb'] = []
                    xref['hmdb'].append(url.replace('http://www.hmdb.ca/cidbolites/', ''))
                if 'http://www.genome.jp/dbget-bin/www_bget?cpd:' in url:
                    if 'kegg.compound' not in xref:
                        xref['kegg.compound'] = []
                    xref['kegg.compound'].append(url.replace('http://www.genome.jp/dbget-bin/www_bget?cpd:', ''))
                if 'http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:' in url:
                    if 'chebi' not in xref:
                        xref['chebi'] = []
                    xref['chebi'].append('CHEBI:'+url.replace('http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:', ''))
                    xref['chebi'].append(url.replace('http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:', ''))
            xref['inchi'] = cid.get('inchi')
            xref['inchi_key'] = cid.get('inchikey')
            xref['smiles'] = cid.get('canonical_smiles')
            cid['xref'] = xref
        self.pubchem_search_cache[query.lower()] = cid
        return cid
    

    def molecule_name_search_inchikey(self, name: str) -> Optional[str]:
        """
        Search for the InChIKey of a molecule based on its name using cached data,
        fuzzy matching with MNX identifiers, and PubChem search.

        Args:
            name (str): The name of the molecule to search.

        Returns:
            Optional[str]: The InChIKey of the molecule if found, otherwise None.
        """
        #try:
        # Attempt lookup in cached PubChem search results
        logging.debug(f'--------- Searching for {name} -----------------')
        if not name:
            logging.debug('The query is empty')
            return None
        fixed_name_mnxm = {
                'nad': 'MNXM8', 
                'nadh': 'MNXM10', 
                'h': 'MNXM1', 
                'h+': 'MNXM1', 
                'h2o': 'WATER',
                'atp': 'MNXM3',
            }
        if name.lower() in fixed_name_mnxm:
            logging.debug('name is in fixed list')
            return self.mnxm_inchikey[fixed_name_mnxm[name.lower()]]
        try:
            cid = self.pubchem_search_cache[name.lower()]
            inchikey = cid.get('inchikey')
            if inchikey:
                return inchikey
        except KeyError:
            pass
        # Fuzzy match against internal MNX mapping
        if name.lower() in self.fuzzy_search_cache:
            logging.debug('query is in fuzzy_search_cache: {fuzzy_mnxm}')
            fuzzy_mnxm = self.fuzzy_search_cache[name.lower()]
            try:
                return self.mnxm_inchikey[fuzzy_mnxm]
            except KeyError:
                pass
        else:
            fuzzy_mnxm = fuzzy_dict_lookup(name.lower(), self.molname_mnxm)
            self.fuzzy_search_cache[name.lower()] = fuzzy_mnxm
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

        #except Exception:
        #    pass

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
                try:
                    mnxm = self.inchikey2_mnxm[inchikey_layer_extract(inchikey, 2)]
                    return self.mnxm_xref(mnxm)
                except KeyError:
                    pass
        logging.debug("Cannot find the follwing inchikey: " + str(inchikey))
        return {}, None

    def mnxm_xref(self, mnxm):
        """Return the cross reference for a MetaNetX molecule ID

        TODO: change this to mnxm_description that includes the xref
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
        cp = self.single_mnxm_prop(_mnxm)
        cp["xref"] = xref
        if not cp:
            raise KeyError(f'Cannot find the cross reference for {mnxm}')
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
        tmp = self.reac_xref[self.reac_xref["ID"] == _mnxr]["source"].to_dict()
        xref = {tmp[i].split(":")[0].lower(): [] for i in tmp if "MNX" not in tmp[i]}
        for i in tmp:
            if "MNX" not in tmp[i]:
                if ":" in tmp[i]:
                    xref[tmp[i].split(":")[0].lower()].append(
                        tmp[i].split(":")[1].lower()
                    )
        #  properties
        cp = self.single_mnxr_prop(_mnxr)
        if cp:
            try:
                xref["ec-code"] = cp["classifs"]
                del cp["classifs"]
            except KeyError:
                pass
        cp["xref"] = xref
        if "metanetx.reaction" not in cp["xref"]:
            cp["xref"]["metanetx.reaction"] = [_mnxr]
        else:
            if _mnxr not in cp["xref"]["mnxr"]:
                cp["xref"]["mnxr"].append(mnxr)
        return cp, _mnxr

    ################  refresh the cache ###############3

    def refresh_cache(
            self, 
            delete_old_files=False, 
            use_progressbar=True,
            parse_brenda_files=True,
            ):
        """Refresh all the cache files"""
        logging.info("Refreshing the cache. This may take a while...")
        self.use_progressbar = use_progressbar
        if delete_old_files:
            cache_dir = os.path.join(self.base_dir, "flatfiles")
            for f in os.listdir(cache_dir):
                logging.info(
                    "Deleting the following file: " + str(os.path.join(cache_dir, f))
                )
                os.remove(os.path.join(cache_dir, f))
        pbar = None
        if use_progressbar:
            if parse_brenda_files:
                pbar = tqdm(total=14)
            else:
                pbar = tqdm(total=15)
            pbar.set_description(f"Processing g_depr_mnxm")
        _ = self.g_depr_mnxm
        if use_progressbar:
            pbar.update(1)
            pbar.set_description(f"Processing g_depr_mnxr")
        _ = self.g_depr_mnxr
        if use_progressbar:
            pbar.update(1)
            pbar.set_description(f"Processing chem_xref")
        _ = self.chem_xref
        if use_progressbar:
            pbar.update(1)
            pbar.set_description(f"Processing reac_xref")
        _ = self.reac_xref
        if use_progressbar:
            pbar.update(1)
            pbar.set_description(f"Processing mnxm_prop")
        _ = self.mnxm_prop
        if use_progressbar:
            pbar.update(1)
            pbar.set_description(f"Processing mnxr_prop")
        _ = self.mnxr_prop
        if use_progressbar:
            pbar.update(1)
            pbar.set_description(f"Processing biggm_mnxm")
        _ = self.biggm_mnxm
        if use_progressbar:
            pbar.update(1)
            pbar.set_description(f"Processing biggr_xref")
        _ = self.biggr_mnxr
        if use_progressbar:
            pbar.update(1)
            pbar.set_description(f"Processing mnxm_inchikey")
        _ = self.mnxm_inchikey
        if use_progressbar:
            pbar.update(1)
            pbar.set_description(f"Processing inchikey_mnxm")
        _ = self.inchikey_mnxm
        if use_progressbar:
            pbar.update(1)
            pbar.set_description(f"Processing inchikey2_mnxm")
        _ = self.inchikey2_mnxm
        if use_progressbar:
            pbar.update(1)
            pbar.set_description(f"Processing keggm_mnxm")
        _ = self.keggm_mnxm
        if use_progressbar:
            pbar.update(1)
            pbar.set_description(f"Processing chebim_mnxm")
        _ = self.chebim_mnxm
        if use_progressbar:
            pbar.update(1)
            pbar.set_description(f"Processing molname_mnxm")
        _ = self.molname_mnxm
        if use_progressbar:
            pbar.update(1)
            pbar.set_description(f"Processing retrorules")
        _ = self.retrorules_prop
        if parse_brenda_files:
            if use_progressbar:
                pbar.update(1)
                pbar.set_description(f"Processing brenda_ec_g")
            _ = self.brenda_ec_g
            if use_progressbar:
                pbar.close()
            _ = self.brenda_ec_inchikey_sa
            _ = self.brenda_ec_inchikey_kcat
        else:
            if use_progressbar:
                pbar.close()

