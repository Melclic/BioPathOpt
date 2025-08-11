from cobra import Model, Reaction, Metabolite
#from cobra.io import read_sbml_model, save_json_model, load_json_model, write_sbml_model
import cobra
from typing import Any, List, Tuple, Dict, Union, Optional, Callable

import pandas as pd
import numpy as np
import copy
import os
import gzip
import itertools
import logging
from biopathopt import Data
from biopathopt import utils
from tqdm import tqdm

"""
This is a collection of functions to build a model
"""


class ModelBuilder(Data):

    def __init__(self, path_to_model=None, use_progressbar=False):
        """Class that inherits Data used to build a cobra model
        """
        super().__init__()
        # Manually created mapping to convert annotation IDs from MetaNetX to BiGG or other standards
        self.mnxm_bigg_annot_convert = {
            'CHEBI': 'chebi',
            'metacyc.compound': 'metacyc.compound',
            'kegg.glycan': 'kegg.compound',
            'metacycM': 'metacyc.compound',
            'seedM': 'seed.compound',
            'kegg.drug': 'kegg.compound',
            'rheaP': 'rhea.compound',
            'chebi': 'chebi',
            'rheaG': 'rhea.compound',
            'hmdb': 'hmdb',
            'envipath': 'envipath',
            'envipathM': 'envipath',
            'keggC': 'kegg.compound',
            'sabiork.compound': 'sabiork',
            'seed.compound': 'seed.compound',
            'reactome': 'reactome.compound',
            'lipidmaps': 'lipidmaps',
            'mnx': 'metanetx.chemical',
            'bigg.metabolite': 'bigg.metabolite',
            'slm': 'slm',
            'biggM': 'bigg.metabolite',
            'keggD': 'kegg.compound',
            'keggG': 'kegg.compound',
            'kegg.compound': 'kegg.compound',
            'reactomeM': 'reactome.compound',
            'sabiorkM': 'sabiork',
            'lipidmapsM': 'lipidmaps',
            'SLM': 'slm',
            'InChI': 'inchi',
            'SMILES': 'smiles',
            'InChIKey': 'inchi_key,'
        }

        # Manually created mapping to convert annotation IDs from MetaNetX to BiGG
        self.mnxr_bigg_annot_convert = {
            'metacyc.reaction': 'metacyc.reaction',
            'rheaR': 'rhea',
            'seedR': 'seed.reaction',
            'mnx': 'metanetx.reaction',
            'mnxr': 'metanetx.reaction',
            'kegg.reaction': 'kegg.reaction',
            'rhea': 'rhea',
            'keggR': 'kegg.reaction',
            'keggr': 'kegg.reaction',
            'sabiork.reaction': 'sabiork',
            'sabiorkR': 'sabiork',
            'seed.reaction': 'seed.reaction',
            'biggR': 'bigg.reaction',
            'biggr': 'bigg.reaction',
            'bigg.reaction': 'bigg.reaction',
            'metacycR': 'metacyc.reaction',
            'ec': 'ec-code',
        }

        self.gene_bigg_annot_convert = {
            'UniProtKB': 'uniprot',
        }

        if path_to_model:
            self._read_model(path_to_model, use_progressbar)
        else:
            self.model = Model()


    def _read_model(self, path_to_model: str, use_progressbar: bool = False):
        #read if gzip , sbml or json
        if path_to_model.endswith('json.gz') or path_to_model.endswith('json.gzip'):
            with gzip.open(path_to_model, 'rt', encoding='utf-8') as gz_file:
                self.model = cobra.io.load_json_model(gz_file)
        elif path_to_model.endswith('json'):
            self.model = cobra.io.load_json_model(path_to_model)
        elif path_to_model.endswith('sbml') or path_to_model.endswith('xml'):
            self.model = cobra.io.read_sbml_model(path_to_model)
        else:
            raise TypeError('Cannot recognise the input file extension')
        if not 'biopathopt_enriched' in self.model.annotation:
            #replace deprecated metanetx ids
            logging.info('Replacing deprecated MXNM')
            self._replace_depr_mnxm()
            self._replace_dpr_mnxr()
            #update the metabolite annotation keys
            #TODO
            #update the reactions annotation keys
            #TODO
            #update the metabolite annotations
            met_iterator = tqdm(self.model.metabolites, desc='Updating the metabolite annotations') if use_progressbar else self.model.metabolites
            for m in met_iterator:
                m.annotation = utils.merge_annot_dicts(m.annotation, self._find_metabolite_xref(m.annotation))
            #BUG: remove nan or None
            for m in self.model.reactions:
                annot = m.annotation
                updated_annot = copy.deepcopy(annot)
                for i in annot:
                    if pd.isna(i) or i is None:
                        _ = updated_annot.pop(i, None)
                    if isinstance(annot[i], list):
                        updated_annot[i] = [y for y in annot[i] if not pd.isna(y)]
                        updated_annot[i] = [y for y in updated_annot[i] if y is not None]
                    elif pd.isna(annot[i]) or i is None: 
                        _ = updated_annot.pop(i, None)
                m.annotation = updated_annot
            #update the reaction annotations
            reac_iterator = tqdm(self.model.reactions, desc='Updating the reaction annotations') if use_progressbar else self.model.reactions
            for r in reac_iterator:
                r.annotation = utils.merge_annot_dicts(r.annotation, self._find_reaction_xref(r.annotation))
            #BUG: remove nan or None
            for r in self.model.reactions:
                annot = r.annotation
                updated_annot = copy.deepcopy(annot)
                for i in annot:
                    if pd.isna(i) or i is None:
                        _ = updated_annot.pop(i, None)
                    if isinstance(annot[i], list):
                        updated_annot[i] = [y for y in annot[i] if not pd.isna(y)]
                        updated_annot[i] = [y for y in updated_annot[i] if y is not None]
                    elif pd.isna(annot[i]) or i is None: 
                        _ = updated_annot.pop(i, None)
                r.annotation = updated_annot
            #update the gene annotations
            #TODO
            self.model.annotation['biopathopt_enriched'] = True


    def save_model(self, file_path: str) -> None:
        """
        Save a COBRApy model to a gzip-compressed JSON file.

        Args:
            file_path (str): Path to output .json.gz file.

        Returns:
            None
        """
        if file_path.endswith('json.gz') or file_path.endswith('json.gzip'):
            json_str = cobra.io.to_json(self.model)  # serialize model to JSON string
            with gzip.open(file_path, 'wt', encoding='utf-8') as f:
                f.write(json_str)    
        elif file_path.endswith('json'):
            cobra.io.save_json_model(self.model, file_path)
        elif file_path.endswith('sbml') or file_path.endswith('.xml'):
            cobra.io.write_sbml_model(self.model, file_path)
        else:
            raise TypeError('Cannot recognize the output file type')


    #################### Parse the model ##################
    # These functions are to be run when reading a model
    # to make sure that models are consistent when using them

    def _find_reaction_xref(self, r_annotation: Optional[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Attempts to find and enrich cross-references for a metabolite using multiple identifiers
        (InChIKey, MetaNetX, KEGG, BiGG) and data sources.

        Args:
            r_annotation (Optional[Dict[str, Any]]): Annotation dictionary of a reaction.

        Returns:
            Dict[str, Any]: A dictionary of enriched cross-references if found; otherwise empty.
        """
        def search_xref(uid: Union[str, list], data_search: Callable[[str], tuple]) -> Dict[str, Any]:
            """
            Searches a single identifier or a list of identifiers using the provided search function.

            Args:
                uid (Union[str, list]): Identifier or list of identifiers.
                data_search (Callable[[str], tuple]): Function returning (xref_dict, matched_id).

            Returns:
                Dict[str, Any]: The xref dictionary, or empty if not found.
            """
            if isinstance(uid, list):
                for i in uid:
                    xref, _ = data_search(i)
                    return xref
            elif isinstance(uid, str):
                xref, _ = data_search(uid)
                return xref
            else:
                logging.warning(f'Cannot find type of {uid}')
            return {}

        if r_annotation:
            mnxr = r_annotation.get('metanetx.reaction')
            if mnxr:
                xref = search_xref(mnxr, self.mnxr_xref)
                if xref:
                    return xref     
            keggr = r_annotation.get('kegg.reaction')
            if keggr:
                xref = search_xref(keggr, self.keggr_xref)
                if xref:
                    return xref   
            biggr = r_annotation.get('bigg.reaction')
            if not biggr:
                biggr = r_annotation.get('biggr')
            if not biggr:
                biggr = r_annotation.get('biggR')
            if biggr:
                xref = search_xref(biggr, self.biggr_xref)
                if xref:
                    return xref
        return {}

    
    def _find_metabolite_xref(self, m_annotation: Optional[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Attempts to find and enrich cross-references for a metabolite using multiple identifiers
        (InChIKey, MetaNetX, KEGG, BiGG) and data sources.

        Args:
            m_annotation (Optional[Dict[str, Any]]): Annotation dictionary of a metabolite.

        Returns:
            Dict[str, Any]: A dictionary of enriched cross-references if found; otherwise empty.
        """
        def search_xref(uid: Union[str, list], data_search: Callable[[str], tuple]) -> Dict[str, Any]:
            """
            Searches a single identifier or a list of identifiers using the provided search function.

            Args:
                uid (Union[str, list]): Identifier or list of identifiers.
                data_search (Callable[[str], tuple]): Function returning (xref_dict, matched_id).

            Returns:
                Dict[str, Any]: The xref dictionary, or empty if not found.
            """
            if isinstance(uid, list):
                for i in uid:
                    xref, _ = data_search(i)
                    return xref
            elif isinstance(uid, str):
                xref, _ = data_search(uid)
                return xref
            else:
                logging.warning(f'Cannot find type of {uid}')
            return {}

        def enrich_xref(xref_dict: Dict[str, Any]) -> Dict[str, Any]:
            """
            Extracts and enriches standard cross-reference fields from a raw xref result.

            Args:
                xref_dict (Dict[str, Any]): Raw xref data.

            Returns:
                Dict[str, Any]: Cleaned and enriched xref dictionary.
            """
            if xref_dict:
                try:
                    to_ret = xref_dict['xref']
                    for i in to_ret:
                        #TODO: remove nan in the annotations instead of doing this
                        if pd.isna(i) or str(i).lower()=='nan':
                            to_ret.pop(i, None)
                        else:
                            if isinstance(to_ret[i], list):
                                to_ret[i] = [y for y in to_ret[i] if not pd.isna(y)]
                                to_ret[i] = [y for y in to_ret[i] if str(y).lower() != 'nan']
                            elif pd.isna(to_ret[i]):
                                to_ret[i] = ''
                except KeyError:
                    return {}
                #TODO: remove any xref entries that are nan
                try:
                    if not pd.isna(xref_dict['SMILES']) and not xref_dict['SMILES'].lower()=='nan':
                        to_ret['smiles'] = xref_dict['SMILES']
                except KeyError:
                    pass
                try:
                    if not pd.isna(xref_dict['InChI']) and not xref_dict['InChI'].lower()=='nan':
                        to_ret['inchi'] = xref_dict['InChI']
                except KeyError:
                    pass
                try:
                    if not pd.isna(xref_dict['InChIKey']) and not xref_dict['InChIKey'].lower()=='nan':
                        to_ret['inchi_key'] = xref_dict['InChIKey']
                except KeyError:
                    pass
                return to_ret
            return {}

        if m_annotation:
            inchikey = m_annotation.get('inchi_key')
            if inchikey:
                xref = enrich_xref(search_xref(inchikey, self.inchikey_xref))
                if xref:
                    return xref
            mnxm = m_annotation.get('metanetx.chemical')
            if mnxm:
                xref = enrich_xref(search_xref(mnxm, self.mnxm_xref))
                if xref:
                    return xref     
            keggm = m_annotation.get('kegg.compound')
            if keggm:
                xref = enrich_xref(search_xref(keggm, self.keggm_xref))
                if xref:
                    return xref   
            biggm = m_annotation.get('bigg.metabolite')
            if biggm:
                xref = enrich_xref(search_xref(biggm, self.biggm_xref))
                if xref:
                    return xref
        return {}


    def _replace_depr_mnxm(self) -> None:
        """
        Replace deprecated MetaNetX chemical annotations for metabolites in the model.

        This function iterates through each metabolite in the model and checks if it has a 
        'metanetx.chemical' annotation. If the annotation is a list, it applies the 
        `single_depr_mnxm` function to each entry. If the annotation is a string, it applies 
        `single_depr_mnxm` to the string directly, updating the annotation with the result.

        Returns:
            None
        """
        # Iterate through metabolites in the model to find a match
        for metabolite in self.model.metabolites:
            if 'metanetx.chemical' in metabolite.annotation:
                # Check if the annotation is a list or a string
                if isinstance(metabolite.annotation['metanetx.chemical'], list):
                    # Update the annotation for each entry in the list
                    metabolite.annotation['metanetx.chemical'] = [
                        self.single_depr_mnxm(i) for i in metabolite.annotation['metanetx.chemical']
                    ]
                elif isinstance(metabolite.annotation['metanetx.chemical'], str):
                    # Update the annotation for a single string
                    metabolite.annotation['metanetx.chemical'] = self.single_depr_mnxm(metabolite.annotation['metanetx.chemical'])


    def _replace_dpr_mnxr(self) -> None:
        """
        Replace deprecated MetaNetX reaction annotations for reactions in the model.

        This function iterates through each reaction in the model and checks if it has a 
        'metanetx.reaction' annotation. If the annotation is a list, it applies the 
        `single_depr_mnxr` function to each entry. If the annotation is a string, it applies 
        `single_depr_mnxr` to the string directly, updating the annotation with the result.

        Returns:
            None
        """
        for reaction in self.model.reactions:
            if 'metanetx.reaction' in reaction.annotation:
                # Check if the annotation is a list or a string
                if isinstance(reaction.annotation['metanetx.reaction'], list):
                    # Update the annotation for each entry in the list
                    reaction.annotation['metanetx.reaction'] = [
                        self.single_depr_mnxr(i) for i in reaction.annotation['metanetx.reaction']
                    ]
                elif isinstance(reaction.annotation['metanetx.reaction'], str):
                    # Update the annotation for a single string
                    reaction.annotation['metanetx.reaction'] = self.single_depr_mnxr(reaction.annotation['metanetx.reaction'])

    #TODO: add a method that loops through all the ec-code of the model and 
    #finds the one that uses the substrates that is described in the model (or the closest one)
    

    ## Query the model ##

    def downstream_reactions(
        self,
        reaction_id: str, 
        ignore_metabolites: List[str] = [], 
        ignore_cofactors: bool = True,
        must_contain_met_ids: List[str] = [],
    ) -> List['Reaction']:
        """Get the downstream reactions of a reaction based on the provided reaction ID.

        Args:
            reaction_id (str): The ID of the reaction for which downstream reactions are to be retrieved.
            ignore_metabolites (List[str], optional): A list of metabolite IDs to ignore. Defaults to an empty list.

        Returns:
            List[Reaction]: A list of downstream reaction objects linked to the reaction's products.

        Raises:
            KeyError: If the given reaction ID is not found in the model.
        """
        
        # Try to fetch the reaction object by the provided ID from the model.
        try:
            r = self.model.reactions.get_by_id(reaction_id)
        except KeyError:
            raise KeyError(f'Cannot recognize the input ID: {reaction_id}')
        
        downstream_reacts = []
        if r.reversibility:
            mols = r.products + r.reactants
        else:
            mols = r.products
            
        # Loop through the products of the reaction if the model is irreversible.
        for react_product in mols:
            # Skip the metabolite if it's in the ignore list.
            if react_product.id in ignore_metabolites:
                continue
            # Check for MetaNetX annotations and ensure they are not in the cofactor list.
            if 'metanetx.chemical' in react_product.annotation:
                if ignore_cofactors:
                    # If the annotation is a list, check if any of them are in the cofactors list.
                    if isinstance(react_product.annotation['metanetx.chemical'], list):
                        if set(react_product.annotation['metanetx.chemical']) & set(self.mnxm_cofactors):
                            continue
                    # If the annotation is a string, check if it is in the cofactors list.
                    elif isinstance(react_product.annotation['metanetx.chemical'], str):
                        if react_product.annotation['metanetx.chemical'] in self.mnxm_cofactors:
                            continue
            # For each reaction linked to the product, check if it is downstream.
            for product_reaction in react_product.reactions:
                # if the reaction is reversible use both direction
                if product_reaction.reversibility:
                    # Check if the product is a reactant in the downstream reaction.
                    r_mols = [i.id for i in product_reaction.reactants + product_reaction.products]
                else:
                    r_mols = [i.id for i in product_reaction.reactants]
                if product_reaction.id != r.id:
                    if react_product.id in r_mols:
                        if must_contain_met_ids:
                            if set(must_contain_met_ids) & set(r_mols):
                                downstream_reacts.append(product_reaction)
                        else:
                            downstream_reacts.append(product_reaction)

        return downstream_reacts

    def does_mnxm_exist(
        self,
        mnxm: str,
        compartment_id: str = None
    ) -> Union[Metabolite, str]:
        """Checks if a MetaNetX metabolite exists in the COBRApy model.

        This function searches for a metabolite in the provided COBRApy model based
        on the MetaNetX ID. It also considers deprecated MetaNetX IDs if a
        deprecated ID DataFrame is provided.

        Args:
            model (Model): The COBRApy model containing metabolites.
            mnxm (str): The MetaNetX metabolite ID to search for.

        Returns:
            Union[Metabolite, str]: The found COBRApy `Metabolite` object, or an
            empty string if the metabolite is not found.
        """

        mnxm_newer = self.single_depr_mnxm(mnxm)
        if mnxm_newer!=mnxm:
            logging.warning(f"The input mnxm ({mnxm}) is depreacted, using {mnxm_newer}")
        mnxm = mnxm_newer

        # Iterate through metabolites in the model to find a match
        for metabolite in self.model.metabolites:
            # If a compartment ID is provided, restrict the search to that compartment
            if compartment_id:
                if metabolite.compartment == compartment_id:
                    if 'metanetx.chemical' in metabolite.annotation:
                        if isinstance(metabolite.annotation['metanetx.chemical'], list):
                            if mnxm in metabolite.annotation['metanetx.chemical']:
                                return metabolite
                        elif isinstance(metabolite.annotation['metanetx.chemical'], str):
                            if metabolite.annotation['metanetx.chemical'] == mnxm:
                                return metabolite
            else:
                if 'metanetx.chemical' in metabolite.annotation:
                    if isinstance(metabolite.annotation['metanetx.chemical'], list):
                        if mnxm in metabolite.annotation['metanetx.chemical']:
                            return metabolite
                    elif isinstance(metabolite.annotation['metanetx.chemical'], str):
                        if metabolite.annotation['metanetx.chemical'] == mnxm:
                            return metabolite
        return None

    def _search_inchikey(
        self,
        inchikey: str,
        inchikey_levels: int = 3,
        compartment_id: str = None
    ) -> Union[Metabolite, str]:
        """
        Searches for a metabolite in the model by its InChIKey. The search can be restricted 
        to a specific compartment if provided, and it compares up to the specified number of 
        levels of the InChIKey.

        Args:
            inchikey (str): The InChIKey of the metabolite to search for.
            inchikey_levels (int, optional): The number of InChIKey segments (separated by '-') to compare (default is 3).
            compartment_id (str, optional): The compartment ID to restrict the search to (default is None, meaning all compartments).

        Returns:
            Union[Metabolite, str]: Returns the matching COBRApy Metabolite object if found, or an empty string if no match is found.
        """
        # Iterate through metabolites in the model to find a match
        for metabolite in self.model.metabolites:
            # If a compartment ID is provided, restrict the search to that compartment
            if compartment_id:
                if metabolite.compartment == compartment_id:
                    if 'inchi_key' in metabolite.annotation:
                        if '-'.join(str(metabolite.annotation['inchi_key']).split('-')[:inchikey_levels]) == '-'.join(str(inchikey).split('-')[:inchikey_levels]):
                            return metabolite
            else:
                if 'inchi_key' in metabolite.annotation:
                    if '-'.join(str(metabolite.annotation['inchi_key']).split('-')[:inchikey_levels]) == '-'.join(str(inchikey).split('-')[:inchikey_levels]):
                        return metabolite

        # Return an empty string if the metabolite is not found
        return ''

    def does_inchikey_exist(
            self,
            inchikey: str,
            compartment_id: str = None,
            inchikey_levels: int = 3,
    ) -> Union[Metabolite, str]:
        """
        Checks if a metabolite with the given InChIKey exists in the model. It first performs a direct search 
        using the InChIKey up to two levels. If no match is found, it tries to map the InChIKey to a MNXM identifier 
        and checks if a metabolite with that MNXM identifier exists.

        Args:
            inchikey (str): The InChIKey of the metabolite to check for.
            compartment_id (str, optional): The compartment ID to restrict the search to (default is None, meaning all compartments).

        Returns:
            Union[Metabolite, str]: Returns the COBRApy Metabolite object if found, or an empty string if no match is found.
        """
        # Perform a direct search for the InChIKey in the model
        direct_search = self._search_inchikey(
            inchikey=inchikey,
            inchikey_levels=inchikey_levels,
            compartment_id=compartment_id,
        )
        if direct_search:
            return direct_search

        # if that fails, try to use use the inchikey to get the MetaNetX ID
        mnxm = None
        if inchikey_levels==3 or inchikey_levels==2:
            try:
                mnxm = self.inchikey_mnxm[inchikey]
            except KeyError:
                pass
        if inchikey_levels==2 and not mnxm:
            inchikey2 = '-'.join(inchikey.split('-')[:2])
            try:
                mnxm = self.inchikey2_mnxm[inchikey2]
            except KeyError:
                pass
        if mnxm: 
            # Check if the MNXM identifier exists in the model
            return self.does_mnxm_exist(mnxm, compartment_id=compartment_id)
        return ''

    def does_mnxr_exist(
        self,
        mnxr: str,
    ) -> Union[Reaction, str]:
        """Checks if a MetaNetX reaction exists in the COBRApy model.

        This function searches for a reaction in the provided COBRApy model based
        on the MetaNetX ID. It also considers deprecated MetaNetX IDs if a
        deprecated ID DataFrame is provided.

        Args:
            model (Model): The COBRApy model containing reactions.
            mnxr (str): The MetaNetX reaction ID to search for.
            compartment_id (str, optional): The compartment ID to restrict the search to. Defaults to None.

        Returns:
            Union[Reaction, str]: The found COBRApy `Reaction` object, or an
            empty string if the reaction is not found.
        """
        mnxr_newer = self.single_depr_mnxr(mnxr)
        if mnxr_newer!=mnxr:
            logging.warning(f"The input mnxr ({mnxr}) is depreacted, using {mnxr_newer}")
        mnxr = mnxr_newer
        # Iterate through reactions in the model to find a match
        for reaction in self.model.reactions:
            if 'metanetx.reaction' in reaction.annotation:
                # Check if the annotation is a list or a string
                # TODO: update the mnxr once you load a model
                if isinstance(reaction.annotation['metanetx.reaction'], list):
                    if mnxr in reaction.annotation['metanetx.reaction']:
                        return reaction
                elif isinstance(reaction.annotation['metanetx.reaction'], str):
                    if mnxr == reaction.annotation['metanetx.reaction']:
                        return reaction

        # Return an empty string if the reaction is not found
        return ''

    def find_mnxr_from_inchikey_reaction(
        self,
        reactants: List[Tuple[str, str]], 
        products: List[Tuple[str, str]], 
        inchikey_levels: int = 3,
        directional: bool = True,
    ) -> Optional[str]:
        """Finds the MetaNetX reaction ID (MNXR) from InChIKey-based reactants and products.

        This function constructs all possible reaction strings from the given reactants and products
        with InChIKeys, then matches these strings against a DataFrame of reaction properties to 
        find the corresponding MetaNetX reaction ID (MNXR).

        Args:
            reactants (List[Tuple[str, str]]): A list of tuples where each tuple contains 
                the coefficient (str) and the InChIKey (str) for a reactant.
            products (List[Tuple[str, str]]): A list of tuples where each tuple contains 
                the coefficient (str) and the InChIKey (str) for a product.
            inchikey_levels (int, optional): The number of levels of the InChIKey to include 
                in the reaction string. Must be 2 or 3. Defaults to 3.

        Returns:
            Optional[str]: The MetaNetX reaction ID (MNXR) if found, otherwise None.

        Example:
            >>> reactants = [("2", "ABCDEF-ABCDEF-ABCDEF"), ("1", "GHIJKL-GHIJKL-GHIJKL")]
            >>> products = [("1", "MNOPQR-MNOPQR-MNOPQR")]
            >>> find_mnxr_from_inchikey_reaction(reactants, products,  inchikey_levels=3)
            'MNXR1'
        """
        all_equation_strings = []

        # Generate all possible permutations of reactants and products and construct reaction strings
        for i in list(itertools.permutations(reactants)):
            for y in list(itertools.permutations(products)):
                all_equation_strings.append(
                    self.construct_reaction_string(
                        i,
                        y,
                        inchikey_levels=inchikey_levels
                    )
                )
                if not directional:
                    # because the reactants and products in mnxr are not in the same order sometimes
                    all_equation_strings.append(
                        self.construct_reaction_string(
                            y,
                            i,
                            inchikey_levels=inchikey_levels
                        )
                    )
        
        #search for the inchi reaction
        if inchikey_levels==3:
            for i in self.reac_prop:
                if self.reac_prop[i]['inchikey_equation'] in all_equation_strings:
                    return i
        elif inchikey_levels==2:
            for i in self.reac_prop:
                if self.reac_prop[i]['inchikey2_equation'] in all_equation_strings:
                    return i
        else:
            raise ValueError('We only handle InChIKey levels of 2 and 3')

        # Return None if no matches are found
        return None

    ##### Model Manipulation

    def remove_reaction_and_metabolites(self, reaction_id: str) -> None:
        """
        Removes a reaction and its associated metabolites from the COBRApy model.
        If any metabolites are only involved in the removed reaction, they are also removed.
        
        Args:
            model (Model): The COBRApy model object.
            reaction_id (str): The ID of the reaction to be removed.

        Returns:
            None
        """
        # Access the reaction by its ID
        if reaction_id not in self.model.reactions:
            logging.warning(f"Reaction '{reaction_id}' not found in the model.")
            return None
        reaction = self.model.reactions.get_by_id(reaction_id)
        # Get the list of metabolites associated with the reaction
        associated_metabolites = set(reaction.metabolites.keys()) 
        # Remove the reaction from the model
        self.model.remove_reactions([reaction_id])
        # Check if any of the metabolites are no longer used in any other reactions
        for metabolite in associated_metabolites:
            # If the metabolite is not part of any other reaction, remove it
            if len(metabolite.reactions) == 0:
                self.model.metabolites.remove(metabolite)
 
    ##### Model Creation #####

    ## Metabolite ##

    def create_metabolite_inchikey(
            self,
            inchikey: str,
            compartment_id: str = 'c',
            annotation: dict = {},
            inchikey_levels: int = 3,
            uid: str = None,
    ) -> Metabolite:
        """
        Creates a COBRApy Metabolite object using the provided InChIKey, chemical properties,
        and cross-reference data. If the InChIKey is found in the provided mapping (`inchikey2_mnxm`), 
        the corresponding metabolite is created; otherwise, a new Metabolite object is generated.

        Args:
            inchikey (str): The full InChIKey of the metabolite.
            compartment_id (str, optional): The compartment identifier (default is 'c').
            annotation (dict, optional): Annotations to be added to the metabolite object. Defaults to an empty dictionary.

        Returns:
            Metabolite: The created COBRApy Metabolite object.
        """
        met = self.does_inchikey_exist(inchikey, compartment_id=compartment_id, inchikey_levels=inchikey_levels)
        if met:
            return met
        ik = '-'.join(inchikey.split('-')[:inchikey_levels])
        if inchikey_levels==2:
            try:
                mnxm = self.inchikey2_mnxm[ik]
                return self.create_metabolite(
                    mnxm,
                    compartment_id,
                )
            except KeyError:
                pass
        elif inchikey_levels==3:
            try:
                mnxm = self.inchikey_mnxm[ik]
                return self.create_metabolite(
                    mnxm,
                    compartment_id,
                )
            except KeyError:
                pass
        else:
            raise KeyError('Can only handle inchikey_levels of 2 or 3')
        logging.info('Creating your own metabolite')
        # Create the COBRApy metabolite object
        if uid:
            met = Metabolite(
                id=uid,
                formula=annotation.get('formula', ''),
                name=annotation.get('name', inchikey),
                compartment=compartment_id,
            )
        else:
            met = Metabolite(
                id=annotation.get('name', 'inchikey'),
                formula=annotation.get('formula', ''),
                name=annotation.get('name', inchikey),
                compartment=compartment_id,
            )
        a = {}
        if 'xref' in annotation:
            a = annotation['xref']
        a['inchikey'] = inchikey
        if 'smiles' in annotation:
            a['smiles'] = annotation['smiles']
        if 'inchi' in annotation:
            a['inchi'] = annotation['inchi']
        # Add annotations to the metabolite object
        met.annotation = a
        return met

    def create_metabolite(
        self,
        mnxm: str,
        compartment_id: str = 'c'
    ) -> Metabolite:
        """Creates a COBRApy metabolite object from MetaNetX data.

        This function generates a COBRApy `Metabolite` object using MetaNetX data from
        the provided dataframes. It converts MetaNetX annotations to BiGG or other
        relevant formats, and it populates the metabolite with relevant information
        such as formula, name, and compartment.

        Args:
            mnxm (str): The MetaNetX metabolite ID.
            compartment_id (str, optional): Compartment ID where the metabolite resides. Defaults to 'c'.

        Returns:
            Metabolite: A COBRApy `Metabolite` object populated with the specified data.
        """

        chem_prop, fetched_mnxm  = self.mnxm_xref(mnxm)
        if not fetched_mnxm==mnxm:
            raise KeyError(f'You are using a deprecated {mnxm}')
        chem_xref = chem_prop['xref']
        # Convert MetaNetX annotations to BiGG or other annotations
        annot_dict: Dict[str, list] = {}
        for prefix in chem_xref:
            if prefix in self.mnxm_bigg_annot_convert:
                key = self.mnxm_bigg_annot_convert[prefix]
                if key not in annot_dict:
                    annot_dict[key] = chem_xref[prefix]
                else:
                    annot_dict[key] += chem_xref[prefix]
            else:
                if prefix not in annot_dict:
                    annot_dict[prefix] = chem_xref[prefix]
                else:
                    annot_dict[prefix] += chem_xref[prefix]

        # Ensure all values in annotation dictionary are unique
        for key in annot_dict:
            annot_dict[key] = list(np.unique(annot_dict[key]))

        # Add InChIKey to the annotation dictionary
        annot_dict['inchi_key'] = chem_prop['InChIKey']

        # Create a unique identifier for the metabolite
        if 'bigg.metabolite' in annot_dict:
            # Use the smallest BiGG name and append the compartment ID
            uid = str(min(annot_dict['bigg.metabolite'], key=len)) + '_' + compartment_id
        else:
            uid = mnxm

        # Create the COBRApy metabolite object
        met = Metabolite(
            id=uid,
            formula=chem_prop['formula'],
            name=chem_prop['name'],
            compartment=compartment_id,
        )

        # Add annotations to the metabolite object
        met.annotation = annot_dict

        return met


    ## Reactions ##

    def create_reaction(
        self,
        mnxr: str,
        name: str = None,
        reverse_reaction: bool = False,
        lower_bound: float = 0.0,
        upper_bound: float = 1000.0,
        compartment_id: str = 'c',
    ) -> Reaction:
        """Creates a COBRApy reaction object from MetaNetX data.

        This function generates a COBRApy `Reaction` object using MetaNetX data from 
        the provided dataframes, with optional reverse reaction direction, lower 
        and upper bounds, and specified compartment ID.

        Args:
            mnxr (str): The MetaNetX reaction ID.
            reverse_reaction (bool, optional): Flag to reverse the reaction direction. Defaults to False.
            lower_bound (float, optional): Lower bound for the reaction flux. Defaults to 0.0.
            upper_bound (float, optional): Upper bound for the reaction flux. Defaults to 1000.0.
            compartment_id (str, optional): Compartment ID where the reaction takes place. Defaults to 'c'.

        Returns:
            Reaction: A COBRApy `Reaction` object populated with the specified data.
        """
    
        # Create a temporary dictionary to store the annotation mapping
        mnxr_prop, fetched_mnxr = self.mnxr_xref(mnxr)
        if not mnxr==fetched_mnxr:
            raise KeyError(f'The fetched MNXR is a deprecated one ({mnxr})')
        mnxr_xref = mnxr_prop['xref']

        # Convert MetaNetX annotations to BiGG annotations
        annot_dict: Dict[str, List[str]] = {}
        for prefix in mnxr_xref:
            if prefix in self.mnxr_bigg_annot_convert:
                key = self.mnxr_bigg_annot_convert[prefix]
                annot_dict[key] = mnxr_xref[prefix]
            else:
                annot_dict[prefix] = mnxr_xref[prefix]

        # Ensure all values in annotation dictionary are unique
        for key in annot_dict:
            annot_dict[key] = list(np.unique(annot_dict[key]))

        # Create a name for the reaction
        if 'bigg.reaction' in annot_dict:
            # Use the smallest BiGG name and add the compartment ID
            uid = str(min(annot_dict['bigg.reaction'], key=len))
        else:
            uid = mnxr

        # Parse the reaction equation to get reactants and products
        reac_str = mnxr_prop['mnx_equation']
        reactants, products = utils.parse_mnxr_equation(reac_str)

        # Create a stoichiometry dictionary
        tmp_stoichio: Dict[str, float] = {}
        for item in reactants:
            coefficient = -float(item[0]) if reverse_reaction else float(item[0])
            tmp_stoichio[item[1]] = coefficient

        for item in products:
            coefficient = float(item[0]) if reverse_reaction else -float(item[0])
            tmp_stoichio[item[1]] = coefficient

        # Replace MetaNetX IDs with model-specific IDs
        stoichio: Dict[str, float] = {}
        for mnxm_id in tmp_stoichio:
            met = self.does_mnxm_exist(mnxm_id)
            if met:
                stoichio[met] = tmp_stoichio[mnxm_id]
            else:
                met = self.create_metabolite(mnxm_id, compartment_id=compartment_id)
                stoichio[met] = tmp_stoichio[mnxm_id]

        # Create and configure the COBRApy reaction object
        reaction = Reaction(uid)
        if name:
            reaction.name = name
        else:
            reaction.name = uid
        reaction.subsystem = None
        reaction.lower_bound = lower_bound
        reaction.upper_bound = upper_bound
        reaction.add_metabolites(stoichio)
        reaction.annotation = annot_dict    

        return reaction

    def create_inchikey_reaction(
        self,
        reactants: List[Tuple[float, str]],
        products: List[Tuple[float, str]],
        uid: str,
        name: str = None,
        inchikey_levels: int = 3,
        lower_bound: float = 0.0,
        upper_bound: float = 1000.0,
        compartment_id: str = 'c',
        annotations: Dict[str, Dict] = {},
    ) -> Reaction:
        """
        Create a reaction using InChIKey identifiers for reactants and products.

        This function builds a stoichiometric reaction from given reactants and products. 
        The stoichiometry is constructed using InChIKey identifiers, and the metabolites 
        are either fetched if they exist or created if they don't.

        Args:
            reactants (List[Tuple[float, str]]): A list of tuples representing reactants, 
                where each tuple contains the stoichiometric coefficient (float) and the 
                InChIKey (str) of the reactant.
            products (List[Tuple[float, str]]): A list of tuples representing products, 
                where each tuple contains the stoichiometric coefficient (float) and the 
                InChIKey (str) of the product.
            uid (str): Unique identifier for the reaction.
            lower_bound (float, optional): Lower bound for the reaction flux. Defaults to 0.0.
            upper_bound (float, optional): Upper bound for the reaction flux. Defaults to 1000.0.
            compartment_id (str, optional): Compartment identifier for the reaction. Defaults to 'c'.
            annotations (Dict[str, Dict], optional): A dictionary of annotations for metabolites, 
                where the keys are InChIKeys and values are annotation dictionaries. Defaults to {}.

        Returns:
            Reaction: A COBRApy `Reaction` object with the specified reactants and products.

        Raises:
            KeyError: If a metabolite cannot be added because it lacks an annotation.

        """
        # Create a temporary dictionary to hold the reaction stoichiometry based on InChIKeys
        tmp_stoichio = {}
        # Add reactants to the stoichiometry with negative coefficients
        for r in reactants:
            tmp_stoichio[r[1]] = -float(r[0])
        # Add products to the stoichiometry with positive coefficients
        for p in products:
            tmp_stoichio[p[1]] = float(p[0])

        # Initialize the final stoichiometry dictionary
        stoichio = {}
        # Iterate through the temporary stoichiometry to check if the metabolite exists or needs to be created
        for i in tmp_stoichio:
            # If metabolite doesn't exist, attempt to create it
            try:
                stoichio[self.create_metabolite_inchikey(
                    i, 
                    compartment_id,
                    annotation=annotations[i],
                    inchikey_levels=inchikey_levels,
                )] = tmp_stoichio[i]
            except KeyError:
                # Raise error if the metabolite cannot be added due to missing annotation
                raise KeyError(f'Cannot add {i}, {tmp_stoichio} ---- {annotations}')

        # Create the reaction object and set its properties
        reaction = Reaction(uid)
        if name:
            reaction.name = name
        else:
            reaction.name = uid
        reaction.subsystem = None
        reaction.lower_bound = lower_bound
        reaction.upper_bound = upper_bound
        # Add the metabolites and their stoichiometry to the reaction
        reaction.add_metabolites(stoichio)
        # Set the reaction annotations to an empty dictionary
        reaction.annotation = {}

        return reaction

    def find_gene_from_uniprot(self, uniprot_id: str) -> str:
        """
        Finds the gene in the model that corresponds to a given UniProt ID.

        Args:
            uniprot_id (str): The UniProt ID to search for in the model's genes.

        Returns:
            str: The model gene ID corresponding to the given UniProt ID.

        Raises:
            TypeError: If no genes are found or if multiple genes match the UniProt ID.
        """
        model_gene_id = []
        # Search for genes in the model with the given UniProt ID
        for g in self.model.genes:
            if 'uniprot' in g.annotation and g.annotation['uniprot'] == uniprot_id:
                model_gene_id.append(g)
        # Handle cases where the gene is not found or multiple matches exist
        if not model_gene_id:
            raise TypeError(f"Could not find UniProt ID: {uniprot_id}")
        elif len(model_gene_id) > 1:
            raise TypeError(f"Found multiple genes for UniProt ID {uniprot_id}: {model_gene_id}")
        return model_gene_id[0]

    def find_genes_from_uniprots(self, uniprot_ids: List[str]) -> List[str]:
        """
        Finds the genes in the model that correspond to a list of UniProt IDs.

        Args:
            uniprot_ids (List[str]): A list of UniProt IDs to search for in the model's genes.

        Returns:
            List[str]: A list of model gene IDs corresponding to the given UniProt IDs.

        Raises:
            TypeError: If `find_gene_from_uniprot` raises an error for any UniProt ID.
        """
        model_gene_ids: List[str] = []
        for u_i in uniprot_ids:
            m_g = self.find_gene_from_uniprot(u_i)
            if m_g not in model_gene_ids:
                model_gene_ids.append(m_g)
        return model_gene_ids    

    def knockout_uniprot(self, uniprot_ids: List[str]) -> List[str]:
        """
        Knocks out reactions associated with given UniProt IDs in the model.

        This function identifies genes in the model by their UniProt annotations,
        validates their uniqueness, and knocks out reactions associated with those genes.

        Args:
            uniprot_ids (List[str]): A list of UniProt IDs for which to knock out reactions.

        Returns:
            List[str]: A list of IDs of reactions that were knocked out.

        Raises:
            TypeError: If a UniProt ID cannot be found in the model or if multiple genes are associated with the same UniProt ID.

        Example:
            ```python
            model = SomeModel()
            knocked_out_reactions = model.knockout_uniprot(["P12345", "Q67890"])
            print(knocked_out_reactions)
            ```
        """
        knocked_out_reacs: List[str] = []
        model_gene_ids = self.find_genes_from_uniprots(uniprot_ids)

        #get the unique reactions
        unique_reactions = []
        for i in model_gene_ids:
            for y in i.reactions:
                if y not in unique_reactions:
                    unique_reactions.append(y)

        #identify the genes that are knocked out
        model_gene_ids_str = [i.id for i in model_gene_ids]
        for reaction in unique_reactions:
            gene_lists = utils.extract_genes_from_gpr(reaction.gpr.to_string())
            # only knockout the reaction if gene is a complex and has no replacement
            rm_genes = [ gene_list for gene_list in gene_lists
                if not any(element in model_gene_ids_str for element in gene_list)
            ]
            if not rm_genes:
                self.model.reactions.get_by_id(reaction.id).knock_out()
                knocked_out_reacs.append(reaction.id)

        return knocked_out_reacs
