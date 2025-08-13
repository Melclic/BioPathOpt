import re
import copy
import time
import os
import pandas as pd
import json
import xmltodict
from bioservices import UniProt
import cobra
from cobra import Model, Reaction
from cobra.util.solver import set_objective
from cobra.flux_analysis import pfba
import requests
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import networkx as nx
from io import StringIO
from collections import defaultdict
from ete3 import NCBITaxa
import numpy as np
import logging
from tqdm import tqdm
from rapidfuzz.distance import Levenshtein
import gzip

from typing import List, Union, Dict, Optional, Tuple, Any

from biopathopt import ModelBuilder
from biopathopt import dlkcat
from biopathopt.utils import inchikey_layer_extract, annotation_contains_value

"""
This is a collection of functions to build and enzyme constrained model
Includes methods to fetch kcat from MIRIAM annotation and sMOMENT definition
and other

TODO: change such that all use taxonomy instead of species name
"""

import warnings
warnings.filterwarnings(
    "ignore",
    category=SyntaxWarning,
    message=".*Malformed gene_reaction_rule.*"
)

class EnzymeConstrainedModel(ModelBuilder):
    
    def __init__(self, path_to_model: str, species_name=None, taxonomy_id=None, use_progressbar: bool = False):
        super().__init__(path_to_model=path_to_model, use_progressbar=use_progressbar)
        #TODO: print the stats of uniprot, structure, etc... so that we can
        # give a warning if the coverage is too low
        self.species_name = species_name
        self.taxonomy_id = taxonomy_id
        self.ec_model = None
        if not self.taxonomy_id:
            if self.model:
                self.taxonomy_id = self.model.annotation.get('taxonomy', None)
        if species_name and not self.taxonomy_id:
            self.taxonomy = self._get_taxid_from_species(species_name)
        if self.taxonomy_id and not species_name:
            self.species_name = self._get_species_name(self.taxonomy_id)
        if self.model:
            if not self.model.annotation.get('taxonomy', None) and self.taxonomy_id:
                self.model.annotation['taxonomy'] = self.taxonomy_id
        #TODO: Add the taxonomy id and species name in the model annotation
        self.brenda_ec_G = self._json_brenda_to_G()
        self.cofactors_inchikey_layers3 = []
        self.cofactors_inchikey_layers2 = []
        path_file = os.path.join(
            self.base_dir, "flatfiles/EC_kcat_max.json"
        )
        self.brenda_kcat_max = json.load(open(path_file))
        for mnxm in self.mnxm_cofactors:
            try:
                inchikey = self.mnxm_inchikey[mnxm]
                self.cofactors_inchikey_layers3.append(inchikey)
                self.cofactors_inchikey_layers2.append('-'.join(inchikey.split('-')[:-1]))
            except KeyError:
                pass
        self.mnxr_cofactor_transport = ['MNXR146072', 'MNXR02', 'MNXR99522', 'MNXR101670', 'MNXR96810', 'MNXR99505', 'MNXR100495', 'MNXR105280', 'MNXR101507', 'MNXR96436', 'MNXR104460', 'MNXR101958', 'MNXR191149', 'MNXR96951', 'MNXR100625', 'MNXR96499', 'MNXR98640', 'MNXR98641', 'MNXR101912', 'MNXR101950', 'MNXR104469', 'MNXR100950', 'MNXR101804', 'MNXR102090', 'MNXR191153', 'MNXR142917', 'MNXR96954', 'MNXR104456', 'MNXR96797', 'MNXR145749', 'MNXR155386', 'MNXR96821']


    #### get the ecModel

    def return_ec_model(
            self,
            substrate_exchange_reaction_id: str,
            substrate_concentration: float = 10.0,
            enzyme_mass_fraction: float = 0.405, 
            protein_fraction: float = 0.56, 
            enzyme_saturation: float = 1.0,
            lowerbound: int = None,
            upperbound: int = None,
            overwrite_kcat: dict  = {}):
        ec_model = copy.deepcopy(self.model)
        try:
            substrate_exchange_reaction = ec_model.reactions.get_by_id(substrate_exchange_reaction_id)
        except KeyError:
            logging.error(f'Cannot recognise input substrate reaction id {substrate_exchange_reaction_id}')
            return None
        #save the original reaction ids
        for r in ec_model.reactions:
            r.annotation['original_id'] = r.id
        #NOTE the order of the two function matters. 
        """
        Option 1: Run the isoenzyme split first and then make the whole system 
        irreversible. This converges to a solution faster when optimizing... 
        not sure why
        Option 2: (default in original ECMpy) Turn reactions irreversible first
        and then generate isoenzyme split. This does not find 0.66 after 50 iterations
        when optimizing 
        """
        self._isoenzyme_split(ec_model)
        self._convert_to_irreversible(ec_model, substrate_exchange_reaction)
        #self._isoenzyme_split(ec_model)
        #set model level paramters
        if not lowerbound:
            lowerbound = 0
        if not upperbound:
            upperbound = round(protein_fraction * enzyme_mass_fraction * enzyme_saturation, 3)
        ec_model.annotation['enzyme_constrain'] = {
            'enzyme_mass_fraction': enzyme_mass_fraction, 
            'total_protein_fraction': protein_fraction,
            'average_saturation': enzyme_saturation, 
            'lowerbound': lowerbound, 
            'upperbound': upperbound,
        }
        #set model coefficients
        coefficients = {}
        for r in ec_model.reactions:
            try:
                if r.id in overwrite_kcat:
                    logging.debug(f'Setting reaction {r.id} kcat from {r.annotation["kcat"]} to {overwrite_kcat[r.id]}')
                    r.annotation['kcat'] = overwrite_kcat[r.id]
                    r.annotation['kcat_mw'] = overwrite_kcat[r.id] * 3600*1000/r.annotation.get('mw', 0.0)
                coefficients[r.forward_variable] = 1.0 / float(r.annotation.get('kcat_mw', 0.0))
                #logging.debug(f'Updated the coefficient for {r.id} - {r.forward_variable}: {coefficients[r.forward_variable]}')
            except (ZeroDivisionError, TypeError) as e:
                pass
        constraint = ec_model.problem.Constraint(0, lb=lowerbound, ub=upperbound)
        ec_model.add_cons_vars(constraint)
        ec_model.solver.update()
        constraint.set_linear_coefficients(coefficients=coefficients)
        self.ec_model = ec_model
        return ec_model


    ### ec model optimize
    def optimize_ec_model(
        self,
        target_growth_rate_h: float,
        #substrate_id: str,
        carbon_exchange_reaction: str,
        model_objective: str = None,
        substrate_concentration: float = 10.0,
        enzyme_mass_fraction: float = 0.405, 
        protein_fraction: float = 0.56, 
        enzyme_saturation: float = 1.0,
        lowerbound: int = None,
        upperbound: int = None,
        maximum_optimization_rounds: int = 50,
        use_progressbar: bool = False,
    ):
        """Optimize EC model to get to a realistic growth value. This method 
        compensates for the shortcomings of this method that may restrict the 
        flux through a reaction such that the growth rate becomes unrealistic.
        Here, we replace those reactions that have the highest enzyme ratio
        with their corresponding EC number highest recorded kcat.
        """
        """
        carbon_reac = self.change_carbon_source(
            metabolite_id=substrate_id,
            extracellular_compartment_id=extracellular_compartment_id,
            substrate_concentration=substrate_concentration,)
        """
        #logging.info(f'The carbon reaction is: {carbon_reac.id}')
        if carbon_exchange_reaction not in self.model.reactions:
            raise KeyError(f'The reaction {carbon_exchange_reaction} does not exist in the model')
        ec_model = self.return_ec_model(
                substrate_exchange_reaction_id=carbon_exchange_reaction,
                substrate_concentration=substrate_concentration,
                enzyme_mass_fraction=enzyme_mass_fraction, 
                protein_fraction=protein_fraction, 
                enzyme_saturation=enzyme_saturation,
                lowerbound=lowerbound,
                upperbound=upperbound,
            )
        if model_objective:
            if model_objective not in ec_model.reactions:
                raise KeyError(f'The objective {model_objective} does not exist')
            ec_model.obj = model_objective

        pfba_sol = pfba(ec_model)
        reaction_E = {}
        #given these constrains run and calculate E
        for reac_id, flux in pfba_sol.fluxes.items():
            r = ec_model.reactions.get_by_id(reac_id)
            try:
                reaction_E[r.id] = float(flux)/float(r.annotation.get('kcat_mw', 0.0))
            except (ZeroDivisionError, TypeError) as e:
                pass
        sum_E = sum([reaction_E[i] for i in reaction_E])
        logging.debug(f'sum_E: {sum_E}')
        enz_ratio = {i: reaction_E[i]/sum_E for i in reaction_E}
        sorted_enz_ratio = sorted(enz_ratio.items(), key=lambda item: item[1], reverse=True)
        logging.debug('sorted_enz_ratio: {sorted_enz_ratio}')
        model_gr = ec_model.slim_optimize()
        logging.debug(f'-------------- {model_gr} ---------------')
        updated_reacs = []
        #while not np.isclose(model_gr, target_growth_rate_h, rtol=rtol, atol=atol):
        if use_progressbar:
            progress = tqdm(total=maximum_optimization_rounds, desc=f'Growth rate is {model_gr}')
        num_round = 0
        new_kcat = {}
        while (model_gr<target_growth_rate_h and num_round<maximum_optimization_rounds):
            #select the next reaction
            select_reaction = None
            for s_r, e_ratio in sorted_enz_ratio:
                if s_r not in updated_reacs:
                    select_reaction = s_r
                    break
            logging.debug(f'select_reaction: {select_reaction}')
            # update the kcat to max
            update_reac_kcat_mw = {}
            reac = ec_model.reactions.get_by_id(select_reaction)
            if "ec-code" in reac.annotation.keys():
                ec_number = reac.annotation["ec-code"]
                kcat_max_list = []
                if isinstance(ec_number, str):
                    ec_number = [ec_number]
                for ec in ec_number:
                    logging.debug(f'===== EC: {ec} =====')
                    logging.debug(f"mw: {reac.annotation.get('mw')}")
                    logging.debug(f"old kcat: {reac.annotation.get('kcat')}")
                    logging.debug(f"old kcat_mw: {reac.annotation.get('kcat_mw')}")
                    reaction_kcat_max = self.brenda_kcat_max.get(ec, {'kcat_max': None})['kcat_max']
                    if reaction_kcat_max:
                        kcat_max_list.append(reaction_kcat_max)
                    if len(kcat_max_list)>0:
                        reaction_kcat_max = np.max(kcat_max_list)
                    else:
                        reaction_kcat_max = 0
                logging.debug(f'reaction_kcat_max: {reaction_kcat_max}')
                logging.debug(reac.annotation.get('kcat', 99999999.0))
                if reac.annotation.get('kcat', 99999999.0) < reaction_kcat_max:
                    reac.annotation['kcat'] = reaction_kcat_max #h_1
                    logging.debug(f'new kcat: {reaction_kcat_max}')
                    try:
                        #TODO: change to just multiply the kcat instead of finding the highest recorded kcat
                        # this is erronous since these include genetically engineered enzymes
                        new_kcat[reac.id] = reaction_kcat_max
                    except ZeroDivisionError:
                        pass
            #### update the coefficients
            ec_model = self.return_ec_model(
                    substrate_exchange_reaction_id=carbon_exchange_reaction,
                    substrate_concentration=substrate_concentration,
                    enzyme_mass_fraction=enzyme_mass_fraction, 
                    protein_fraction=protein_fraction, 
                    enzyme_saturation=enzyme_saturation,
                    lowerbound=lowerbound,
                    upperbound=upperbound,
                    overwrite_kcat=new_kcat,
                )
            ### update the results
            updated_reacs.append(select_reaction)
            model_gr = ec_model.slim_optimize()
            num_round += 1
            if use_progressbar:
                progress.set_description(f"Growth Rate is {model_gr}")
                progress.update(1)
            logging.debug(f'-------------- {model_gr} ---------------')
        self.ec_model = ec_model
        return updated_reacs


    def enzyme_enrich_model(
            self, 
            species_name = None,
            taxonomy_id = None,
            use_progressbar: bool = False):
        # set the subunit and split the model to only forward reactions
        self._get_subunit_data(use_progressbar=use_progressbar)
        self._get_uniprot_aaseq_mw(use_progressbar=use_progressbar)
        self.calculate_dlkcat(use_progressbar=use_progressbar)
        self.database_kcat(species_name=species_name, taxonomy_id=taxonomy_id, use_progressbar=use_progressbar)
        self._set_kcat()



    ##### DLKcat ######


    def calculate_dlkcat(self, use_progressbar: bool = False):
        dl_kcat = dlkcat.KcatPredictor()
        #kcat_pred.model_predict_kcat(self.model)
        reac_iterator = tqdm(self.model.reactions, desc='Deep learning estimation of kcat') if use_progressbar else self.model.reactions
        for r in reac_iterator:
            #if you already have some
            if 'DLkcat' in r.annotation:
                if r.annotation['DLkcat'].get('kcat'):
                    continue
            kcat = 0
            kcat_sec_value = 0
            if r.genes and len(r.reactants) > 1:
                for m in r.reactants:
                    if 'metanetx.chemical' in m.annotation:
                        if isinstance(m.annotation['metanetx.chemical'], list):
                            is_cof = False
                            for i in m.annotation['metanetx.chemical']:
                                if i in self.mnxm_cofactors:
                                    is_cof = True
                            if is_cof:
                                continue
                        elif isinstance(m.annotation['metanetx.chemical'], str):
                            if m.annotation['metanetx.chemical'] in self.mnxm_cofactors:
                                continue
                    #if you don't have smiles use 
                    smiles = m.annotation.get('smiles')
                    if not smiles:
                        smiles = m.annotation.get('SMILES')
                    inchi = m.annotation.get('inchi')
                    if not inchi:
                        inchi = m.annotation.get('InChI')
                    for g in r.genes:
                        sequence = g.annotation.get('aaseq')
                        kcat_sec_value = dl_kcat.predict_kcat(
                                smiles=smiles, 
                                sequence=sequence, 
                                inchi=inchi)
                        #keep only the largest (i.e. slowest) reaction subunit
                        if kcat_sec_value:
                            if kcat_sec_value>kcat:
                                kcat = kcat_sec_value
                        logging.info(f"Predicted kcat for reaction {r.id} and molecule {m.name} and gene {g.name}: {kcat_sec_value}")
            if kcat==0:
                kcat = None
            r.annotation['DLkcat'] = {}
            r.annotation['DLkcat']['kcat'] = kcat
        self._calculate_mw_kcat()


    def database_kcat(
            self, 
            species_name=None,
            taxonomy_id=None,
            inchikey_levels: int = 2,
            join_mode='mean',
            closest_taxonomy=True,
            use_progressbar: bool = False):
        if species_name:
            self.species_name = species_name
        #TODO: remove this since its redundant with the constructor
        if taxonomy_id:
            self.taxonomy_id = taxonomy_id
        if not self.species_name and not self.taxonomy_id:
            taxid = self.model.annotation.get('taxonomy', None)
            if taxid:
                self.taxonomy_id = taxid
                self.species_name = self._get_species_name(taxid)
                if not self.species_name:
                    raise NameError(f'species_name cannot be found automatically, please define manually')
            else:
                raise NameError(f'species_name cannot be found automatically, please define manually')
        self._search_kcat(
                self.species_name,
                inchikey_levels=inchikey_levels,
                join_mode=join_mode,
                closest_taxonomy=closest_taxonomy,
                use_progressbar=use_progressbar,
            )
        self._calculate_mw_kcat()


    def calculate_enzyme_mass_fraction(self, gene_abundance_csv_path, abundance_colname, name_colname, model_gene_name_feature):
        '''
        Calculate the enzyme mass fraction (f) based on protein abundance.

        Arguments:
        * gene_abundance_file: str - File path of the gene abundance data file.
        * gene_abundance_colname: str - Name of the column in the gene abundance file containing the abundance values.
        * taxonom_id: str - Taxonomy ID of the organism.

        Returns:
        * f: float - The enzyme mass fraction.

        '''
        model_gene_list = [eachg.id for eachg in model.genes]
        uni_model_gene_list = list(set(model_gene_list))

        gene_abundance_df = pd.read_csv(gene_abundance_csv_path, index_col=None)
        #gene_names = gene_abundance_df[name_colname].to_list()
        gene_name_model_id = {}
        for g_name in gene_abundance_df[name_colname]:
            for g in self.model.genes:
                #check if any of the entries have the gene name in the annotation
                #WARNING does not check if there are multiple hits
                if annotation_contains_value(g.annotation, g_name):
                    gene_name_model_id[g_name] = g.id
        enzy_abundance = 0
        pro_abundance = 0
        gene_name_abundance = gene_abundance_df[[name_colname, abundance_colname]].set_index(name_colname).to_dict()[abundance_colname]
        for g_name in gene_abundance_df[name_colname]:
            abundance = None
            try:
                abundance = gene_name_abundance[g_name] * self.model.genes.get_by_id(gene_name_model_id[g_name]).annotation.get('total_subunits_mw', 0)
                pro_abundance += abundance
                enzy_abundance += abundance
            except KeyError:
                uniprot_query_url = f"https://rest.uniprot.org/uniprotkb/search?query={g_name}+AND+organism_id:{self.taxonomy_id}&format=tsv&fields=accession,mass"
                uniprot_mass = requests.get(uniprot_query_url).text.split("\n")[1].split("\t")[1]
                abundance = gene_name_abundance[g_name] * int(uniprot_mass)
                pro_abundance += abundance
        return enzy_abundance / pro_abundance

    #### metabolic engineering strategies

    def run_FSEOF(
        self,
        substrate_reaction_id: str,
        biomass_reaction_id: str,
        obj_reaction_id: str,
        carbon_exchange_reaction: str,
        substrate_reaction_flux: float = 10.0,
    ):
        #self.change_carbon_source()
        ec_model = ecm.gen_ec_model(carbon_exchange_reaction)
        ec_model.objective = biomass_id
        model_pfba_solution = cobra.flux_analysis.pfba(ec_model)
        wt_sol = model_pfba_solution.fluxes.to_dict()
        wt_biomass_sol = model_pfba_solution.fluxes[biomass_id]

        exlist = [float(round(i, 3)) for i in np.linspace(0.5 * wt_biomass_sol, wt_biomass_sol * 0.9, 10)]

        all_cond_sol = {}
        all_fc_sol = {}
        for cond in exlist:
            tmp_model = copy.deepcopy(ec_model)
            self.change_carbon_source()
            tmp_model.reactions.get_by_id(biomass_reaction_id).bounds = (cond, cond)
            tmp_model.objective = obj_reaction_id
            model_pfba_solution = cobra.flux_analysis.pfba(tmp_model)
            cond_sol = model_pfba_solution.fluxes.to_dict()
            all_cond_sol[cond] = cond_sol
            fc_sol = {i: cond_sol[i]/wt_sol[i] for i in cond_sol.keys()}
            all_fc_sol[cond] = fc_sol
        fc_mean = {i: [] for i in wt_sol}
        for i in all_fc_sol:
            for y in all_fc_sol[i]:
                if not np.isinf(all_fc_sol[i][y]):
                    fc_mean[y].append(all_fc_sol[i][y])
                else:
                    fc_mean[y].append(1000)
        fc_mean = {i: np.mean([y for y in fc_mean[i] if y]) for i in fc_mean}
        fc_status = {}
        for i in fc_mean:
            if all(np.array(fc_mean[i])<=1.0) and all(np.array(fc_mean[i])>=0.95):
                fc_status[i] = 'unchanged'
            elif all(np.array(fc_mean[i])<0.95):
                fc_status[i] = 'down'
            elif all(np.array(fc_mean[i])>1.0):
                fc_status[i] = 'up'
        fc_mean = {i: np.mean([y for y in fc_mean[i] if y]) for i in fc_mean}
        return fc_mean

    ################################
    ####### Private Functions ######
    ################################

        
    def _set_kcat(self):
        #select the kcat and kcat_MW values
        reac_no_kcat = []
        for r in self.model.reactions:
            kcat = r.annotation.get('kcat')
            if not kcat:
                for i in ['DBkcat', 'DBsa', 'DLkcat']:
                    kcat_entry = r.annotation.get(i)
                    if kcat_entry:
                        if not kcat:
                            kcat = kcat_entry.get('kcat')
                            if kcat:
                                r.annotation['kcat'] = kcat
                                break
            kcat_mw = r.annotation.get('kcat_mw')
            if not kcat_mw:
                for i in ['DBkcat', 'DBsa', 'DLkcat']:
                    kcat_entry = r.annotation.get(i)
                    if kcat_entry:
                        if not kcat_mw:
                            kcat_mw = kcat_entry.get('kcat_mw')
                            if kcat_mw:
                                r.annotation['kcat_mw'] = kcat_mw
                                break
            if not kcat or not kcat_mw:
                reac_no_kcat.append(r.id)
        logging.info(f'There are {len(reac_no_kcat)} reactions with no kcat values')

    def _calculate_mw_kcat(self):
        for r in self.model.reactions:
            totalmass = 0
            for g in r.genes:
                totalmass += g.annotation.get('total_subunits_mw', 0)
            r.annotation['mw'] = totalmass
            try:
                dlkcat = r.annotation.get('DLkcat').get('kcat')
                if dlkcat:
                    r.annotation['DLkcat']['kcat_mw'] = dlkcat * 3600 * 1000 / totalmass
            except ZeroDivisionError:
                r.annotation['DLkcat']['kcat_mw'] = None
            except AttributeError:
                pass
            try:
                dbkcat = r.annotation.get('DBkcat').get('kcat')
                if dbkcat:
                    r.annotation['DBkcat']['kcat_mw'] = dbkcat * 3600 * 1000 / totalmass
            except ZeroDivisionError:
                r.annotation['DBkcat']['kcat_mw'] = None
            except AttributeError:
                pass
            try:
                dbsa = r.annotation.get('DBsa').get('kcat')
                if dbsa:
                    r.annotation['DBsa']['kcat_mw'] = dbsa * 3600 * 1000 / totalmass
            except ZeroDivisionError:
                r.annotation['DBsa']['kcat_mw'] = None
            except AttributeError:
                pass

    ### Subunit split and protein search and Molecular Weight (MW)

    def _get_subunit_num(self, protein_names: List[str], gene_names: List[str], txt_subunit: str) -> str:
        """
        Infers the number of subunits from protein and gene names and subunit description text.
        THIS IS A MESS AND NEEDS A REWRITE AND UNITTEST

        Args:
            protein_names (List[str]): List of protein names associated with the gene.
            gene_names (List[str]): List of gene names associated with the gene.
            txt_subunit (str): Subunit description text retrieved from UniProt.

        Returns:
            str: The inferred subunit number, as a string. Values include:
                 - Exact integers (e.g., '1', '2', etc.)
                 - 'manual' if inferred heuristically
                 - 'vacuum' if the subunit text is minimal or empty
        """
        lPossibleSubunitsFromName = []
        lAllNames = protein_names + gene_names
        lAllNames = [i for i in lAllNames if not pd.isna(i)]
        lnameEnzyme = []
        possibleWords = ['subunit', 'subunits', 'component', 'alpha chain', 'beta chain', 'gamma chain',
                         '30S ribosomal protein', '50S ribosomal protein', 'binding protein', 'large chain',
                         'small chain', 'permease protein', 'insertase', 'translocase protein', 'accessory protein',
                         'UvrABC system protein', 'Chaperonin', 'Co-chaperonin', 'assembly factor', 'recombinase',
                         'flavoprotein']
        namesWithSubunitsIndication = []
        
        for p in possibleWords:
            if not lAllNames or all(x is None or pd.isna(x) for x in lAllNames):
                namesWithSubunitsIndication = []
            else:
                namesWithSubunitsIndication = [el for el in lAllNames if p in el]
            
            if re.search('RNA-binding protein|methyltransferase|Nucleotide-binding protein|Phosphotransferase|\
                          Transport|AMP-forming|Two-component|SsrA|Ribonuclease', str(namesWithSubunitsIndication)):
                continue
            
            if len(namesWithSubunitsIndication) != 0:
                for n in namesWithSubunitsIndication:
                    if n.endswith(p) or p + ',' in n:
                        try:
                            lPossibleSubunitsFromName.append(n.split(p)[0].split()[-1].strip(','))
                            lnameEnzyme.append(' '.join(n.split(p)[0].split()[:-1]).lower().strip())
                            break
                        except:
                            logging.info(n)
                    elif 'alpha-ketoacid' in namesWithSubunitsIndication:
                        try:
                            lPossibleSubunitsFromName.append(n.split(p)[0].split()[-1].strip(','))
                            lnameEnzyme.append(' '.join(n.split(p)[0].split()[:-1]).lower().strip())
                            break
                        except:
                            logging.info(n)
                    else:
                        try:
                            lPossibleSubunitsFromName.append(n.split(p)[1].split()[0].strip(','))
                            lnameEnzyme.append(n.split(p)[0].lower().strip())
                            break
                        except:
                            logging.info(n)
        
        sub_dict = {
            'Monomer.': '1','Homodimer.': '2','Homodimerizes.':'2','Homotrimer.':'3','Homotetramer.':'4','Homopentamer.':'5','Homohexamer.':'6','Homoheptamer.':'7',
            'Homooctamer.':'8','Octamer.':'8','Homodecamer.':'10','Homododecamer.':'12',
            'Active as a monomer':'1','Binds DNA as a monomer':'1',"Monomer \(in vitro\) \(PubMed:":'1',"Monomer \(PubMed:":'1','Monomer \(Ref.':'1',
            'Monomer in solution':'1','Monomer;':'1','Monomeric in solution.':'1','Binds DNA as monomer':'1','Forms monomers in solution':'1','Monomer \(disintegrin\)':'1',
            'Monomer \(G-actin\)':'1','Dimer of': '1', 'Heterodimer': '1','heterodimer': '1','Heterotrimer': '1','Composed of two chains': '2',
            'Forms head-to-head homodimers':'2','Forms head-to-tail homodimers':'2','Forms homodimer in solution':'2',
            'Forms homodimers':'2','Homodimer \(PubMed:':'2','Homodimer \(Ref.':'2','Homodimer \(via':'2','Homodimer in solution':'2','Homodimer of':'2','Homodimer,':'2',
            'Homodimer;':'2','Forms a homodimer':'2','Active as a homodimer.':'2','Acts as homodimer':'2','Binds DNA as a homodimer':'2','Binds DNA as homodimer':'2',
            'Can form homodimer':'2','Can homodimerize':'2','Forms a functional homodimer':'2','Forms an asymmetric homodimer':'2','Forms disulfide-linked homodimers':'2',
            'Forms homodimer':'2','Head to tail homodimer':'2','Headphone-shaped homodimer':'2','Homodimer \(in vitro\)':'2','Homodimer formed by':'2','Homodimer in':'2',
            'Homodimer that':'2','Homodimers.':'2','homodimers': '2', 'homodimer around DNA': '2','Tetramer of': '2','heterotetramer': '2',  
            'Homotrimer': '3','Homotrimer,':'3','Homotrimer;':'3','Forms homotrimers \(PubMed:':'3','Homotrimer formed of':'3','Homotrimer in solution':'3',
            'Can form homotrimer':'3','Forms homotrimers':'3','Homotrimer \(PubMed:':'3','Homotrimers of':'3','Forms homotetramers':'4','Homotetramer,':'4',
            'Homotetramer composed of':'4','Homotetramer consisting of':'4','Homotetramer formed':'4','Homotetramer in solution':'4','Homotetramer:':'4','Homotetramer;':'4',
            'A homotetramer formed':'4','Binds DNA as a homotetramer':'4','Homotetramer \(in vitro\)':'4','Homotetramer \(MAT-I\)':'4','Homotetramer \(PubMed:':'4','Heterooctamer': '4', 
            'Homotetramer \(Ref.':'4','Homotetramer in':'4','Homopentamer \(PubMed:':'5','Homopentamer with':'5','Homopentamer arranged in':'5',
            'Homopentamer;':'5','Forms a homopentamer':'5','Homopentamer \(in vitro\)':'5','Homopentamer,':'5','Homohexameric': '6', 'Homohexamer': '6', 
            'Homohexamer \(PubMed:':'6','Homohexamer composed of':'6','Homohexamer in solution':'6','Homohexamer with':'6','Homohexamer,':'6','Homohexamer;':'6',
            'Homohexameric ring arranged as':'6','A double ring-shaped homohexamer of':'6','Forms a homohexamer':'6','Forms homohexameric rings':'6','Forms homohexamers':'6',
            'Homohexamer \(dimer of homotrimers\)':'6','Homoheptamer arranged in':'7','Homoheptamer;':'7','Heterooligomer': '7', 'Forms only homooctamers':'8',
            'Homooctamer composed of':'8','Homooctamer,':'8','Homooctamer;':'8','Homooctamer \(isoform 2\)':'8','Homooctamer formed by':'8','Homooctamer of':'8',
            'Homooctomer \(PubMed:':'8', 'Homodecamer;':'10','Homodecamer composed of':'10','Forms an asymmetric tunnel-fold homodecamer':'10','Homodecamer,':'10',
            'Homodecamer consisting of':'10','Homodecamer of':'10','Homodecamer; composed of':'10', 'Homododecamer \(PubMed:':'12','Homododecamer composed of':'12',
            'Homododecamer;':'12', 'Heterotetramer': '2', 'Dodecamer': '12','Oligomer of 12 subunits': '12','24-polypeptide': '24','24 subunits': '24',
            'ClpP subunits': '7', '7 subunits': '7','Tat system': '1', 'homodecamer,': '10', 'dodecamer': '12','octamer': '8', 'cylinder': '7', 
            'cytochrome bc1': '1', 'UreD, UreF and UreG': '2','spirosomes':'40','The complex forms trimers':'3'
        }
        
        for eachkey in sub_dict.keys():
            sub_num = 'manual'
            
            if re.search(eachkey, txt_subunit):
                sub_num = sub_dict[eachkey]
                break
            
            if len(txt_subunit) == 1:
                sub_num = 'vacuum'
                break
            
            # Extract subunit number from complex txt_subunit sentences
            if sub_num == 'manual':
                possibleWords = list(set(lPossibleSubunitsFromName))
                
                if len(possibleWords) > 0:
                    for p in possibleWords:
                        namesWithSubunitsIndication = []
                        
                        try:
                            namesWithSubunitsIndication = txt_subunit[0]
                        except:
                            namesWithSubunitsIndication = txt_subunit

                        try:
                            if re.search(p, str(namesWithSubunitsIndication)):
                                if re.search('FGAM|RNAP', namesWithSubunitsIndication):
                                    try:
                                        sub_num = namesWithSubunitsIndication.split(p)[0].split()[-1].strip(',')
                                    except:
                                        sub_num == 'manual'
                                elif 'CF' in namesWithSubunitsIndication:
                                    sub_num = namesWithSubunitsIndication.split(p)[1].split()[0].strip(',')
                                    if p == 'a':
                                        sub_num = '1'
                                        break
                                else:
                                    sub_num = '1'
                                
                                if re.search('UreD', p):
                                    sub_num = '1'
                            else:
                                sub_num = '1'
                        except:
                            sub_num = '1'
                else:
                    sub_num = '1'
        
        namesWithSubunitsIndication = txt_subunit
        
        if re.search('ATP-binding proteins', namesWithSubunitsIndication):
            sub_num = '2'
        if re.search('transmembrane proteins', namesWithSubunitsIndication):
            sub_num = '2'
        if re.search('solute-binding protein', namesWithSubunitsIndication):
            sub_num = '1'
        if re.search('ribosom', namesWithSubunitsIndication):
            sub_num = '1'
        
        sub_num = sub_num.strip('.').strip('()')
        return sub_num


    def _get_subunit_data(self, batch_size: int = 10, use_progressbar = False) -> None:
        """
        Query UniProt for subunit information and annotate genes in the model
        with the number of subunits based on protein name, gene name, and comment fields.

        Args:
            batch_size (int): Number of UniProt IDs to include in each batch API call.
            use_progressbar (bool): Show the progressbar (default: False)
        
        Returns:
            None
        """
        if batch_size>500:
            raise ValueError('The batch_size cannot be larger than 500')
        uniprot_gene_id = {}
        uniprot_ids = []
        for g in self.model.genes:
            u = g.annotation.get('uniprot', None)
            if u:
                if isinstance(u, list):
                    # WARNING: we deal with multiple possible uniprot entries 
                    #by selecting the first that works.... not great
                    for y in u:
                        uniprot_gene_id[y] = g.id
                        uniprot_ids.append(y)
                elif isinstance(u, str):
                    uniprot_gene_id[u] = g.id
                    uniprot_ids.append(u)
                else:
                    logging.warning(f'Cannot recognise uniprot entry {u}')
        
        # Gather the UniProt IDs only
        logging.debug(f'uniprot_ids: {uniprot_ids}')

        # Split into batches
        batch_iterator = tqdm(range(0, len(uniprot_ids), batch_size), desc='Fetching subunit information using UniProt') if use_progressbar else range(0, len(uniprot_ids), batch_size)
        for b_it in batch_iterator:
            batch = uniprot_ids[b_it:b_it + batch_size]
            logging.debug(f'Processing batch: {batch}')
            query = " OR ".join([f'accession:{uid}' for uid in batch])
            logging.debug(f'query: {query}')
            base_url = "https://rest.uniprot.org/uniprotkb/search"
            params_reviewed = {
                "query": query,
                "format": "json",
                "fields": "protein_name,gene_names,cc_subunit",
                "size": batch_size+10,
            }

            response = requests.get(base_url, params=params_reviewed)
            data = response.json()

            uniprot_subunit_num = {}
            for uniprot_entry in data.get('results', []):
                logging.debug(f'-------------------')
                protein_names = []
                prot_desc = uniprot_entry.get('proteinDescription', {})
                for section in ['recommendedName', 'alternativeNames']:
                    entries = prot_desc.get(section, [])
                    entries = [entries] if isinstance(entries, dict) else entries
                    for entry in entries:
                        for key in ['fullName', 'shortNames']:
                            val = entry.get(key)
                            vals = val if isinstance(val, list) else [val]
                            for v in vals:
                                if v and 'value' in v:
                                    protein_names.append(v['value'])
                
                gene_names = [i.get('geneName', {}).get('value') for i in uniprot_entry.get('genes', [])]

                txt_subunit = []
                for i in uniprot_entry.get('comments', []):
                    if i.get('commentType', '').lower() == 'subunit':
                        for y in i.get('texts', []):
                            txt_subunit.append(y.get('value', ''))
                txt_subunit = ' '.join(txt_subunit)

                uniprot_id = uniprot_entry.get('primaryAccession')
                logging.debug(f'uniprot_id: {uniprot_id}')
                if uniprot_id:
                    logging.debug(f'protein_names: {protein_names}')
                    logging.debug(f'gene_names: {gene_names}')
                    logging.debug(f'txt_subunit: {txt_subunit}')
                    # Get the gene and add the annotation
                    try:
                        gene = self.model.genes.get_by_id(uniprot_gene_id[uniprot_id])
                        logging.debug(f'Found gene {gene.id}')
                        if ('number_of_subunits' in gene.annotation and gene.annotation.get('number_of_subunits')) or not 'number_of_subunits' in gene.annotation:
                            logging.debug(f'-> Adding subunit number')
                            gene.annotation['number_of_subunits'] = self._get_subunit_num(
                                protein_names, gene_names, txt_subunit
                            )
                    except KeyError:
                        logging.warning(f'Cannot find the gene for the following uniprot: {uniprot}')


    def _fetch_protein_sequence(self, uniprot_id: str) -> Optional[str]:
        """
        Fetches the amino acid sequence for a given UniProt ID.

        Args:
            uniprot_id (str): The UniProt accession ID.

        Returns:
            Optional[str]: The amino acid sequence, or None if not found.
        """
        try:
            url = (
                f"https://rest.uniprot.org/uniprotkb/search?query=accession:{uniprot_id}"
                "&format=tsv&fields=accession,sequence"
            )
            response = requests.get(url)
            response.raise_for_status()
            lines = response.text.strip().split("\n")
            if len(lines) > 1:
                return lines[1].split("\t")[1]
        except Exception:
            pass
        return None


    def _compute_mass(self, aaseq: str) -> float:
        """
        Compute the molecular mass of a protein sequence.

        Args:
            aaseq (str): Amino acid sequence.

        Returns:
            float: The molecular mass in Daltons.
        """
        cleaned_seq = aaseq.replace('X', '')
        return ProteinAnalysis(cleaned_seq).molecular_weight()


    def _get_uniprot_aaseq_mw(self, batch_size: int = 10, use_progressbar: bool = False) -> None:
        """
        Fetch amino acid sequences and molecular weights from UniProt and annotate COBRA model genes.

        Args:
            model (Any): COBRA model containing genes with UniProt annotations.
            batch_size (int): Number of UniProt entries per API batch request.
            use_progressbar (bool): Whether to display a tqdm progress bar during execution.

        Returns:
            None
        """
        if batch_size>500:
            raise ValueError('The batch_size cannot be larger than 500')
        uniprot_gene_id = {}
        uniprot_ids = []
        for g in self.model.genes:
            u = g.annotation.get('uniprot', None)
            if u:
                if isinstance(u, list):
                    # WARNING: we deal with multiple possible uniprot entries 
                    #by selecting the first that works.... not great
                    for y in u:
                        uniprot_gene_id[y] = g.id
                        uniprot_ids.append(y)
                elif isinstance(u, str):
                    uniprot_gene_id[u] = g.id
                    uniprot_ids.append(u)
                else:
                    logging.warning(f'Cannot recognise uniprot entry {u}')
        
        # Gather the UniProt IDs only
        logging.debug(f'uniprot_ids: {uniprot_ids}')

        base_url = "https://rest.uniprot.org/uniprotkb/search"
        batch_iterator = tqdm(
            range(0, len(uniprot_ids), batch_size),
            desc='Fetching protein information using UniProt'
        ) if use_progressbar else range(0, len(uniprot_ids), batch_size)

        # Split into batches
        for b_it in batch_iterator:
            batch = uniprot_ids[b_it:b_it + batch_size]
            query = " OR ".join([f'accession:{uid}' for uid in batch]) 
            # Send UniProt API request for sequence and molecular weight
            params_reviewed = {
                "query": query,
                "format": "json",
                "fields": "accession,sequence,mass",
                "size": batch_size + 10
            }

            response = requests.get(base_url, params=params_reviewed)
            data = response.json()

            for res_entry in data.get("results", []):
                uniprot_id = res_entry.get('primaryAccession')
                sequence = res_entry.get('sequence', {}).get('value')
                mw = res_entry.get('sequence', {}).get('molWeight')

                # If molWeight is not given, compute it
                if not isinstance(mw, (int, float)) and sequence:
                    mw = self._compute_mass(sequence)

                if not isinstance(mw, (int, float)):
                    logging.warning(f'The unitprot {uniprot_id} does not have valid mw {mw}')
                    continue

                if uniprot_id:
                    try:
                        gene = self.model.genes.get_by_id(uniprot_gene_id[uniprot_id])
                        if ('aaseq' in gene.annotation and gene.annotation.get('number_of_subunits')) or not 'aaseq' in gene.annotation:
                            gene.annotation['aaseq'] = sequence
                        if ('mw' in gene.annotation and gene.annotation.get('mw')) or not 'mw' in gene.annotation:
                            gene.annotation['mw'] = mw
                        # Optionally compute total molecular weight based on subunit number
                        subunit_number = gene.annotation.get('number_of_subunits')
                        if subunit_number:
                            if ('total_subunits_mw' in gene.annotation and gene.annotation.get('total_subunits_mw')) or not 'total_subunits_mw' in gene.annotation:
                                try:
                                    gene.annotation["total_subunits_mw"] = float(mw) * float(subunit_number)
                                except ValueError:
                                    logging.warning(f'Could not convert {mw} or {subunit_number}')
                    except KeyError:
                        logging.warning(f'UniProt ID {uniprot_id} not found in gene map')


    def _convert_to_irreversible(
            self, 
            input_model,
            substrate_exchange_reaction: Reaction,
            use_progressbar: bool = False) -> None:
        """
        Convert all reversible reactions in a COBRA model to irreversible format.

        This function modifies the model in place by:
        - Splitting fully reversible reactions into forward and reverse irreversible reactions.
        - Adjusting reaction bounds to ensure only positive flux.
        - Tagging reverse reactions using the `reflection` note.
        - Copying annotations, gene-reaction rules, and subsystems.

        Args:
            model (Model): The COBRA model to convert. This object is modified in place.

        Returns:
            None
        """
        reactions_to_add = []
        coefficients = {}
        for reaction in input_model.reactions:
            if reaction.lower_bound < 0 and reaction.upper_bound == 0:
                for metabolite in reaction.metabolites:
                    original_coefficient = reaction.get_coefficient(metabolite)
                    #QUESTION: original has -2 not sure why we would double that?
                    reaction.add_metabolites({metabolite: -2 * original_coefficient})
                try:
                    reaction.id += "_reverse"
                    reaction.upper_bound = -reaction.lower_bound
                    reaction.lower_bound = 0
                except ValueError:
                    pass

            if reaction.lower_bound < 0 and reaction.upper_bound > 0:
                reverse_reaction = Reaction(reaction.id + "_reverse")
                if substrate_exchange_reaction==reaction:
                    reverse_reaction.lower_bound = max(0, 0)
                else:
                    reverse_reaction.lower_bound = max(0, -reaction.upper_bound)
                reverse_reaction.upper_bound = -reaction.lower_bound
                coefficients[reverse_reaction] = reaction.objective_coefficient * -1
                reaction.lower_bound = max(0, reaction.lower_bound)
                reaction.upper_bound = max(0, reaction.upper_bound)
                reaction.notes["reflection"] = reverse_reaction.id
                reverse_reaction.notes["reflection"] = reaction.id
                reaction_dict = {k: v * -1 for k, v in reaction._metabolites.items()}
                reverse_reaction.add_metabolites(reaction_dict)
                reverse_reaction._model = reaction._model
                reverse_reaction._genes = reaction._genes
                for gene in reaction._genes:
                    gene._reaction.add(reverse_reaction)
                reverse_reaction.subsystem = reaction.subsystem
                reverse_reaction.gene_reaction_rule = reaction.gene_reaction_rule
                try:
                    reaction.annotation
                except:
                    pass
                else:
                    reverse_reaction.annotation = reaction.annotation
                reactions_to_add.append(reverse_reaction)

        input_model.add_reactions(reactions_to_add)
        #QUESTION: why set the objective to add the coefficients?
        #logging.debug(coefficients)
        set_objective(input_model, coefficients, additive=True)


    def _isoenzyme_split(self, input_model):
        """Split isoenzyme reaction to mutiple reaction

        Arguments
        ----------
        * model: cobra.Model.
        
        :return: new cobra.Model.
        """  
        for r in input_model.reactions:
            if re.search(" or ", r.gene_reaction_rule):
                rea = r.copy()
                gene = r.gene_reaction_rule.split(" or ")
                subunits_ids = []
                for index, value in enumerate(gene):
                    if index == 0:
                        r.id = r.id + "_num1"
                        subunits_ids.append(r.id)
                        r.gene_reaction_rule = value
                    else:
                        r_add = rea.copy()
                        r_add.id = rea.id + "_num" + str(index+1)
                        r_add.gene_reaction_rule = value
                        subunits_ids.append(r_add.id)
                        #model.add_reaction(r_add)#3.7
                        input_model.add_reactions([r_add])#3.8
                for sub_id in subunits_ids:
                   tmp_r = input_model.reactions.get_by_id(sub_id)
                   tmp_r.notes['isoenzymes'] = [i for i in subunits_ids if i!=sub_id]
        for r in input_model.reactions:
            r.gene_reaction_rule = r.gene_reaction_rule.strip("( )")
            #TODO try this
            #r.gene_reaction_rule = r.gene_reaction_rule.removeprefix('(').removesuffix(')')


    ####### Database kcat search ####o

    def _search_kcat(
        self,
        species_name: str, 
        inchikey_levels: int = 2,
        join_mode='mean',
        closest_taxonomy=True,
        use_progressbar: bool = False
    ) -> Dict[str, Any]:
        """Search for kcat values for each reaction in a model using BRENDA and SABIORK databases.

        This function attempts to retrieve catalytic constants (kcat) by matching the reaction's
        EC number, substrates (via InChIKey), and associated genes (via UniProt) to database entries.
        Matching is done with progressive relaxation of input constraints.

        Args:
            model: COBRApy model object containing reactions and annotated genes.
            species_name (str): Organism name used to prioritize matches in BRENDA and SABIORK.
            inchikey_levels (int): Number of InChIKey layers to match (default is 2).
            join_kcat_mode (str): How to join the kcat values (default is mean). Can be max, mean, median

        Returns:
            Dict[str, Any]: Mapping of reaction IDs to matched kcat values and metadata.
        """
        reac_iterator = tqdm(self.model.reactions, desc='Database lookup of reaction kcat') if use_progressbar else self.model.reactions
        for r in reac_iterator:
            logging.debug(f'------ {r.id} ------')
            ### Gather the information
            uniprot = []
            for g in r.genes:
                u = g.annotation.get('uniprot', [])
                if isinstance(u, list):
                    uniprot += u
                elif isinstance(u, str):
                    uniprot.append(u)
                else:
                    logging.warning(f'Cannot find out the type of {u}')
            uniprot = [x for x in uniprot if x is not pd.isna(x)]
            uniprot = list(set(uniprot))
            brenda_ec_code = None
            model_ec_code = r.annotation.get('ec-code')
            if model_ec_code:
                if isinstance(model_ec_code, list):
                    #WARNING: only using one that works if there are multiple
                    for e in model_ec_code:
                        try:
                            brenda_ec_code = self._find_source_node(e)
                            if brenda_ec_code:
                                break
                        except KeyError:
                            pass
                elif isinstance(model_ec_code, str):
                    try:
                        brenda_ec_code = self._find_source_node(model_ec_code)
                    except KeyError:
                        pass
                else:
                    logging.warning(f'Cannot detect the type of {model_ec_code}: {type(model_ec_code)}')
            ### get the substrates
            substrate_inchikey = []
            inchikey = None
            for m in r.reactants:
                inchikey = m.annotation.get('inchikey')
                if not inchikey:
                    mnxm = m.annotation.get('metanetx.chemical')
                    if isinstance(mnxm, list):
                        for i in mnxm:
                            try:
                                inchikey = self.mnxm_inchikey[i]
                            except KeyError:
                                inchikey = None
                    elif isinstance(mnxm, str): 
                        try:
                            inchikey = self.mnxm_inchikey[mnxm]
                        except KeyError:
                            inchikey = None
                if inchikey:
                    if not self._is_inchikey_cofactor(inchikey):
                        substrate_inchikey.append(inchikey)
            substrate_inchikey = tuple(sorted(substrate_inchikey))
            logging.debug(f'uniprot: {uniprot}')
            logging.debug(f'brenda_ec_code: {brenda_ec_code}')
            logging.debug(f'substrate_inchikey: {substrate_inchikey}')
            logging.debug(f'species_name: {self.species_name}')

            #these are all the combinations that we will use to try to fetch the
            #closest kcat entry in the brenda database
            param_combinations = [
                # Full set
                {"species_name": self.species_name, "uniprot": uniprot, "substrate_inchikey": substrate_inchikey, "closest_taxonomy": False},
                # Drop substrate_inchikey
                {"species_name": self.species_name, "uniprot": uniprot, "substrate_inchikey": None, "closest_taxonomy": False},
                # Drop uniprot
                {"species_name": self.species_name, "uniprot": None, "substrate_inchikey": substrate_inchikey, "closest_taxonomy": False},
                # Drop both uniprot and substrate_inchikey
                {"species_name": self.species_name, "uniprot": None, "substrate_inchikey": None, "closest_taxonomy": False},
                #### cloesest_taxonmy
                {"species_name": None, "uniprot": uniprot, "substrate_inchikey": substrate_inchikey, "closest_taxonomy": True},
                # Drop substrate_inchikey
                {"species_name": None, "uniprot": uniprot, "substrate_inchikey": None, "closest_taxonomy": True},
                # Drop uniprot
                {"species_name": None, "uniprot": None, "substrate_inchikey": substrate_inchikey, "closest_taxonomy": True},
                # Drop both uniprot and substrate_inchikey
                {"species_name": None, "uniprot": None, "substrate_inchikey": None, "closest_taxonomy": True},
                # Drop species_name, keep others
                {"species_name": None, "uniprot": uniprot, "substrate_inchikey": substrate_inchikey, "closest_taxonomy":  False},
                # Drop species_name and substrate
                {"species_name": None, "uniprot": uniprot, "substrate_inchikey": None, "closest_taxonomy": False},
                # Drop species_name and uniprot
                {"species_name": None, "uniprot": None, "substrate_inchikey": substrate_inchikey, "closest_taxonomy": False},
                # All None
                {"species_name": None, "uniprot": None, "substrate_inchikey": None, "closest_taxonomy": False},
            ]

            if brenda_ec_code:
                # KCAT
                #if you already have some
                kcat_res = {}
                search_for_kcat = True
                if search_for_kcat:
                    ## BRENDA
                    for params in param_combinations:
                        if search_for_kcat:
                            try:
                                logging.debug(params)
                                logging.debug(brenda_ec_code)
                                kcat_res = self._filter_kinetic_entry(
                                    self.brenda_ec_inchikey_kcat[brenda_ec_code], 
                                    species_name=params['species_name'],
                                    substrate_inchikey=params['substrate_inchikey'],
                                    uniprot=params['uniprot'],
                                    closest_taxonomy=params['closest_taxonomy'],
                                    inchikey_levels=inchikey_levels,
                                )
                                if kcat_res:
                                    logging.debug(f'BRENDA kcat_res: {kcat_res}')
                                    search_for_kcat = False
                            except KeyError:
                                pass
                    ## SabioRK
                    if not kcat_res:
                        #WARNING: only returning one 
                        #TODO: find the best EC number based on the list
                        sabiork_entry = self._search_kcat_sabiork(ec_number=brenda_ec_code)
                        if sabiork_entry:
                            for params in param_combinations:
                                try:
                                    logging.debug(params)
                                    kcat_res = self._filter_kinetic_entry(
                                        sabiork_entry,
                                        species_name=params['species_name'],
                                        substrate_inchikey=params['substrate_inchikey'],
                                        uniprot=params['uniprot'],
                                        inchikey_levels=inchikey_levels,
                                        closest_taxonomy=params['closest_taxonomy'],
                                    )
                                    if kcat_res:
                                        logging.debug(f'SABIORK kcat_res: {kcat_res}')
                                        search_for_kcat = False
                                        break
                                except (KeyError, TypeError) as e:
                                    pass
                if kcat_res:
                    #find the closest taxonomic member
                    kcat, kcat_std = self._summarize_values(kcat_res, mode=join_mode)
                    if kcat:
                        r.annotation['DBkcat'] = {'kcat': kcat, 'kcat_std': kcat_std}
                # SA
                search_for_sa = True
                sa_res = {}
                if search_for_sa:
                    for params in param_combinations:
                        if search_for_sa:
                            try:
                                sa_res = self._filter_kinetic_entry(
                                    self.brenda_ec_inchikey_sa[brenda_ec_code], 
                                    species_name=params['species_name'],
                                    substrate_inchikey=params['substrate_inchikey'],
                                    uniprot=params['uniprot'],
                                    inchikey_levels=inchikey_levels,
                                    closest_taxonomy=params['closest_taxonomy'],
                                )
                                logging.debug(f'BRENDA sa_res: {sa_res}')
                                if sa_res:
                                    search_for_sa = False
                            except KeyError:
                                pass
                if sa_res:
                    sa, sa_std = self._summarize_values(sa_res, mode=join_mode)
                    if sa:
                        r.annotation['DBsa'] = {'sa': sa, 'sa_std': sa_std}

    ## parse brenda db

    def _json_brenda_to_G(self):
        def extract_history_ec(text):
            match = re.search(r'EC\s+(\d+\.\d+\.\d+\.\d+)', text)
            if match:
                ec_number = match.group(1)
                return ec_number
            return None
        
        G = nx.DiGraph()
        
        for i in self.json_brenda['data'].keys():
            G.add_node(i)
            
        for i in self.json_brenda['data'].keys():
            if 'history' in self.json_brenda['data'][i]:
                to_ec = extract_history_ec(
                    self.json_brenda['data'][i]['history']
                )
                if to_ec:
                    G.add_edge(to_ec, i)
        return G


    def _find_source_node(self, node) -> str:
        """
        Find the furthest upstream source (node with in-degree 0) from a given node.

        Args:

            node: The starting node.

        Returns:
            The source node if reachable, otherwise None.
        """
        current = node
        visited = set()
        if node not in self.brenda_ec_G:
            raise KeyError(f'{node} is not in the graph')
        while True:
            preds = list(self.brenda_ec_G.predecessors(current))
            if not preds:
                return current
            if current in visited:
                return None  # Avoid infinite loops (e.g., in cycles)
            visited.add(current)
            current = preds[0]  # Assumes a tree or DAG with one upstream path


    def _is_inchikey_cofactor(self, inchikey):
        if inchikey in self.cofactors_inchikey_layers3:
            return True
        if inchikey_layer_extract(inchikey, 2) in self.cofactors_inchikey_layers2:
            return True
        return False


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
                        if not self._is_inchikey_cofactor(inchikey):
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
                    reactions_list = [y for y in self.json_brenda['data'][ec_number]['substrates_products'] if set(protein_ids) & set(y['proteins'])]
                else:
                    reactions_list = self.json_brenda['data'][ec_number]['substrates_products']
        if not reactions_list:
            if 'reaction' in self.json_brenda['data'][ec_number]:
                if protein_ids:
                    reactions_list = [y for y in self.json_brenda['data'][ec_number]['reaction'] if set(protein_ids) & set(y['proteins'])]
                else:
                    reactions_list = self.json_brenda['data'][ec_number]['reaction']
        if 'cofactor' in self.json_brenda['data'][ec_number]:
            if protein_ids:
                cofactor_list = [y for y in self.json_brenda['data'][ec_number]['cofactor'] if set(protein_ids) & set(y['proteins'])]
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


    def _find_matching_keys(
        self,
        input_keys: List[Tuple], 
        input_query: Any, 
        inchikey_layers=3) -> List[Tuple]:
        """
        Find all keys that contain or are equal to the query.

        Args:
            keys: List of tuple or nested tuple keys.
            query: A string, tuple, or tuple of tuples to search for.

        Returns:
            A list of matching keys.
        """
        matches = []
        query = None
        if isinstance(input_query, str):
            query = inchikey_layer_extract(input_query, inchikey_layers)
        elif isinstance(input_query, tuple):
            query = tuple(sorted([inchikey_layer_extract(i, inchikey_layers) for i in input_query]))
        elif isinstance(query, (list, tuple)) and all(isinstance(t, tuple) for t in query):
            query = tuple([tuple(sorted([inchikey_layer_extract(y, inchikey_layers) for y in i if y])) for i in input_query if i])
        keys = []
        for k in input_keys:
            keys.append(
                tuple([tuple(sorted([inchikey_layer_extract(y, inchikey_layers) for y in i if y])) for i in k if i]))
        
        for key, ori_key in zip(keys, input_keys):
            # Exact match
            if key == query:
                matches.append(ori_key)
                continue
            
            # Flatten the nested tuple structure for search
            flattened = [item for sub in key for item in (sub if isinstance(sub, tuple) else (sub,))]
            
            if isinstance(query, str):
                if query in flattened:
                    matches.append(ori_key)
            elif isinstance(query, tuple):
                if all(elem in flattened for elem in query):
                    matches.append(ori_key)
            elif isinstance(query, (list, tuple)) and all(isinstance(t, tuple) for t in query):
                # query is a tuple of tuples (nested), compare by flattening
                query_flat = [item for sub in query for item in sub]
                if all(elem in flattened for elem in query_flat):
                    matches.append(ori_key)
        return matches


    def _filter_kinetic_entry(
        self,
        input_entry_ec_number, 
        species_name=None, 
        closest_taxonomy=False,
        substrate_inchikey=None,
        uniprot=None,
        inchikey_levels=3,
    ):
        entry_ec_number = copy.deepcopy(input_entry_ec_number)
        #### check the species name 
        if species_name:
            new_entry_ec_number = {}
            for substrates in entry_ec_number:
                if species_name in entry_ec_number[substrates]:
                    if substrates not in new_entry_ec_number:
                        new_entry_ec_number[substrates] = {}
                    new_entry_ec_number[substrates][species_name] = entry_ec_number[substrates][species_name]
            entry_ec_number = new_entry_ec_number
        #find the closest taxonomic member
        if closest_taxonomy:
            other_orgs = list(set([y for i in entry_ec_number for y in entry_ec_number[i].keys()]))
            if other_orgs:
                try:
                    to_keep_org = self._rank_taxonomic_proximity(
                        target=self.species_name,
                        others=other_orgs,
                    )[0][0]
                    if to_keep_org:
                        entry_ec_number = {i: {y: entry_ec_number[i][y]} for i in entry_ec_number for y in entry_ec_number[i] if y==to_keep_org}
                except IndexError:
                    pass
        #### if substrate is input check
        if substrate_inchikey:
            closest_substrates = self._find_matching_keys(
                input_keys=entry_ec_number.keys(), 
                input_query=substrate_inchikey, 
                inchikey_layers=inchikey_levels)
            if closest_substrates:
                logging.debug(f'substrate_inchikey: {substrate_inchikey}')
                logging.debug(f'closest_substrates: {closest_substrates}')
                entry_ec_number = {key: entry_ec_number[key] for key in closest_substrates if key in entry_ec_number}
        ### check the uniprot ####
        if uniprot:
            new_entry_ec_number = copy.deepcopy(entry_ec_number)
            for substrates in entry_ec_number:
                for species_n in entry_ec_number[substrates]: 
                    if isinstance(uniprot, str):
                        if not uniprot in entry_ec_number[substrates][species_n]['uniprot']:
                            del new_entry_ec_number[substrates][species_n]
                    elif isinstance(uniprot, list):
                        if not set(uniprot) & set(entry_ec_number[substrates][species_n]['uniprot']):
                            del new_entry_ec_number[substrates][species_n]
            entry_ec_number = new_entry_ec_number
        #### remove empty values
        res = {}
        for i in entry_ec_number:
            if entry_ec_number[i]:
                res[i] = entry_ec_number[i]
        return res


    def _search_kcat_sabiork(
        self,
        ec_number,
        kegg_id=None
    ) -> dict:
        sabiork_wait_time = 1.5
        sabiork_query_url = "http://sabiork.h-its.org/sabioRestWebServices/kineticlawsExportTsv"
        sabiork_unit_multiplier: Dict[str, float] = {
            "s^(-1)": 1.0,
            "min^(-1)": 1/60,
            "h^(-1)": 1/(60*60),
        }
 
        if ec_number:
            q = {'fields[]': ['ECNumber', 'KeggReactionID', 'Organism', 'Parameter', 'Substrate', 'uniprotID'], 'q': f'(ECNumber:{ec_number} AND Parametertype:kcat AND EnzymeType:wildtype)'}
        elif kegg_id:
            q = {'fields[]': ['ECNumber', 'KeggReactionID', 'Organism', 'Parameter', 'Substrate', 'uniprotID'], 'q': f'(KeggReactionID:{kegg_id} AND Parametertype:kcat AND EnzymeType:wildtype)'}
        else:
            raise KeyError('Must input ec_number {ec_number} or kegg_id {kegg_id}')

        try:
            #request = requests.post(sabiork_query_url, params=q, verify=False)
            request = requests.post(sabiork_query_url, params=q)
        except Exception:
            time.sleep(sabiork_wait_time)
            return {}
        try:
            request.raise_for_status()
        except Exception:
            #logging.info("SABIO-RK API error with query:")
            #logging.info(query_string)
            time.sleep(sabiork_wait_time)
            return {}
        time.sleep(sabiork_wait_time)
        ### Parse the results
        df = pd.read_csv(StringIO(request.text), sep="\t")
        df_kcat = df[df['parameter.type']=='kcat'].copy()
        df_kcat.loc[:, "adjusted_kcat"] = df_kcat["parameter.startValue"] * df_kcat["parameter.unit"].map(sabiork_unit_multiplier)
        res = []
        for i in df_kcat['Substrate']:
            r = []
            for y in i.split(';'):
                inchikey = self.molecule_name_search_inchikey(y.lower())
                if inchikey:
                    r.append(inchikey)
            if r:
                res.append(tuple(sorted(r)))
            else:
                res.append('no_identifiable_substrate')
        df_kcat.loc[:, 'inchikey_substrate'] = res

        #return it as a dict
        result = {}
        for _, row in df_kcat.iterrows():
            key = (row['inchikey_substrate'],)  # Make it a tuple of tuples
            if not key in result:
                result[key] = {}
            organism = row['Organism']
            if not organism in result[key]:
                result[key][organism] = {'uniprot': [], 'values': [], 'comments': 'sabiork'}
            result[key][organism]['uniprot'].append(row['uniprotID'])
            result[key][organism]['values'].append(row['adjusted_kcat'])
            value = {
                organism: {
                    'uniprot': [row['uniprotID']],
                    'values': [row['adjusted_kcat']],
                    'comments': []  # You can populate this if data available
                }
            }
            result[key] = value

        return result


    def _get_species_name(
        self,
        taxid: int,
    ) -> str:
        """Return the species-level name from a taxonomy ID."""
        ncbi = NCBITaxa()
        lineage = ncbi.get_lineage(taxid)
        names = ncbi.get_taxid_translator(lineage)
        ranks = ncbi.get_rank(lineage)

        for tid in lineage:
            if ranks[tid] == "species":
                return names[tid]
        return names.get(taxid, "Unknown")


    def _rank_taxonomic_proximity(
        self,
        target: str, 
        others: List[str]
    ) -> List[Tuple[str, int]]:
        """
        Rank species_names by taxonomic proximity to a target species_name.

        Args:
            target: Scientific name of the target species_name (e.g., "Homo sapiens").
            others: List of other species_name names to compare.

        Returns:
            A list of tuples (species_name_name, distance), sorted by ascending distance.
        """
        ncbi = NCBITaxa()
        # Get taxids
        try:
            target_taxid = ncbi.get_name_translator([target])[target][0]
        except KeyError:
            raise ValueError(f"Target species_name '{target}' not found.")

        other_taxids = {}
        for name in others:
            try:
                if name not in other_taxids:
                    other_taxids[name] = ncbi.get_name_translator([name])[name][0]
            except KeyError as e:
                logging.warning(f"Organism not found: {e}")

        # Get full lineages
        target_lineage = ncbi.get_lineage(target_taxid)
        lineages = {
            name: ncbi.get_lineage(taxid)
            for name, taxid in other_taxids.items()
        }

        # Compute distance as inverse depth of common ancestor
        ranked = []
        for name, lineage in lineages.items():
            lca = set(target_lineage).intersection(lineage)
            distance = len(set(target_lineage + lineage)) - 2 * len(lca)
            ranked.append((name, distance))

        # Sort by ascending distance (closer means smaller distance)
        ranked.sort(key=lambda x: x[1])
        return ranked


    def _summarize_values(
        self,
        substrate_data: Dict,
        mode: str = 'mean',
    ) -> Dict:
        """Summarize values values for each reaction entry in a nested dictionary.

        Args:
            data: Nested dictionary of the format {reaction: {inchikeys: {species_name: {'values': [...]}}}}
            mode: Statistic to compute from ['mean', 'median', 'max']

        Returns:
            A dictionary of the form:
            {
                'reaction_id': {
                    'kcat_value': float,
                    'std_dev': float
                },
                ...
            }
        """
        assert mode in {'mean', 'median', 'max'}, "Mode must be 'mean', 'median', or 'max'"

        all_values = []

        for substrates, org_dict in substrate_data.items():
            for org, val_dict in org_dict.items():
                values = val_dict.get('values', [])
                all_values.extend(values)

        if not all_values:
            return None, None

        all_values = [i for i in all_values if i]
        all_values = [i for i in all_values if i>=0.0]
        if not all_values:
            return None, None
        all_values_array = np.array(all_values, dtype=float)
        if mode == 'mean':
            value = np.mean(all_values_array)
        elif mode == 'median':
            value = np.median(all_values_array)
        elif mode == 'max':
            value = np.max(all_values_array)

        if pd.isna(value):
            return None, None

        return float(value), float(np.std(all_values_array))


