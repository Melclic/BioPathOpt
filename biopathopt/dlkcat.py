import math
import torch
import pickle
import logging
from rdkit import Chem
from typing import Any
import os
import sys
import numpy as np
from collections import defaultdict

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from sklearn.metrics import mean_squared_error,r2_score

from biopathopt.utils import load_pickle
from biopathopt import Data

class KcatPredictor(Data):
    def __init__(
        self,
        radius: int = 2,
        ngram: int = 3,
        dim: int = 20,
        layer_gnn: int = 3,
        window: int = 11,
        layer_cnn: int = 3,
        layer_output: int = 3,
        lr: float = 1e-3,
        lr_decay: float = 0.5,
        decay_interval: int = 10,
        weight_decay: float = 1e-6,
        iteration: int = 50,
    ):
        """
        Initializes the KcatPredictor with necessary models, dictionaries, and utility functions.

        Args:
            model_module (Any): Module containing KcatPrediction class and load_pickle function.
            Predictor (Any): Predictor wrapper for model inference.
            get_smiles (Callable): Function to get SMILES from metabolite name.
            create_atoms (Callable): Function to extract atom features.
            create_ijbonddict (Callable): Function to build bond dictionary.
            extract_fingerprints (Callable): Function to extract molecular fingerprints.
            create_adjacency (Callable): Function to create adjacency matrix.
            split_sequence (Callable): Function to tokenize sequence into n-grams.
            dim (int): Feature dimensionality.
            layer_gnn (int): Number of GNN layers.
            window (int): Window size for CNN.
            layer_cnn (int): Number of CNN layers.
            layer_output (int): Number of output layers.
            radius (int): Radius for fingerprint extraction.
            ngram (int): Size of n-grams for sequence splitting.
            weights_path (str): Path to the pretrained weights.
        """

        super().__init__()
        self.edge_dict = load_pickle(f'./data/edge_dict.pickle')
        self.fingerprint_dict = load_pickle(f'./data/fingerprint_dict.pickle')
        self.atom_dict = load_pickle(f'./data/atom_dict.pickle')
        self.bond_dict = load_pickle(f'./data/bond_dict.pickle')
        self.word_dict = load_pickle(f'./data/sequence_dict.pickle')

        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        torch.manual_seed(1234)
        torch_model_path = f'./data/all--radius{radius}--ngram{ngram}--dim{dim}--layer_gnn{layer_gnn}--window{window}--layer_cnn{layer_cnn}--layer_output{layer_output}--lr{"{:.0e}".format(lr).replace("e-0", "e-").replace("e+0", "e+") }--lr_decay{lr_decay}--decay_interval{decay_interval}--weight_decay{"{:.0e}".format(weight_decay).replace("e-0", "e-").replace("e+0", "e+")}--iteration{iteration}'
        #kcat_model.load_state_dict(torch.load(f'{run_file}data/all--radius2--ngram3--dim20--layer_gnn3--window11--layer_cnn3--layer_output3--lr1e-3--lr_decay0.5--decay_interval10--weight_decay1e-6--iteration50', map_location=device))
        if not os.path.isfile(torch_model_path):
            raise FileNotFoundError(f'Could not find {torch_model_path}')

        self.kcat_gnn_model = KcatPrediction(
            self.device,
            len(self.fingerprint_dict),
            len(self.word_dict),
            dim,
            layer_gnn,
            window,
            layer_cnn,
            layer_output
        ).to(self.device)
        self.kcat_gnn_model.load_state_dict(
            torch.load(
                torch_model_path,
                map_location=self.device
            )
        )

        self.radius = radius
        self.ngram = ngram


    def split_sequence(self, sequence, ngram,word_dict):
        sequence = '-' + sequence + '='
        # logging.debug(sequence)
        # words = [word_dict[sequence[i:i+ngram]] for i in range(len(sequence)-ngram+1)]
    
        words = list()
        for i in range(len(sequence)-ngram+1) :
            try :
                words.append(word_dict[sequence[i:i+ngram]])
            except :
                word_dict[sequence[i:i+ngram]] = 0
                words.append(word_dict[sequence[i:i+ngram]])
    
        return np.array(words)
        # return word_dict
    
    def create_atoms(self, mol,atom_dict):
        """Create a list of atom (e.g., hydrogen and oxygen) IDs
        considering the aromaticity."""
        # atom_dict = defaultdict(lambda: len(atom_dict))
        atoms = [a.GetSymbol() for a in mol.GetAtoms()]
        # logging.debug(atoms)
        for a in mol.GetAromaticAtoms():
            i = a.GetIdx()
            atoms[i] = (atoms[i], 'aromatic')
        atoms = [atom_dict[a] for a in atoms]
        # atoms = list()
        # for a in atoms :
        #     try: 
        #         atoms.append(atom_dict[a])
        #     except :
        #         atom_dict[a] = 0
        #         atoms.append(atom_dict[a])
    
        return np.array(atoms)
    
    def create_ijbonddict(self, mol,bond_dict):
        """Create a dictionary, which each key is a node ID
        and each value is the tuples of its neighboring node
        and bond (e.g., single and double) IDs."""
        # bond_dict = defaultdict(lambda: len(bond_dict))
        i_jbond_dict = defaultdict(lambda: [])
        for b in mol.GetBonds():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            bond = bond_dict[str(b.GetBondType())]
            i_jbond_dict[i].append((j, bond))
            i_jbond_dict[j].append((i, bond))
        return i_jbond_dict
    
    def extract_fingerprints(self, atoms, i_jbond_dict, radius,fingerprint_dict,edge_dict):
        """Extract the r-radius subgraphs (i.e., fingerprints)
        from a molecular graph using Weisfeiler-Lehman algorithm."""
    
        # fingerprint_dict = defaultdict(lambda: len(fingerprint_dict))
        # edge_dict = defaultdict(lambda: len(edge_dict))
    
        if (len(atoms) == 1) or (radius == 0):
            logging.debug(f'atoms: {atoms}')
            fingerprints = [fingerprint_dict[a] for a in atoms]
    
        else:
            nodes = atoms
            i_jedge_dict = i_jbond_dict
    
            for _ in range(radius):
    
                """Update each node ID considering its neighboring nodes and edges
                (i.e., r-radius subgraphs or fingerprints)."""
                fingerprints = []
                for i, j_edge in i_jedge_dict.items():
                    neighbors = [(nodes[j], edge) for j, edge in j_edge]
                    fingerprint = (nodes[i], tuple(sorted(neighbors)))
                    # fingerprints.append(fingerprint_dict[fingerprint])
                    # fingerprints.append(fingerprint_dict.get(fingerprint))
                    try :
                        fingerprints.append(fingerprint_dict[fingerprint])
                    except :
                        fingerprint_dict[fingerprint] = 0
                        fingerprints.append(fingerprint_dict[fingerprint])
    
                nodes = fingerprints
    
                """Also update each edge ID considering two nodes
                on its both sides."""
                _i_jedge_dict = defaultdict(lambda: [])
                for i, j_edge in i_jedge_dict.items():
                    for j, edge in j_edge:
                        both_side = tuple(sorted((nodes[i], nodes[j])))
                        # edge = edge_dict[(both_side, edge)]
                        # edge = edge_dict.get((both_side, edge))
                        try :
                            edge = edge_dict[(both_side, edge)]
                        except :
                            edge_dict[(both_side, edge)] = 0
                            edge = edge_dict[(both_side, edge)]
    
                        _i_jedge_dict[i].append((j, edge))
                i_jedge_dict = _i_jedge_dict
    
        return np.array(fingerprints)
    
    def create_adjacency(self, mol):
        adjacency = Chem.GetAdjacencyMatrix(mol)
        return np.array(adjacency)
    
    def predict_kcat(self, smiles, sequence, inchi=None):
        if not sequence:
            logging.warning(f'Sequence is empty: {sequence}, cannot predict kcat')
            return None
        if not smiles or not inchi:
            logging.warning(f'Need at least one smiles or inchi, cannot predict kcat')
            return None
        try:
            mol = Chem.AddHs(Chem.inchi.MolFromInchi(inchi))
        except TypeError:
            logging.warning('Cannot use inchi, switch to smiles')
            try:
                mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
            except TypeError:
                logging.error('Cannot use the passed structures')
                return None

        try:
            atoms = self.create_atoms(mol, self.atom_dict)
        except KeyError:
            logging.warning(f'Cannot create atoms: {self.atom_dict}')
            return None
        i_jbond_dict = self.create_ijbonddict(mol, self.bond_dict)
        try:
            fingerprints = self.extract_fingerprints(
                atoms, i_jbond_dict, self.radius, self.fingerprint_dict, self.edge_dict
            )
        except IndexError as e:
            logging.warning('Cannot make fingerprint for {smiles}, {inchi}, {sequence}')
            return None
        adjacency = self.create_adjacency(mol)
        words = self.split_sequence(sequence, self.ngram, self.word_dict)
        inputs = [
            torch.LongTensor(fingerprints),
            torch.FloatTensor(adjacency),
            torch.LongTensor(words)
        ]
        try:
            prediction = self.kcat_gnn_model.forward(inputs)
        except RuntimeError:
            logging.warning('Cannot get forward prediction for {inputs}')
            return None
        kcat_log_value = prediction.item()
        kcat_sec_value = float('%.4f' % math.pow(2, kcat_log_value))
        #keep only the largest (i.e. slowest) reaction subunit
        #if kcat_sec_value>kcat:
        #    kcat = kcat_sec_value
        return kcat_sec_value

    '''
    def model_predict_kcat(self, model: Any) -> None:
        """
        Predict Kcat values and assign them to each reaction's annotation in the model.

        Args:
            model (Any): COBRA model to process.

        Returns:
            None
        """
        for r in model.reactions:
            if r.gene_reaction_rule and len(r.reactants) > 1:
                kcat = 0
                totalmass = 0
                for g in r.genes:
                    totalmass += g.annotation.get('subunitmass', 0)
                    for m in r.reactants:
                        if 'metanetx.chemical' in m.annotation:
                            if m.annotation['metanetx.chemical'] not in self.mnxm_cofactors:
                                continue
                        smiles = m.annotation.get('smiles')
                        sequence = g.annotation.get('aaseq')
                        if not smiles or not sequence:
                            continue
                        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
                        atoms = self.create_atoms(mol, self.atom_dict)
                        i_jbond_dict = self.create_ijbonddict(mol, self.bond_dict)
                        fingerprints = self.extract_fingerprints(
                            atoms, i_jbond_dict, self.radius, self.fingerprint_dict, self.edge_dict
                        )
                        adjacency = self.create_adjacency(mol)
                        words = self.split_sequence(sequence, self.ngram, self.word_dict)
                        inputs = [
                            torch.LongTensor(fingerprints),
                            torch.FloatTensor(adjacency),
                            torch.LongTensor(words)
                        ]
                        prediction = self.kcat_gnn_model.forward(inputs)
                        kcat_log_value = prediction.item()
                        kcat_sec_value = float('%.4f' % math.pow(2, kcat_log_value))
                        #keep only the largest (i.e. slowest) reaction subunit
                        if kcat_sec_value>kcat:
                            kcat = kcat_sec_value
                        logging.info(f"Predicted kcat for reaction {r.id} and molecule {m.name} and gene {g.name}: {kcat_sec_value}")
                if kcat!=0:
                    r.annotation['DLkcat'] = {}
                    r.annotation['DLkcat']['kcat'] = kcat
                    try:
                        kcat_mw = kcat * 3600 * 1000 / totalmass
                        r.annotation['DLkcat']['kcat_mw'] = kcat_mw
                    except ZeroDivisionError:
                        pass
    '''




class KcatPrediction(nn.Module):
    """Taken directly from ECMpy
    """
    def __init__(self, device, n_fingerprint, n_word, dim, layer_gnn, window, layer_cnn, layer_output):
        super(KcatPrediction, self).__init__()
        self.embed_fingerprint = nn.Embedding(n_fingerprint, dim)
        self.embed_word = nn.Embedding(n_word, dim)
        self.W_gnn = nn.ModuleList([nn.Linear(dim, dim)
                                    for _ in range(layer_gnn)])
        self.W_cnn = nn.ModuleList([nn.Conv2d(
                     in_channels=1, out_channels=1, kernel_size=2*window+1,
                     stride=1, padding=window) for _ in range(layer_cnn)])
        self.W_attention = nn.Linear(dim, dim)
        self.W_out = nn.ModuleList([nn.Linear(2*dim, 2*dim)
                                    for _ in range(layer_output)])
        # self.W_interaction = nn.Linear(2*dim, 2)
        self.W_interaction = nn.Linear(2*dim, 1)

        self.device = device
        self.dim = dim
        self.layer_gnn = layer_gnn
        self.window = window
        self.layer_cnn = layer_cnn
        self.layer_output = layer_output

    def gnn(self, xs, A, layer):
        for i in range(layer):
            hs = torch.relu(self.W_gnn[i](xs))
            xs = xs + torch.matmul(A, hs)
        # return torch.unsqueeze(torch.sum(xs, 0), 0)
        return torch.unsqueeze(torch.mean(xs, 0), 0)

    def attention_cnn(self, x, xs, layer):
        """The attention mechanism is applied to the last layer of CNN."""

        xs = torch.unsqueeze(torch.unsqueeze(xs, 0), 0)
        for i in range(layer):
            xs = torch.relu(self.W_cnn[i](xs))
        xs = torch.squeeze(torch.squeeze(xs, 0), 0)

        h = torch.relu(self.W_attention(x))
        hs = torch.relu(self.W_attention(xs))
        weights = torch.tanh(F.linear(h, hs))
        ys = torch.t(weights) * hs

        # return torch.unsqueeze(torch.sum(ys, 0), 0)
        return torch.unsqueeze(torch.mean(ys, 0), 0)

    def forward(self, inputs):

        fingerprints, adjacency, words = inputs

        layer_gnn = 3
        layer_cnn = 3
        layer_output = 3

        """Compound vector with GNN."""
        fingerprint_vectors = self.embed_fingerprint(fingerprints)
        compound_vector = self.gnn(fingerprint_vectors, adjacency, layer_gnn)

        """Protein vector with attention-CNN."""
        word_vectors = self.embed_word(words)
        protein_vector = self.attention_cnn(compound_vector,
                                            word_vectors, layer_cnn)

        """Concatenate the above two vectors and output the interaction."""
        cat_vector = torch.cat((compound_vector, protein_vector), 1)
        for j in range(layer_output):
            cat_vector = torch.relu(self.W_out[j](cat_vector))
        interaction = self.W_interaction(cat_vector)
        # logging.debug(interaction)

        return interaction

    def __call__(self, data, train=True):

        inputs, correct_interaction = data[:-1], data[-1]
        predicted_interaction = self.forward(inputs)
        logging.debug(predicted_interaction)

        if train:
            loss = F.mse_loss(predicted_interaction, correct_interaction)
            return loss
        else:
            correct_values = correct_interaction.to('cpu').data.numpy()
            predicted_values = predicted_interaction.to('cpu').data.numpy()[0]
            # correct_values = np.concatenate(correct_values)
            # predicted_values = np.concatenate(predicted_values)
            # ys = F.softmax(predicted_interaction, 1).to('cpu').data.numpy()
            # predicted_values = list(map(lambda x: np.argmax(x), ys))
            logging.info(correct_values)
            logging.info(predicted_values)
            # predicted_scores = list(map(lambda x: x[1], ys))
            return correct_values, predicted_values
