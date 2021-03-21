import os
import re
import logging
import pickle
import warnings
import Bio
from Bio import pairwise2
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio import BiopythonWarning
import numpy as np
import networkx as nx
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
warnings.simplefilter('ignore', BiopythonWarning)

#TODO support and test multiple protein chains chains = [chain for chain in model]

amino_d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'HSD':'H', 'HSE':'H', 'LYR':'X'}#TODO: remove ret

def create_logger(folder):
    logger = logging.getLogger('cgraph')
    if not len(logger.handlers):
        logger.setLevel(logging.INFO)
        fh = logging.FileHandler(folder+'cgraph_logs.log')
        with open(folder+'cgraph_logs.log', 'w'):
            pass
        fh.setLevel(logging.DEBUG)
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        fh_form = logging.Formatter('%(asctime)s - %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
        ch_form = logging.Formatter('%(levelname)s: %(message)s')
        fh.setFormatter(fh_form)
        ch.setFormatter(ch_form)
        logger.addHandler(fh)
        logger.addHandler(ch)
    return logger

def create_directory(directory):
    if not os.path.isdir(directory):
        os.makedirs(directory)
    return directory

def get_pdb_files(folder):
    return [file for file in os.listdir(folder) if file.endswith('.pdb')]

def get_files(folder, endswith):
    return [file for file in os.listdir(folder) if file.endswith(endswith)]

def pickle_write_file(path, obj):
    with open(path, 'wb') as fp:
        pickle.dump(obj, fp)

def get_node_name(node):
    chain, res, ind = node.split('-')
    #FIX water id issue from mdhbond --> issue from MDAnalysis
    # if res == 'HOH' and int(ind) > 10000: node = res+'-'+str(int(ind)-10000)
    # else: node = res+'-'+ind
    return res+'-'+ind

def concatenate_arrays(arrays):
    concatenated = []
    for arr in arrays:
        if arr.ndim > 1:
            for row in arr:
                concatenated.append(row)
        elif arr.size != 0:
            concatenated.append(arr)
    return np.array(concatenated)

def load_pdb_structure(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('model', pdb_file)
    return structure


def water_in_pdb(pdb_file):
    structure = load_pdb_structure(pdb_file)
    residues = list(structure[0].get_residues())
    waters = [res for res in residues if res.get_id()[0] == 'W']
    return waters


def water_coordinates(pdb_file):
    waters = water_in_pdb(pdb_file)
    water_coord = [water['O'].get_coord() for water in waters]
    return np.array(water_coord)


def get_sequence(pdb_file):
    sequence = []
    structure = load_pdb_structure(pdb_file)
    ppb = PPBuilder()
    for pp in ppb.build_peptides(structure):
        sequence.append(pp.get_sequence())
    seq = ''
    for s in sequence:
        seq += s
    return seq

def align_sequence(logger, pdb_ref, pdb_move, threshold=0.75):
    ref_sequence = get_sequence(pdb_ref)
    move_sequence = get_sequence(pdb_move)
    alignments = pairwise2.align.globalxx(ref_sequence, move_sequence)
    pdb_name = pdb_move.split('/')[-1]

    best_score = 0
    best_i = 0
    for i, alignment in enumerate(alignments):
        if best_score <= alignment.score:
            best_score = alignment.score
            best_i = i

    best_alginment = alignments[best_i]
    if len(best_alginment.seqA) != len(best_alginment.seqB):
        logger.warning('Aligned sequences have different lenght')
        logger.info('Thus '+pdb_name+' is excluded from further analysis.')
        return None, None
    if best_alginment.score / len(alignments[best_i].seqA) <= threshold:
        logger.warning('Sequences of '+pdb_name+' has lower sequence identity than the threshold value ('+str(threshold*100)+'%) compared to the reference structure')
        logger.info('Thus '+pdb_name+' is excluded from further analysis.')
        return None, None
    if best_alginment.score / len(alignments[best_i].seqB) <= threshold:
        logger.warning('Sequences of '+pdb_name+' has lower sequence identity than the threshold value ('+str(threshold*100)+'%) compared to the reference structure')
        logger.info('Thus '+pdb_name+' s excluded from further analysis.')
        return None, None
    return best_alginment.seqA, best_alginment.seqB

def superimpose_aligned_atoms(logger, seq_ref, pdb_ref, seq_move, pdb_move, save_file_to='', save=True):
    if save_file_to == '': save_file_to = pdb_move.split('/')[-1].split('.pdb')[0]
    else: save_file_to = save_file_to.split('.pdb')[0]
    #TODO: maybe creae regex or parameter to filnave OR retihnik this filename conscept
    pdb_name = pdb_move.split('/')[-1]
    ref_atoms = []
    move_atoms = []
    ref_struct = load_pdb_structure(pdb_ref)
    move_struct = load_pdb_structure(pdb_move)

    #TODO
#     assert len(ref_struct) == len(move_struct) == 1, 'The reference structure and '+pdb_move+'have more than one protein chain in the provided pdb file.'
    ref_residues = list(ref_struct[0].get_residues())
    move_residues = list(move_struct[0].get_residues())

    for i, r in enumerate(seq_ref):
        for j, m in enumerate(seq_move):
            if ((r == m) and r in amino_d.values() and move_residues[j].get_resname() in amino_d.keys() and ref_residues[i].get_resname() in amino_d.keys() and ref_residues[i].get_id()[1] == move_residues[j].get_id()[1]):
                ref_atoms.append(ref_residues[i]['CA'])
                move_atoms.append(move_residues[j]['CA'])

    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, move_atoms)
    all_atoms = []
    for model in move_struct:
        all_atoms += list(model.get_atoms())
    super_imposer.apply(all_atoms)
    logger.debug('Superimposition RMS value of '+pdb_name+' to the reference structure is: '+str(super_imposer.rms))
    if super_imposer.rms > 5:
        logger.warning('Automatic superimposition of '+pdb_name+' was not sucessful, please provide a pdb file superimposed to the reference structure. This structure is excluded from further analysis.')
        return
    io = Bio.PDB.PDBIO()
    io.set_structure(move_struct)
    if save: io.save(save_file_to+'_superimposed.pdb')
    logger.debug('Superimposed file is saved as: '+str(save_file_to+'_superimposed.pdb'))
    return move_struct

def get_connected_components(graph):
    return list(nx.connected_components(graph))

def get_water_coordinates(protein_chain, res_index):
    #FIX water id issue from mdhbond --> issue from MDAnalysis
    if int(res_index) > 10000: res_index = int(res_index) - 10000
    return protein_chain[('W', int(res_index), ' ')]['O'].get_coord()

def calculate_connected_compontents_coordinates(connected_components, protein_chain):
    all_chains = []
    for connected_chain in connected_components:
        chain_details = []
        for res_name in list(connected_chain):
            res_index = int(res_name.split('-')[-1])
            if re.search('HOH', res_name): coords = get_water_coordinates(protein_chain, res_index)
            else: coords = protein_chain[res_index]['CA'].get_coord()
            chain_details.append(np.array([res_name, coords], dtype='object'))
        all_chains.append(chain_details)

    return [c for c in sorted(all_chains, key=len, reverse=True)]

def calculate_pca_positions(coordinates):
        pca_positions = {}
        XY = [i[0:2] for i in coordinates.values()]
        pca = PCA(n_components=1)
        xy = pca.fit_transform(XY)

        for i, (key, value) in enumerate(coordinates.items()):
            pca_positions[key] = [xy[i][0], value[2]]
        return pca_positions

def check_projection_sign(projection, reference):
    for i in range(len(reference.keys())):
        _k = list(reference.keys())[i]
        if _k in projection.keys() and reference[_k][0] > 5:
            if np.sign(reference[_k][0]) != np.sign(projection[_k][0]):
                t_projection = {}
                for k, v in projection.items():
                    t_projection.update({k: [v[0]*-1, v[1]]})
                projection = t_projection
                return projection
    return projection

#TODO set back plot size from git
def create_plot(figsize=(15,16), title='', xlabel='', ylabel=''):
    fig, ax = plt.subplots(figsize=figsize)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_title(title, fontsize=20)
    ax.set_xlabel(xlabel, fontsize=36)
    ax.set_ylabel(ylabel, fontsize=36)
    ax.tick_params(axis='x', labelsize=33)
    ax.tick_params(axis='y', labelsize=33)
    return fig, ax
