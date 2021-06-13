import os
import re
import logging
import json
import pickle
import warnings
import MDAnalysis as _mda
import MDAnalysis.analysis.align as md_align
import MDAnalysis.analysis.rms  as rms
from Bio.SVDSuperimposer import SVDSuperimposer

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
warnings.filterwarnings('ignore')
#TODO support and test multiple protein chains chains = [chain for chain in model]

amino_d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'HSD':'H', 'HSE':'H'}
water_def = '(resname TIP3 and name OH2) or (resname HOH and name O)'


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

def pickle_load_file(path):
    with open(path, 'rb') as fp:
        obj = pickle.load(fp)
    return obj

def json_write_file(path, obj):
    with open(path, 'w', encoding='utf-8') as fp:
        json.dump(obj, fp)

def get_node_name(node):
    return node

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
    u = _mda.Universe(pdb_file)
    return u

def water_in_pdb(pdb_file):
    structure = load_pdb_structure(pdb_file)
    waters = structure.select_atoms(water_def)
    return waters


def water_coordinates(pdb_file):
    waters = water_in_pdb(pdb_file)
    water_coord = waters.positions
    return np.array(water_coord)


def get_sequence(pdb_file):
    structure = load_pdb_structure(pdb_file)
    protein = structure.select_atoms('protein and name CA')
    seq = ''
    for res in protein:
        seq += amino_d[res.resname]
    print(seq)
    return seq

def align_sequence(logger, pdb_ref, pdb_move, threshold=0.75):
    ref_sequence = get_sequence(pdb_ref)
    move_sequence = get_sequence(pdb_move)
    alignments = pairwise2.align.globalxx(ref_sequence, move_sequence)
    pdb_name = pdb_move.split('/')[-1]

    best_score = 0
    best_i = 0
    for i, alignment in enumerate(alignments):
        print(alignment)
        print(alignment.score)
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
    if save_file_to == '': save_file_to = retrieve_pdb_code(pdb_move, '.pdb')
    else: save_file_to = save_file_to.split('.pdb')[0]
    #TODO: maybe creae regex or parameter to filnave OR retihnik this filename conscept
    pdb_name = pdb_move.split('/')[-1]
    seq_move = seq_move.replace('-', '')
    seq_ref = seq_ref.replace('-', '')

    ref_struct = load_pdb_structure(pdb_ref)
    move_struct = load_pdb_structure(pdb_move)
    move_atoms_pos=[]
    ref_atoms_pos =[]
    ref_atoms = ref_struct.select_atoms('protein and name CA')
    move_atoms = move_struct.select_atoms('protein and name CA')
    unique_seg_move, unique_seg_ref = np.unique(move_atoms.segids), np.unique(ref_atoms.segids)
    if len(unique_seg_move)==1 and len(unique_seg_ref)==1 and unique_seg_move != unique_seg_ref:
        logger.warning('Chains must have the same ID. To compare the conserved graph of multiple structures, same segments has to have the same chain ID.')
        logger.info('Chins of '+pdb_name+' has different chain ID than the reference structure. Thus excluded from further analysis.')
        return None

    print('len(ref_atoms)',len(ref_atoms))
    print('len(seq_ref)',len(seq_ref))
    print('len(ref_atoms)',len(ref_atoms))
    print('len(seq_move)',len(seq_move))
    for i, r in enumerate(seq_ref):
        for j, m in enumerate(seq_move):
            if ((r == m) and r in amino_d.values() and move_atoms[j].resname in amino_d.keys() and ref_atoms[i].resname in amino_d.keys() and ref_atoms[i].segid == move_atoms[j].segid):
                ref_atoms_pos.append(ref_atoms[i].position)
                move_atoms_pos.append(move_atoms[j].position)

    move_atoms_pos = np.array(move_atoms_pos,  dtype='float64').reshape(-1, 3)
    ref_atoms_pos = np.array(ref_atoms_pos,  dtype='float64').reshape(-1, 3)

    sup = SVDSuperimposer()
    sup.set(ref_atoms_pos, move_atoms_pos)
    sup.run()
    rot, tran = sup.get_rotran()

    rot = rot.astype('f')
    tran = tran.astype('f')
    move_struct.atoms.positions = np.dot(move_struct.atoms.positions, rot) + tran
    move_struct.atoms.write(save_file_to+'_superimposed.pdb')

    logger.info('Superimposition RMS value of '+pdb_name+' to the reference structure is: '+str(sup.get_rms()))
    logger.debug('Superimposed file is saved as: '+str(save_file_to+'_superimposed.pdb'))
    return move_struct

def get_connected_components(graph):
    return list(nx.connected_components(graph))

def get_water_coordinates(protein_chain, res_index):
    #FIX water id issue from mdhbond --> issue from MDAnalysis
    # if int(res_index) > 10000: res_index = int(res_index) - 10000
    # return protein_chain[('W', int(res_index), ' ')]['O'].get_coord()
    sel = protein_chain.select_atoms(water_def+' and resid '+ res_index)
    if len(sel.positions):
        return sel.positions[0]
    else:
        print('Water '+res_index+ ' not found in the PDB file. INFO: https://github.com/evabertalan/cgraph/blob/main/README.md')
        return None

def calculate_connected_compontents_coordinates(connected_components, struct_object, option='pdb'):
    all_chains = []
    for connected_chain in connected_components:
        chain_details = []
        for node in list(connected_chain):
            chain_id, res_name, res_id  = node.split('-')[0], node.split('-')[1], node.split('-')[2]
            if option  == 'pdb':
                chain = struct_object.select_atoms('segid '+ chain_id)
                if res_name in ['HOH', 'TIP3']: coords = get_water_coordinates(chain, res_id)
                else: coords = chain.select_atoms('name CA and resid '+ res_id).positions[0]

            else: coords = struct_object.select_atoms('resid '+ res_id).positions[0]
            if coords is not None:
                chain_details.append(np.array([node, coords], dtype='object'))
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

def is_conserved_edge(conserved_edges, e0, e1):
    conserved_edge = (len(np.where((conserved_edges == [e0, e1]).all(axis=1))[0]) != 0 or len(np.where((conserved_edges == [e1, e0]).all(axis=1))[0]) != 0)
    conserved_edge_with_water = False
    for edge in conserved_edges:
        if (e0 in edge and e1.split('-')[1] in ['HOH', 'TIP3'] and (edge[0].startswith('X-w') or edge[1].startswith('X-w'))) or (e1 in edge and e0.split('-')[1] in ['HOH', 'TIP3'] and (edge[0].startswith('X-w') or edge[1].startswith('X-w'))):
            conserved_edge_with_water = True
    return conserved_edge or conserved_edge_with_water

def retrieve_pdb_code(file_path, split_by):
    """ split_by e.g.: '.pdb' """
    return file_path.split('/')[-1].split(split_by)[0]

def _average_timeseries(hbond_dict):
    return {key: np.mean(hbond_dict[key]) for key in hbond_dict}

def get_edge_params(wba, edges):
    average_water_per_wire = wba.compute_average_water_per_wire()
    occupancy_per_wire = _average_timeseries(wba.filtered_results)

    keys = []
    waters = []
    occ_per_wire = []
    for edge in edges:
        key = str(edge[0])+':'+str(edge[1])
        keys.append(key)

        if key in average_water_per_wire:
            waters.append(average_water_per_wire[key])
        else:
            key = str(edge[1])+':'+str(edge[0])
            waters.append(average_water_per_wire[key])

        if key in occupancy_per_wire:
            occ_per_wire.append(occupancy_per_wire[key])
        else:
            key = str(edge[1])+':'+str(edge[0])
            occ_per_wire.append(occupancy_per_wire[key])

    return waters, occ_per_wire, keys

def edge_info(wba, edges):
    waters, occ_per_wire, keys = get_edge_params(wba, edges)
    edge_info = {}

    for w, o, k in zip(waters, occ_per_wire, keys):
        edge_info.update({k: {'waters': np.round(w,1), 'occupancy': o }})
    return edge_info

def write_text_file(file_path, text_content, logger=None):
    if file_path.endswith('.txt') and isinstance(text_content, list):
        f = open(file_path, 'w')
        f.writelines(text_content)
        f.close()
    elif logger:
        logger.error('The file name has to end to .txt and the text content has to be passed as a list.')


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
