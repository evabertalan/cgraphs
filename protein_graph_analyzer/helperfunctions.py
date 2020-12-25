import os
import logging
import Bio
import pickle
import numpy as np
from Bio import SeqIO
from Bio import pairwise2
from Bio.PDB import Selection
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

amino_d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'HSD':'H', 'HSE':'H'}


def get_pdb_files(folder):
    return [file for file in os.listdir(folder) if file.endswith('.pdb')]

def get_files(folder, endswith):
    return [file for file in os.listdir(folder) if file.endswith(endswith)]

def pickle_write_file(path, obj):
    with open(path, 'wb') as fp:
        pickle.dump(obj, fp)

def concatenate_arrays(arrays):
    concatenated = []
    for arr in arrays:  
        if arr.ndim > 1:
            for row in arr:
                concatenated.append(row)
        elif arr.size != 0:
            concatenated.append(arr)
    return np.array(concatenated) 

def load_pdb_model(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('model', pdb_file)
    return structure


def water_in_pdb(pdb_file):
    structure = load_pdb_model(pdb_file)
    residues = list(structure[0].get_residues())
    waters = [res for res in residues if res.get_id()[0] == 'W']
    return waters


def water_coordinates(pdb_file):
    waters = water_in_pdb(pdb_file)
    water_coord = [water['O'].get_coord() for water in waters]
    return np.array(water_coord)


def get_sequence(pdb_file):
    sequence = []
    structure = load_pdb_model(pdb_file)
    ppb = PPBuilder()
    for pp in ppb.build_peptides(structure):
        sequence.append(pp.get_sequence())
    seq = ''
    for s in sequence:
        seq += s
    return seq

def align_sequence(pdb_ref, pdb_move, threshold=0.75):
    ref_sequence = get_sequence(pdb_ref)
    move_sequence = get_sequence(pdb_move)
    alignments = pairwise2.align.globalxx(ref_sequence, move_sequence)
    
    best_score = 0
    best_i = 0
    for i, alignment in enumerate(alignments):
        if best_score <= alignment.score:
            best_score = alignment.score
            best_i = i

    best_alginment = alignments[best_i]
    if len(best_alginment.seqA) != len(best_alginment.seqB):
        print('Aligned sequences have different lenght') #TODO: change all print to log
        return None, None
    if best_alginment.score / len(alignments[best_i].seqA) <= threshold:
        print('Sequences of '+pdb_move+' differ more than the threshold value ('+str(threshold*100)+'%) from the reference structure')
        return None, None
    if best_alginment.score / len(alignments[best_i].seqB) <= threshold:
        print('Sequences of '+pdb_move+' differ more than the threshold value ('+str(threshold*100)+'%) from the reference structure')
        return None, None
    return best_alginment.seqA, best_alginment.seqB

def superimpose_aligned_atoms(seq_ref, pdb_ref, seq_move, pdb_move, file_name='', save=True):
    if file_name == '': file_name = pdb_move.split('/')[-1].split('.pdb')[0] 
    else: file_name = file_name.split('.pdb')[0]
    #TODO: maybe creae regex or parameter to filnave OR retihnik this filename conscept
    ref_atoms = []
    move_atoms = []
    ref_struct = load_pdb_model(pdb_ref)
    move_struct = load_pdb_model(pdb_move)
    
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
#     print(super_imposer.rms)
    if super_imposer.rms > 5:
        print('Automatic superimposition of '+file_name+' was not sucessful, please provide a pdb file superimposed to the reference structure')
        return
    logging.info('RMS value of superimposed '+file_name+'to the reference structure is '+str(super_imposer.rms))
    io = Bio.PDB.PDBIO()
    io.set_structure(move_struct)
    if save: io.save(file_name+'_superimposed.pdb')
    return move_struct
    
def calculate_pca_positions(coordinates):
        pca_positions = {}
        XY = [i[0:2] for i in coordinates.values()]
        pca = PCA(n_components=1)
        xy = pca.fit_transform(XY)

        for i, (key, value) in enumerate(coordinates.items()):
            pca_positions[key] = [xy[i][0], value[2]]
        return pca_positions
    
def create_plot(figsize=(10,15), title='', xlabel='', ylabel=''):
    fig, ax = plt.subplots(figsize=figsize)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_title(title, fontsize=20)
    ax.set_xlabel(xlabel, fontsize=19)
    ax.set_ylabel(ylabel, fontsize=19)
    ax.tick_params(axis='x', labelsize=17)
    ax.tick_params(axis='y', labelsize=17)
    return fig, ax