import os
import logging
import Bio
from Bio import SeqIO
from Bio import pairwise2
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder

amino_d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def superimpose_structures(self):
  pass

def get_pdb_files(folder):
    return [file for file in os.listdir(folder) if file.endswith('.pdb')]

def water_in_pdb(pdb_file):
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure('model', pdb_file)
    residues = list(struct[0].get_residues())
    waters = [res for res in residues if res.get_id()[0] == 'W']
    return waters

def get_sequence(pdb_file):
    sequence = []
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('model', pdb_file)
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
    assert len(best_alginment.seqA) == len(best_alginment.seqB), 'Aligned sequences have different lenght'
    assert best_alginment.score / len(alignments[best_i].seqA) >= threshold, 'Sequences of '+pdb_move+' differ more than the threshold value ('+str(threshold*100)+'%) from the reference structure'
    assert best_alginment.score / len(alignments[best_i].seqB) >= threshold, 'Sequences of '+pdb_move+' differ more than the threshold value ('+str(threshold*100)+'%) from the reference structure'

    return best_alginment.seqA, best_alginment.seqB

def superimpose_aligned_atoms(seq_ref, pdb_ref, seq_move, pdb_move, file_name=''):
    if file_name == '': file_name = pdb_move.split('/')[-1][:-4] #TODO: maybe creae regex or parameter to filnave
    parser = PDBParser(QUIET=True)
    ref_atoms = []
    move_atoms = []
    ref_struct = parser.get_structure('model', pdb_ref)
    move_struct = parser.get_structure('model', pdb_move)
    
    ref_residues = list(ref_struct[0].get_residues())
    move_residues = list(move_struct[0].get_residues())
    for i,(r, m) in enumerate(zip(seq_ref, seq_move)):
        if (r == m) and r in amino_d.values() and move_residues[i].get_resname() in amino_d.keys():
            ref_atoms.append(ref_residues[i]['CA'])
            move_atoms.append(move_residues[i]['CA'])
    
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, move_atoms)
    super_imposer.apply(move_struct[0].get_atoms())
    logging.info('RMS value of superimposed '+pdb_move.split('/')[-1][:-4]+'to the reference structure is '+str(super_imposer.rms))
    io = Bio.PDB.PDBIO()
    io.set_structure(move_struct)
    io.save(file_name+'_superimposed.pdb')