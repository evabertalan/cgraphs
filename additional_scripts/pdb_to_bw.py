import requests
import re
import json
import copy
import MDAnalysis as _mda

amino_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'HSD':'H', 'HSE':'H'}

def get_uniport_from_pdb(pdb_id, entity_id='1'):
    pdb_url =' https://data.rcsb.org/rest/v1/core/polymer_entity/'+str(pdb_id)+'/'+str(entity_id)
    response = requests.get(pdb_url)
    json_pdb_res = response.json()
    for i in json_pdb_res['rcsb_polymer_entity_align']:
        if 'reference_database_accession' in i.keys():
            uniprot_id = i['reference_database_accession']
            return uniprot_id
    print('UNIPROT ID not found. Please insert it manually to the next function call.')
    return

def get_detailed_sequence_numbering(uniprot_id, save_to='', pdb_id=None,):
    uniprot_url = f'https://rest.uniprot.org/uniprotkb/{str(uniprot_id)}/?format=tsv'
    response = requests.get(uniprot_url)
    entry_name = response.text.split('\n')[1].split('\t')[1]
    entry_url = 'https://gpcrdb.org/services/residues/extended/'+entry_name.lower()+'/'
    response = requests.get(entry_url)
    jsonRes = response.json()
    
    if save_to:
        code = str(pdb_id) if pdb_id else str(uniprot_id)
        with open(str(save_to)+'/'+code+'.json', 'w') as outfile:
            json.dump(jsonRes, outfile)
    return jsonRes

def write_renumbered_pdb_file(pdb_file, pdb_id, seq_details, target_folder):
    protein_struct = _mda.Universe(pdb_file)
    og_struct = _mda.Universe(pdb_file)
    for i, res  in enumerate(og_struct.atoms):
        if res.resid-1 < len(seq_details):
            generic_numbers = seq_details[res.resid-1]['alternative_generic_numbers']
            if res.resname in amino_dict.keys() and seq_details[res.resid-1]['amino_acid'] == amino_dict[res.resname] and len(generic_numbers):
                if generic_numbers[0]['scheme'] == 'BW':
                    bw_number = int(generic_numbers[0]['label'].replace('.', ''))
                    protein_struct.atoms[i].residue.resid = bw_number
                    protein_struct.atoms[i].residue.resname = 'BWX'
    protein_struct.atoms.write(f'{target_folder}/{str(pdb_id)}_bw.pdb')
    

def renumber_pdb_to_bw(pdb_id, pdb_file, target_folder):
    uniprot_id = get_uniport_from_pdb(pdb_id)
    seq_details = get_detailed_sequence_numbering(uniprot_id, save_to=target_folder, pdb_id=pdb_id)
    write_renumbered_pdb_file(pdb_file, pdb_id, seq_details, target_folder)


def assign_generic_numbers(pdb_code, pdb_file):
    post_url = 'https://services/structure/assign_generic_numbers'
    response = requests.post(post_url, data=pdb_file)
    jsonRes = response.json()
    print(jsonRes)

def write_renumbered_pdb_file_from_generic_PDB(pdb_id, generic_pdb, target_folder):
    protein_struct = _mda.Universe(generic_pdb)
    og_struct = _mda.Universe(generic_pdb)
    bw_numbers = {}
    for i, res  in enumerate(og_struct.residues):
        bw_numbers.update({f'{res.resname}-{res.resid}': round(res.atoms[0].tempfactor*100)})

    for i, atom  in enumerate(og_struct.atoms):
        bw_number = bw_numbers[f'{atom.resname}-{atom.resid}']
        if atom.resname in amino_dict.keys() and bw_number < 900:
            protein_struct.atoms[i].residue.resid = bw_number
            protein_struct.atoms[i].residue.resname = 'BWX'
    protein_struct.atoms.write(f'{target_folder}/{str(pdb_id)}_bw.pdb')

