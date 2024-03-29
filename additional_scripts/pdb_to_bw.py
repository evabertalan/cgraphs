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
    uniprot_url = 'https://www.uniprot.org/uniprot/?query=id:'+str(uniprot_id)+'&format=tab&columns=entry name'
    response = requests.get(uniprot_url)
    entry_name = re.split(r'\n+', response.text)[1]
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
    protein_struct.atoms.write(target_folder+str(pdb_id)+'_bw.pdb')
    

def renumber_pdb_to_bw(pdb_id, pdb_file, target_folder):
    uniprot_id = get_uniport_from_pdb(pdb_id)
    seq_details = get_detailed_sequence_numbering(uniprot_id, save_to=target_folder, pdb_id=pdb_id)
    write_renumbered_pdb_file(pdb_file, pdb_id, seq_details, target_folder)