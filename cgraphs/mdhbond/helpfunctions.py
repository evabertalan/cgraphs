#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#    Author: Malte Siemers, Freie Universität Berlin
#
#    If you use this software or anything it produces for work to be published,
#    please cite:
#
#    Malte Siemers, Michalis Lazaratos, Konstantina Karathanou,
#    Federico Guerra, Leonid Brown, and Ana-Nicoleta Bondar.
#    Bridge: A graph-based algorithm to analyze dynamic H-bond networks
#    in membrane proteins, Journal of Chemical Theory and Computation, 2019.

import numpy as np
import MDAnalysis
from collections import defaultdict
import networkx as nx
import itertools

r_covalent = defaultdict(lambda: 1.5, N=1.31, O=1.31, P=1.58, S=1.55)
donor_names_global = {'OH2', 'OW', 'NE', 'NH1', 'NH2', 'ND2', 'SG', 'NE2', 'ND1', 'NZ', 'OG', 'OG1', 'NE1', 'OH', 'OE1', 'OE2', 'N16', 'OD1', 'OD2'}
acceptor_names_global = {'OH2', 'OW', 'OD1', 'OD2', 'SG', 'OE1', 'OE2', 'ND1', 'NE2', 'SD', 'OG', 'OG1', 'OH'}
aa_three2one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'TIP3': 'O', 'LYR':'RET', 'HSE':'H', 'HSD':'H', 'HSP':'H'}
water_definition = '(resname TIP3 and name OH2) or (resname HOH and name O)'

class EmptyGroup:
    def __init__(self):
        self.positions = []

    def __len__(self):
        return 0


class Selection:
    def __init__(self, mda_residue, donor_names, acceptor_names, add_donors_without_hydrogen):

        self.donors = []
        self.donor_info = []
        self.acceptors = []
        self.acceptor_info = []
        self.hydrogens = []
        self.hydrogen_info = []
        self.donor2hydrogens = []

        if mda_residue:
            for atom in mda_residue.atoms:
                if (atom.name in donor_names)  or ((len(atom.name) > 1) and atom.name[0] in ['N', 'O']):
                    if add_donors_without_hydrogen:
                        self.donors.append(atom)
                        self.donor_info.append(
                            str(atom.segid) + '-' + atom.resname + '-' + str(atom.resid) + '-' + atom.name)
                    else:
                        c = 0
                        for otheratom in atom.residue.atoms:
                            if (otheratom.name[0] == 'H' or otheratom.name[:2] in ['1H', '2H', '3H'] or otheratom.type == 'H') and np.sqrt(((atom.position - otheratom.position) ** 2).sum()) < r_covalent[atom.name[0]]:
                                self.hydrogens.append(otheratom)
                                self.hydrogen_info.append(str(otheratom.index))
                                c += 1
                        if c>0:
                            self.donor2hydrogens.append(range(len(self.hydrogens) - c, len(self.hydrogens)))
                            self.donors.append(atom)
                            self.donor_info.append(
                                str(atom.segid) + '-' + atom.resname + '-' + str(atom.resid) + '-' + atom.name)
                if (atom.name in acceptor_names) or ((len(atom.name) > 1) and atom.name[0] in ['N', 'O']):
                    self.acceptors.append(atom)
                    self.acceptor_info.append(
                        str(atom.segid) + '-' + atom.resname + '-' + str(atom.resid) + '-' + atom.name)

            if self.donors: self.donors = MDAnalysis.core.groups.AtomGroup(self.donors)
            if self.acceptors: self.acceptors = MDAnalysis.core.groups.AtomGroup(self.acceptors)
            if self.hydrogens: self.hydrogens = MDAnalysis.core.groups.AtomGroup(self.hydrogens)
            self.count = len(self.donors) + len(self.acceptors)
            if not self.count>0: raise AssertionError('neither donors nor acceptors in the selection')


def angle(p1, p2, p3):
    """angle(p1,p2,p3) takes three numpy arrays of shape (N,3) as arguments.
    p1, p2 and p3 have to represent coordinates of the following form:
        p1
       a  \
      -----p2---p3
    The returned array of shape (N,1) contains degree angles, representing
    positive angles a between the lines p2--p3 and p1--p2 as indicated in the
    above scheme.
    """
    v1s = p2 - p1
    v2s = p3 - p2

    dot_v1_v2 = np.einsum('ij,ij->i', v1s, v2s)
    dot_v1_v1 = np.einsum('ij,ij->i', v1s, v1s)
    dot_v2_v2 = np.einsum('ij,ij->i', v2s, v2s)
    cos_arg = dot_v1_v2/(np.sqrt(dot_v1_v1)*np.sqrt(dot_v2_v2))
    cos_arg[cos_arg<-1.0] = -1.0
    return np.rad2deg(np.arccos(cos_arg))


def MDA_info_list(group, detailed_info=False, special_naming=[]):
    result = []
    for atom in group:
        sid = str(atom.segid)
        rname = str(atom.resname)
        rid = str(atom.resid)
        if rname in special_naming:
            rname = str(atom.name)
        if detailed_info: result.append(sid + '-' + rname + '-' + rid + '-' + atom.name)
        else: result.append(sid + '-' + rname + '-' + rid)
    return result


def dict2graph(bonds, residuewise=True):
    g = nx.Graph()
    temp = []
    for bond in bonds:
        node_a, node_b = bond.split(':')
        if residuewise: cut = 3
        else: cut = 4
        node_a = '-'.join(node_a.split('-')[:cut])
        node_b = '-'.join(node_b.split('-')[:cut])
        if int(node_a.split('-')[2]) > int(node_b.split('-')[2]):
            node_a, node_b = node_b, node_a
        if node_a != node_b:
            temp.append((node_a, node_b))
    g.add_edges_from(temp)
    return g


def check_angle(atoms_in_distance, heavy2hydrogen, local_coordinates, hydrogen_coordinates, cut_angle):
    pairs = np.asarray(atoms_in_distance)
    angle_check_index = []
    a_index = []
    b_index = []
    hydrogen_index = []
    for i, (atom_a, atom_b) in enumerate(atoms_in_distance):
        temp_hydrogen = list(heavy2hydrogen[atom_a]) + list(heavy2hydrogen[atom_b])
        hydrogen_index += temp_hydrogen
        nb_hydrogen = len(temp_hydrogen)
        angle_check_index += [i]*nb_hydrogen
        a_index += [atom_a]*nb_hydrogen
        b_index += [atom_b]*nb_hydrogen
    temp_b = local_coordinates[b_index]
    temp_h = hydrogen_coordinates[hydrogen_index]
    temp_a = local_coordinates[a_index]

    angles = angle(temp_b, temp_h, temp_a)
    angle_check = angles <= cut_angle
    angle_check_index = np.array(angle_check_index)
    bond_index = np.asarray(angle_check_index[angle_check.flatten()], dtype=np.int)
    hbond_pairs = pairs[bond_index]
    return hbond_pairs


def check_angle_water(atoms_in_distance, oxygen_coordinates, hydrogen_coordinates, cut_angle):
    pairs = np.asarray(atoms_in_distance)
    a_index = np.repeat(pairs[:,0], 4)
    b_index = np.repeat(pairs[:,1], 4)
    hydrogen_index = np.zeros(a_index.size, dtype=np.int)
    hydrogen_index[::4] = pairs[:,0] * 2
    hydrogen_index[1::4] = pairs[:,0] * 2 + 1
    hydrogen_index[2::4] = pairs[:,1] * 2
    hydrogen_index[3::4] = pairs[:,1] * 2 + 1
    temp_b = oxygen_coordinates[b_index]
    temp_h = hydrogen_coordinates[hydrogen_index]
    temp_a = oxygen_coordinates[a_index]
    angles = angle(temp_b, temp_h, temp_a)
    angle_check = angles <= cut_angle
    bond_index = angle_check.reshape(-1,4).any(axis=1)
    hbond_pairs = pairs[bond_index]
    return hbond_pairs


def intervals(timeseries):
    ts = np.array(timeseries)
    changes = np.nonzero(np.diff(ts))[0]+1
    changes = np.concatenate(([0],changes,[len(ts)]))
    intervals = [(changes[i],changes[i+1]-1) for i in range(len(changes)-1) if changes[i+1]-changes[i]>1]
    return intervals

def intervals_binary(timeseries):
    ts = np.array(timeseries, dtype=bool)
    if ts[0]:
        if ts[-1]: changes = np.concatenate(([0],np.nonzero(np.diff(ts))[0]+1,[ts.size]))
        else: changes = np.concatenate(([0],np.nonzero(np.diff(ts))[0]+1))
    else:
        if ts[-1]: changes = np.concatenate((np.nonzero(np.diff(ts))[0]+1,[ts.size]))
        else: changes = np.nonzero(np.diff(ts))[0]+1
    return changes.reshape((-1,2))

def block_analysis(dictionary, step_size, block_size, conv_cut):
    min_interval = range(0,int((len(dictionary[dictionary.keys()[0]])-block_size)/step_size), step_size)
    values = np.array(list(dictionary.values()), dtype=np.float)
    blockupancies = np.empty((values.shape[0],len(min_interval)))
    for ii,i in enumerate(min_interval):
        blockupancies[:,ii] = np.mean(values[:,i:i+block_size], axis=1)
    stds = np.std(blockupancies, axis=1)
    conv_index = stds < conv_cut
    result = dict(zip(np.array(dictionary.keys())[conv_index],np.array(dictionary.values())[conv_index]))
    return result


def filter_occupancy(dictionary, min_occupancy):
    for key in dictionary:
        dtype = dictionary[key].dtype
        break
    if dtype == np.int: values = np.array(np.array(list(dictionary.values())) < np.inf, dtype=np.float)
    else: values = np.array(list(dictionary.values()), dtype=np.float)
    filter_index = values.mean(axis=1) > min_occupancy
    keys = [key for key in dictionary]
    values = [dictionary[key] for key in dictionary]
    result = dict(zip(np.array(keys)[filter_index],np.array(values)[filter_index]))
    return result


def pca_2d_projection(pos3d):
    m = pos3d.mean(axis=0)
    pos3dm = pos3d - m
    S = pos3dm.T.dot(pos3dm)
    try:
        eig_val, eig_vec = np.linalg.eigh(S)
        eig_val, eig_vec = eig_val[::-1][:2], eig_vec.T[::-1][:2]
        return eig_vec.dot(pos3dm.T).T, eig_vec
    except:
        return np.array([1.0,1.0]), np.array([[1.0,0,0],[0,1.0,0]])


def predecessor_recursive(d,pred,start,stop):
    if d==0: return pred[start,stop]
    else: return pred[start,predecessor_recursive(d-1,pred,start,stop)]


def predecessor_recursive_1d(d,pred,start):
    if d==0: return pred[start]
    else: return pred[predecessor_recursive_1d(d-1,pred,start)]


def complete_subgraphs(edges_matrix):
    edges_matrix = edges_matrix | np.eye(edges_matrix.shape[0], dtype=bool)
    subgraphs = np.zeros(edges_matrix.shape, dtype=bool)
    for i in range(edges_matrix.shape[0]):
        indices = edges_matrix[i].nonzero()[0]
        if indices.size < 2: continue
        check_matrix = edges_matrix[indices][:, indices]
        subgraphs[i][indices[check_matrix.all(axis=0)]] = True
    subgraphs = np.unique(subgraphs, axis=0)
    subgraphs = subgraphs[np.logical_not([((subgraphs * p) == p).all(axis=1).sum() > 1 for p in subgraphs])]
    return subgraphs


def pairwise(iterable):
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


def connected_component_subgraph(graph, node):
    for component in nx.connected_component_subgraphs(graph):
        if node in component.nodes(): return component
    return None


def deconst_key(key, residuewise):
    a, b = key.split(':')
    if residuewise:
        segna, resna, resia = a.split('-')
        segnb, resnb, resib = b.split('-')

    else:
        segna, resna, resia, atomna = a.split('-')
        segnb, resnb, resib, atomnb = b.split('-')
    return segna, resna, int(resia), segnb, resnb, int(resib)


def remove_by_index(l, indices):
    for i,index in enumerate(indices):
        l.pop(index-i)


def points_in_slice(pt1, pt2, q):
    vec = pt2 - pt1
    return (np.dot(q - pt1, vec) >= 0) & (np.dot(q - pt2, vec) <= 0)


def points_in_cylinder(pt1, pt2, r, q):
    vec = pt2 - pt1
    const = r * np.linalg.norm(vec)
    return (np.dot(q - pt1, vec) >= 0) & (np.dot(q - pt2, vec) <= 0) & (np.linalg.norm(np.cross(q - pt1, vec), axis=1) <= const)

def msd_fft_broadcast(r):
  N=len(r)
  F = np.fft.fft(r, n=2*N, axis=0)
  PSD = F * F.conjugate()
  res = np.fft.ifft(PSD, axis=0)
  res= (res[:N]).real
  n=np.arange(0,N)[::-1]+1
  S2 = (res/n.reshape(-1,1,1)).sum(-1)
  D=np.square(r).sum(-1)
  D=np.vstack((D,np.zeros(D.shape[1])))
  Q=2*D.sum(0)
  S1=np.zeros((N,Q.size))
  for m in range(N):
      Q=Q-D[m-1]-D[N-m]
      S1[m]=Q/(N-m)
  return (S1-2*S2).mean(-1)

def msd_fft(r):
  N=len(r)
  F = np.fft.fft(r, n=2*N, axis=0)
  PSD = F * F.conjugate()
  res = np.fft.ifft(PSD, axis=0)
  res= (res[:N]).real
  n=np.arange(0,N)[::-1]+1
  S2 = (res/n.reshape(-1,1)).sum(-1)
  D=np.square(r).sum(-1)
  D=np.append(D,0)
  Q=2*D.sum(0)
  S1=np.zeros(N)
  for m in range(N):
      Q=Q-D[m-1]-D[N-m]
      S1[m]=Q/(N-m)
  return S1-2*S2

def invperm(p):
    q = np.empty_like(p)
    q[p] = np.arange(len(p))
    return q

def find_map(arr1, arr2):
    o1 = np.argsort(arr1)
    o2 = np.argsort(arr2)
    return o2[invperm(o1)]

def string_in_columns(s):
    lines = s.split('\n')
    if lines[-1] == '': lines = lines[:-1]
    elements = []
    for line in lines:
        elements.append(line.split(' '))
    rotated_elements = np.array(elements, dtype=np.str).T
    return_string = ''
    for rot_element in rotated_elements:
        return_string += ' '.join(rot_element)+'\n'
    return return_string
