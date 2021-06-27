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

from . import helpfunctions as _hf
from .basic import BasicFunctionality
import numpy as _np
import networkx as _nx
import MDAnalysis as _MDAnalysis
from collections import OrderedDict as _odict
from itertools import combinations
from matplotlib.ticker import MaxNLocator
import matplotlib
#matplotlib.use('TKAgg', warn=False)
import matplotlib.pyplot as _plt

class NetworkAnalysis(BasicFunctionality):

    def __init__(self, selection=None, structure=None, trajectories=None, distance=3.5, cut_angle=60.,
                 start=None, stop=None, step=1, additional_donors=[],
                 additional_acceptors=[], exclude_donors=[], exclude_acceptors=[],
                 ions=[], check_angle=True, residuewise=True, add_donors_without_hydrogen=False, restore_filename=None):

        super(NetworkAnalysis, self).__init__(selection=selection, structure=structure, trajectories=trajectories,
             start=start, stop=stop, step=step, ions=ions, restore_filename=restore_filename)
        if restore_filename != None: return
        self.donor_names = _hf.donor_names_global.union(additional_donors)-set(exclude_donors)
        self.acceptor_names = _hf.acceptor_names_global.union(additional_acceptors)-set(exclude_acceptors)
        self.check_angle = check_angle
        self.distance = distance
        self.cut_angle = cut_angle
        self._add_donors_without_hydrogen = add_donors_without_hydrogen
        sorted_selection = _hf.Selection(self._mda_selection, self.donor_names, self.acceptor_names, add_donors_without_hydrogen)
        if not sorted_selection.donors: da_selection = sorted_selection.acceptors
        elif not sorted_selection.acceptors: da_selection = sorted_selection.donors
        else: da_selection = _MDAnalysis.core.groups.AtomGroup(sorted_selection.donors + sorted_selection.acceptors)
        self._da_selection = da_selection
        if sorted_selection.donors: self._donors = _MDAnalysis.core.groups.AtomGroup(sorted_selection.donors)
        else: self._donors = _hf.EmptyGroup()
        self._nb_donors = len(self._donors)
        if sorted_selection.acceptors: self._acceptors = _MDAnalysis.core.groups.AtomGroup(sorted_selection.acceptors)
        else: self._acceptors = _hf.EmptyGroup()
        self._nb_acceptors = len(self._acceptors)
        da_ids = _hf.MDA_info_list(da_selection, detailed_info=not residuewise)
        da_ids_atomwise = _hf.MDA_info_list(da_selection, detailed_info=True)
        self._first_water_id = len(da_selection)
        self._first_water_hydrogen_id = len(sorted_selection.hydrogens)
        if not add_donors_without_hydrogen:
            try:
                water_hydrogen = _MDAnalysis.core.groups.AtomGroup(self._water[0].residue.atoms[1:])
                for l in self._water[1:]:
                    water_hydrogen += l.residue.atoms[1:]
            except:
                water_hydrogen = []
        else:
            water_hydrogen = []
        #water_hydrogen = [h for l in self._water for h in l.residue.atoms[1:]]
        if not sorted_selection.hydrogens and not water_hydrogen:
            if check_angle: raise AssertionError('There are no possible hbond donors in the selection and no water. Since check_angle is True, hydrogen is needed for the calculations!')
            else: hydrogen = _hf.EmptyGroup()
        elif not sorted_selection.hydrogens: hydrogen = _MDAnalysis.core.groups.AtomGroup(water_hydrogen)
        elif not water_hydrogen: hydrogen = sorted_selection.hydrogens
        else: hydrogen = sorted_selection.hydrogens + _MDAnalysis.core.groups.AtomGroup(water_hydrogen)
        self._hydrogen = hydrogen
        self.heavy2hydrogen = sorted_selection.donor2hydrogens + [[] for i in sorted_selection.acceptors] + [[self._first_water_hydrogen_id+i, self._first_water_hydrogen_id+i+1] for i in range(0, len(water_hydrogen), 2)]
        self._all_ids = da_ids+self._water_ids
        self._all_ids_atomwise = da_ids_atomwise + self._water_ids_atomwise
        self._resids = _np.array([int(ids.split('-')[2]) for ids in self._all_ids])
        self.initial_graph = _nx.Graph()
        self.filtered_graph = self.initial_graph
        self.joint_occupancy_series = None
        self.joint_occupancy_frames = None
        self.residuewise=residuewise
        self._exclude_backbone_backbone = True
        self.add_missing_residues = 0
        self._connection_position = None

    def filter_occupancy(self, min_occupancy, use_filtered=True):
        if use_filtered: results = self.filtered_results
        else: results = self.initial_results
        if len(results) == 0: raise AssertionError('nothing to filter!')
        filtered_result = {key:results[key] for key in results if _np.mean(results[key])>min_occupancy}
        if self.applied_filters['occupancy'] != None: self.applied_filters['occupancy']=max(self.applied_filters['occupancy'], min_occupancy)
        else: self.applied_filters['occupancy'] = min_occupancy
        self.filtered_results = filtered_result
        self._generate_filtered_graph_from_filtered_results()

    def filter_connected_component(self, root, atomwise_whole_residue=False, use_filtered=True):
        if use_filtered: graph = self.filtered_graph
        else: graph = self.initial_graph
        if len(graph.nodes()) == 0: raise AssertionError('nothing to filter!')
        if (not self.residuewise) and atomwise_whole_residue:
            start_points = []
            residue = root.split('-')[:3]
            for node in graph.nodes():
                if node.split('-')[:3] == residue: start_points.append(node)
        if root not in graph.nodes(): raise AssertionError('The root node is not in the current graph')
        if (not self.residuewise) and atomwise_whole_residue:
            components = []
            for component in _nx.connected_component_subgraphs(graph):
                for start_point in start_points:
                    if start_point in component.nodes(): components.append(component)
            component = _nx.compose_all(components)
        else:
            for component in _nx.connected_component_subgraphs(graph):
                if root in component.nodes(): break
        self.applied_filters['connected_component'] = root
        self.filtered_graph = component
        self._generate_filtered_results_from_filtered_graph()

    def filter_all_paths(self, start, goal, max_len=_np.inf, only_shortest=True, use_filtered=True):
        if use_filtered: graph = self.filtered_graph
        else: graph = self.initial_graph
        if len(graph.nodes()) == 0: raise AssertionError('nothing to filter!')
        if start not in graph.nodes(): raise AssertionError('The start node is not in the graph')
        if goal not in graph.nodes(): raise AssertionError('The goal node is not in the graph')
        for component in _nx.connected_component_subgraphs(graph):
            if start in component.nodes(): break
        try:
            if only_shortest: paths = _nx.all_shortest_paths(component, start, goal)
            else: paths = _nx.all_simple_paths(component, start, goal, cutoff=max_len)
        except: raise AssertionError('start and goal nodes are not connected')
        shortest_graph = _nx.Graph()
        for path in paths:
            shortest_graph.add_edges_from(_hf.pairwise(path))
        self.applied_filters['shortest_paths']=(start, goal)
        self.filtered_graph = shortest_graph
        self._generate_filtered_results_from_filtered_graph()

    def filter_single_path(self, *nodes, use_filtered=True):
        if use_filtered: graph = self.filtered_graph
        else: graph = self.initial_graph
        if len(graph.nodes()) == 0: raise AssertionError('nothing to filter!')
        keep_edges = []
        for resa, resb in _hf.pairwise(nodes):
            edge, edge_check = (resa, resb), (resb, resa)
            if edge in graph.edges():
                keep_edges.append(edge)
            elif edge_check in graph.edges():
                keep_edges.append(edge_check)
            else:
                raise AssertionError('There is no connection between {} and {}'.format(resa, resb))
        self.applied_filters['single_path']=nodes
        self.filtered_graph = _nx.Graph()
        self.filtered_graph.add_edges_from(keep_edges)
        self._generate_filtered_results_from_filtered_graph()

    def filter_between_segnames(self, segna, segnb=None, use_filtered=True):
        if use_filtered: results = self.filtered_results
        else: results = self.initial_results
        if len(results) == 0: raise AssertionError('nothing to filter!')
        segnames = [segna, segnb]
        keep_bonds = []
        for key in results:
            sa,_,_, sb,_,_ = _hf.deconst_key(key, self.residuewise)
            if ((sa in segnames) and (sb in segnames) and (sa!=sb) and (segnb!=None)) or (_np.logical_xor((sa in segnames),(sb in segnames)) and (segnb==None)):
                keep_bonds.append(key)
        self.applied_filters['segnames']=(segna,segnb)
        self.filtered_results = {key:results[key] for key in keep_bonds}
        self._generate_filtered_graph_from_filtered_results()

    def filter_between_resnames(self, resna, resnb=None, use_filtered=True):
        if use_filtered: results = self.filtered_results
        else: results = self.initial_results
        if len(results) == 0: raise AssertionError('nothing to filter!')
        resnames = [resna, resnb]
        keep_bonds = []
        for key in results:
            _,ra,_, _,rb,_ = _hf.deconst_key(key, self.residuewise)
            if ((ra in resnames) and (rb in resnames) and (ra!=rb) and (resnb!=None)) or (_np.logical_xor((ra in resnames),(rb in resnames)) and (resnb==None)):
                keep_bonds.append(key)
        self.applied_filters['resnames']=(resna,resnb)
        self.filtered_results = {key:results[key] for key in keep_bonds}
        self._generate_filtered_graph_from_filtered_results()

    def filter_to_frame(self, frame, use_filtered=True):
        _ = self._compute_graph_in_frame(frame, True, use_filtered)

    def _lipid_analysis(self, max_length, correction_hba, between_nodes=None, filter_highest_occupancy=False, use_filtered=True):
        if use_filtered:
            graph = self.filtered_graph
            if filter_highest_occupancy: results = self.filtered_results
        else:
            graph = self.initial_graph
            if filter_highest_occupancy: results = self.initial_results
        correction_dict = {}
        for key in correction_hba.filtered_results:
            keya, keyb = key.split(':')
            _, resna, resida, _, resnb, residb = _hf.deconst_key(key, True)
            if resna == 'POPC' :
                correction_dict[keya] = _np.zeros(self.nb_frames, dtype=_np.bool)
            elif resnb == 'POPC' :
                correction_dict[keyb] = _np.zeros(self.nb_frames, dtype=_np.bool)
        for key in correction_hba.filtered_results:
            keya, keyb = key.split(':')
            _, resna, resida, _, resnb, residb = _hf.deconst_key(key, True)
            if resna == 'POPC' :
                correction_dict[keya] |= correction_hba.filtered_results[key]
            elif resnb == 'POPC' :
                correction_dict[keyb] |= correction_hba.filtered_results[key]
        if between_nodes != None: nodes = between_nodes
        else: nodes = [key for key in correction_dict]
        res_length = {}
        res_path = {}
        res_occ = {}
        plot_table = _np.zeros((max_length,self.nb_frames), dtype=int)
        for keya, keyb in combinations(nodes, 2):
            key = keya+':'+keyb
            temp_len=[]
            temp_paths = []
            combined_occupancy = _np.zeros(self.nb_frames, dtype=_np.bool)
            try:
                _nx.shortest_path(graph, keya, keyb)
                for frame in range(self.nb_frames):
                    if not correction_dict[keya][frame]: continue
                    if not correction_dict[keyb][frame]: continue
                    frame_graph = self._compute_graph_in_frame(frame, use_filtered=use_filtered)
                    try:
                        path = _nx.shortest_path(frame_graph, keya, keyb)
                        combined_occupancy[frame]=True
                        temp_len.append(len(path))
                        if path not in temp_paths: temp_paths.append(path)
                        try:
                            plot_table[len(path)-2][frame] += 1
                        except:
                            pass
                    except: continue
                if not temp_len: continue
                if filter_highest_occupancy:
                    combined_occupancy = []
                    for path in temp_paths:
                        combined_occupancy.append(_np.ones(self.nb_frames, dtype=_np.bool))
                        for nodea, nodeb in _hf.pairwise(path):
                            try: combined_occupancy[-1] &= results[':'.join((nodea, nodeb))]
                            except KeyError: combined_occupancy[-1] &= results[':'.join((nodeb, nodea))]
                    for i,c in enumerate(combined_occupancy):
                        combined_occupancy[i]=c.mean()
                    max_index = _np.argmax(combined_occupancy)
                    res_occ[key] = combined_occupancy[max_index]
                    res_path[key] = [int(node.split('-')[2]) for node in temp_paths[max_index]]
                    res_length[key] = len(res_path[key])
                else:
                    occ = combined_occupancy.mean()
                    res_occ[key]=occ
                    res_path[key]= [[int(node.split('-')[2]) for node in path] for path in temp_paths]
                    res_length[key] = _np.array(temp_len).mean()
            except: pass
        return res_length, res_path, res_occ, plot_table

    def compute_all_pairs_shortest_path_per_frame(self, max_length=5, between_nodes=None, filter_highest_occupancy=False, use_filtered=True):
        if use_filtered:
            graph = self.filtered_graph
            if filter_highest_occupancy: results = self.filtered_results
        else:
            graph = self.initial_graph
            if filter_highest_occupancy: results = self.initial_results
        if between_nodes != None: nodes = between_nodes
        else: nodes = graph.nodes()
        res_length = {}
        res_path = {}
        res_occ = {}
        plot_table = _np.zeros((max_length,self.nb_frames), dtype=int)
        for keya, keyb in combinations(nodes, 2):
            key = keya+':'+keyb
            temp_len=[]
            temp_paths = []
            combined_occupancy = _np.zeros(self.nb_frames, dtype=_np.bool)
            try:
                _nx.shortest_path(graph, keya, keyb)
                for frame in range(self.nb_frames):
                    frame_graph = self._compute_graph_in_frame(frame, use_filtered=use_filtered)
                    try:
                        path = _nx.shortest_path(frame_graph, keya, keyb)
                        combined_occupancy[frame]=True
                        temp_len.append(len(path))
                        if path not in temp_paths: temp_paths.append(path)
                        try:
                            plot_table[len(path)-2][frame] += 1
                        except:
                            pass
                    except: continue
                if not temp_len: continue
                if filter_highest_occupancy:
                    combined_occupancy = []
                    for path in temp_paths:
                        combined_occupancy.append(_np.ones(self.nb_frames, dtype=_np.bool))
                        for nodea, nodeb in _hf.pairwise(path):
                            try: combined_occupancy[-1] &= results[':'.join((nodea, nodeb))]
                            except KeyError: combined_occupancy[-1] &= results[':'.join((nodeb, nodea))]
                    for i,c in enumerate(combined_occupancy):
                        combined_occupancy[i]=c.mean()
                    max_index = _np.argmax(combined_occupancy)
                    res_occ[key] = combined_occupancy[max_index]
                    res_path[key] = [int(node.split('-')[2]) for node in temp_paths[max_index]]
                    res_length[key] = len(res_path[key])
                else:
                    occ = combined_occupancy.mean()
                    res_occ[key]=occ
                    res_path[key]= [[int(node.split('-')[2]) for node in path] for path in temp_paths]
                    res_length[key] = _np.array(temp_len).mean()
            except: pass
        return res_length, res_path, res_occ, plot_table

    def compute_all_pairs_shortest_path(self, max_length=5, between_nodes=None, use_filtered=True):
        if use_filtered:
            graph = self.filtered_graph
            results = self.filtered_results
        else:
            graph = self.initial_graph
            results = self.initial_results
        if between_nodes != None: nodes = between_nodes
        else: nodes = graph.nodes()
        res_length = {}
        res_path = {}
        res_occ = {}
        plot_table = _np.zeros((max_length,self.nb_frames))
        for keya, keyb in combinations(nodes, 2):
            try:
                key = keya+':'+keyb
                path = _nx.shortest_path(graph, keya, keyb)
                combined_occupancy = _np.ones(self.nb_frames, dtype=_np.bool)
                for nodea, nodeb in _hf.pairwise(path):
                    try: combined_occupancy &= results[':'.join((nodea, nodeb))]
                    except KeyError: combined_occupancy &= results[':'.join((nodeb, nodea))]
                res_occ[key]=combined_occupancy.mean()
                res_path[key]=[int(node.split('-')[2]) for node in path]
                res_length[key] = len(res_path[key])
                try:
                    plot_table[res_length[key]-2] += combined_occupancy
                except:
                    pass
            except _nx.NetworkXNoPath:
                pass
        return res_length, res_path, res_occ, plot_table

    def compute_total_occupancy(self, use_filtered=True):
        if use_filtered: results = self.filtered_results
        else: results = self.initial_results
        return_res = {}
        for key in results:
            resa, resb = key.split(':')
            try: return_res[resa] |= results[key]
            except KeyError: return_res[resa] = results[key]
            try: return_res[resb] |= results[key]
            except KeyError: return_res[resb] = results[key]
        return return_res

    def compute_joint_occupancy(self, use_filtered=True):
        if use_filtered: results = self.filtered_results
        else: results = self.initial_results
        combined_occupancy = _np.ones(self.nb_frames, dtype=_np.bool)
        for res in results: combined_occupancy &= results[res]
        self.joint_occupancy_series = combined_occupancy
        self.joint_occupancy_frames = _np.nonzero(combined_occupancy)[0]
        return combined_occupancy.mean()

    def compute_centrality(self, centrality_type='betweenness', normalize=True, weight=None, filename=None, use_filtered=True):
        if use_filtered: graph = self.filtered_graph
        else: graph = self.initial_graph
        centralities = {node:_np.zeros(self.nb_frames) for node in graph.nodes()}
        for i in range(self.nb_frames):
            g_i = self._compute_graph_in_frame(i, use_filtered=use_filtered)
            if centrality_type == 'betweenness': centrality_i = _nx.betweenness_centrality(g_i, normalized=normalize)
            elif centrality_type == 'degree': centrality_i = _nx.degree_centrality(g_i)
            else: raise AssertionError("centrality_type has to be 'betweenness' or 'degree'")
            for node in centrality_i:
                if (centrality_type == 'degree') and not normalize: centralities[node][i] = centrality_i[node] * (g_i.number_of_nodes()-1)
                else: centralities[node][i] = centrality_i[node]
        if weight is None:
            for node in centralities: centralities[node] = _np.round(centralities[node].mean(),2)
        else:
            ma = max(centralities.values())
            for node in centralities: centralities[node] = _np.round((centralities[node].mean()/ma)/weight,2)
        if filename is not None:
            string = centrality_type + ' centrality values:\n\n'
            mc = max(centralities.values())
            for c in centralities:
                string += c + ' ' + str(centralities[c]) + '\n'
            string += '\nnormalized to [0,1]:\n'
            for c in centralities:
                string += c + ' ' + str(_np.round(centralities[c]/mc, 2)) + '\n'
            with open(filename, 'w') as f:
                f.write(string)
        return centralities

    def compute_average_position_of_connection(self, use_filtered=True):
        if use_filtered: results = self.filtered_results
        else: results = self.initial_results
        position = {key:_np.zeros(3) for key in results}
        if self.residuewise: all_id = _np.array(self._all_ids+self._ions_ids)
        else: all_id = _np.array(self._all_ids_atomwise+self._ions_ids_atomwise)
        all_coordinates = _np.vstack((self._da_selection.positions, self._water.positions, _np.array(self._ions.positions).reshape((-1,3))))
        residues = _np.array([(key.split(':')[0], key.split(':')[1]) for key in results])
        #bonds_binary = np.array([results[key] for key in results]).reshape((len(results), -1))
        for i,ts in enumerate(self._universe.trajectory[self._trajectory_slice]):
            for ii, bond in enumerate(results):
                if results[bond][i]: position[bond] += (all_coordinates[all_id == residues[ii][0]].reshape((-1, 3)).mean(0) + all_coordinates[all_id == residues[ii][1]].reshape((-1, 3)).mean(0)) / 2
        for key in position:
            position[key] /= results[key].sum()
        self._connection_position = position
        return position

    def draw_multi_segment_connection_timeseries(self, segnames=None, colors=None, use_filtered=True, filename=None, return_figure=False):
        if use_filtered: results = self.filtered_results
        else: results = self.initial_results
        if segnames == None:
            segnames = _np.unique([[_hf.deconst_key(key, self.residuewise)[0], _hf.deconst_key(key, self.residuewise)[3]] for key in results])
        if colors != None:
            segnames = [s for s in colors]
        segname_results = {segname:{} for segname in segnames}
        for segnamea, segnameb in combinations(segnames, 2): segname_results[segnamea+'-'+segnameb]={}
        for key in results:
            res1, res2 = key.split(':')
            sn1, sn2 = res1.split('-')[0], res2.split('-')[0]
            comb, comb_check = '-'.join((sn1,sn2)), '-'.join((sn2,sn1))
            try: segname_results[comb][key] = results[key]
            except:
                try: segname_results[comb_check][key] = results[key]
                except:
                    try: segname_results[sn1][key] = results[key]
                    except:
                        try: segname_results[sn2][key] = results[key]
                        except: pass

        if colors == None: colors = {}
        cmap = matplotlib.cm.get_cmap('Spectral')
        for i,s in enumerate(segname_results):
            if s in colors: continue
            colors[s] = cmap(_np.linspace(0,1,len(segname_results))[i])
        fig, ax = _plt.subplots()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        _plt.xlabel('Frame' , fontsize = 16)
        _plt.ylabel('Count' , fontsize = 16)
        return_string = ''
        for i,seg in enumerate(segname_results):
            return_string += seg + ' '
            res_temp = _np.zeros(self.nb_frames, dtype=_np.int)
            for bond in segname_results[seg]: res_temp += segname_results[seg][bond]
            return_string += ' '.join(_np.array(res_temp, dtype=_np.str)) + '\n'
            _plt.plot(_np.arange(self.nb_frames), res_temp, color=colors[seg], label=seg)
        _plt.legend()
        if return_figure:
            _plt.close()
            return fig, _hf.string_in_columns(return_string)
        self._save_or_draw(filename, return_figure=return_figure)


    def draw_multi_segname_occupancy_histogram(self, min_occupancy=0.0, max_occupancy=1.0, occupancy_step=0.1, segnames=None, colors=None, use_filtered=True, filename=None, return_figure=False):
        if use_filtered: results = self.filtered_results
        else: results = self.initial_results
        if segnames == None:
            segnames = _np.unique([[_hf.deconst_key(key, self.residuewise)[0], _hf.deconst_key(key, self.residuewise)[3]] for key in results])
        if colors != None:
            segnames = [s for s in colors]
        occupancies = _np.arange(min_occupancy, max_occupancy, occupancy_step)
        segname_results = {segname:{} for segname in segnames}
        temp = []
        for segnamea, segnameb in combinations(segname_results, 2): temp.append(segnamea+'-'+segnameb)
        for comb in temp: segname_results[comb]={}
        for key in results:
            res1, res2 = key.split(':')
            sn1, sn2 = res1.split('-')[0], res2.split('-')[0]
            comb, comb_check = '-'.join((sn1,sn2)), '-'.join((sn2,sn1))
            try: segname_results[comb][key] = results[key]
            except:
                try: segname_results[comb_check][key] = results[key]
                except:
                    try: segname_results[sn1][key] = results[key]
                    except:
                        try: segname_results[sn2][key] = results[key]
                        except: pass

        res = _np.zeros((len(segname_results), len(occupancies)), dtype=_np.int)
        if colors == None: colors = {}
        cmap = matplotlib.cm.get_cmap('Spectral')
        for i,s in enumerate(segname_results):
            if s in colors: continue
            colors[s] = cmap(_np.linspace(0,1,len(segname_results))[i])
        labels = ['{0:.{1}f}'.format(occupancy, 2) for occupancy in occupancies]
        fig, ax = _plt.subplots()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        _plt.xticks(range(len(occupancies)), labels)
        _plt.xlabel('Occupancy' , fontsize = 16)
        _plt.ylabel('Count' , fontsize = 16)
        handles = []
        return_string = 'occupancies ' + ' '.join(_np.array(occupancies.round(2), dtype=_np.str)) + '\n'
        for i,segname in enumerate(segnames):
            try:
                res[i] = _np.array([len(_hf.filter_occupancy(segname_results[segname], occupancy)) for occupancy in occupancies])
            except:
                res[i] = 0
            return_string += segname + ' ' + ' '.join(_np.array(res[i], dtype=_np.str)) + '\n'
            handles.append(_plt.bar(range(len(res[i])), res[i], width=0.8, label=segname, color=colors[segname], bottom=res[:i].sum(0)))
        for segname in segname_results:
            if segname in segnames: continue
            i+=1
            try:
                res[i] = _np.array([len(_hf.filter_occupancy(segname_results[segname], occupancy)) for occupancy in occupancies])
            except:
                res[i] = 0
            return_string += segname + ' ' + ' '.join(_np.array(res[i], dtype=_np.str))+ '\n'
            handles.append(_plt.bar(range(len(res[i])), res[i], width=0.8, label=segname, color=colors[segname], bottom=res[:i].sum(0)))
        _plt.legend(handles=handles)
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        if return_figure:
            _plt.close()
            return fig, _hf.string_in_columns(return_string)
        self._save_or_draw(filename, return_figure=return_figure)

    def draw_sum_of_connections_timeseries(self, compare_to=None,legend_text=None, use_filtered=True, filename=None, return_figure=False):
        if use_filtered: results = self.filtered_results
        else: results = self.initial_results
        res = _np.zeros(self.nb_frames, dtype=_np.int)
        for key in results:
            res += results[key]
        if compare_to != None:
            if use_filtered: results_compare = compare_to.filtered_results
            else: results_compare = compare_to.initial_results
            res_c = _np.zeros(compare_to.nb_frames, dtype=_np.int)
            for key in results_compare:
                res_c += results_compare[key]
        fig, ax = _plt.subplots()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        _plt.xlabel('Frame' , fontsize = 16)
        _plt.ylabel('Count' , fontsize = 16)
        return_string = ''
        if legend_text != None:
            _plt.plot(_np.arange(self.nb_frames), res, label=legend_text[0])
            return_string += legend_text[0].replace(' ', '_') + ' '.join(_np.array(res, dtype=_np.str)) + '\n'
            if compare_to != None:
                _plt.plot(_np.arange(compare_to.nb_frames), res_c, label=legend_text[1])
                return_string += legend_text[0].replace(' ', '_') + ' '.join(_np.array(res_c, dtype=_np.str)) + '\n'
            _plt.legend()
        else:
            _plt.plot(_np.arange(self.nb_frames), res)
            return_string += 'connections_per_frame ' + ' '.join(_np.array(res, dtype=_np.str)) + '\n'
            if compare_to != None:
                _plt.plot(_np.arange(compare_to.nb_frames), res_c)
                return_string += 'compared_connections_per_frame ' + ' '.join(_np.array(res_c, dtype=_np.str)) + '\n'
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        if return_figure:
            _plt.close()
            return fig, _hf.string_in_columns(return_string)
        self._save_or_draw(filename, return_figure=return_figure)

    def draw_connections_per_residue(self, average=False, residues_to_plot=None, compare_to=None, xtick_resnames=True, legend_text=None, segnames=None, use_filtered=True, filename=None, return_figure=False):
        if use_filtered: results = self.filtered_results
        else: results = self.initial_results
        if segnames == None: segnames = _np.unique([[_hf.deconst_key(key, self.residuewise)[0], _hf.deconst_key(key, self.residuewise)[3]] for key in results])
        segnames = _np.array(segnames)
        resnames = {}
        for key in results:
            _, resna, resia, _, resnb, resib = _hf.deconst_key(key, self.residuewise)
            resnames[resia]= resna
            resnames[resib]= resnb
        residues = _np.unique([_hf.deconst_key(key, self.residuewise)[2] for key in results]+[_hf.deconst_key(key, self.residuewise)[5] for key in results])
        if average: results_per_residue = {segname:_odict({i:_np.zeros(self.nb_frames) for i in residues}) for segname in segnames}
        else: results_per_residue = {segname:_odict({i:0 for i in residues}) for segname in segnames}
        if compare_to != None:
            if use_filtered: results_compare = compare_to.filtered_results
            else: results_compare = compare_to.initial_results
            residues_compare = _np.unique([_hf.deconst_key(key, compare_to.residuewise)[2] for key in results_compare]+[_hf.deconst_key(key, compare_to.residuewise)[5] for key in results_compare])
            for key in results_compare:
                _, resna, resia, _, resnb, resib = _hf.deconst_key(key, compare_to.residuewise)
                resnames[resia]= resna
                resnames[resib]= resnb
            segnames_compare = _np.unique([[_hf.deconst_key(key, compare_to.residuewise)[0], _hf.deconst_key(key, compare_to.residuewise)[3]] for key in results_compare])
            if average: results_per_residue_compare = {segname:_odict({i:_np.zeros(compare_to.nb_frames) for i in residues_compare}) for segname in segnames_compare}
            else: results_per_residue_compare = {segname:_odict({i:0 for i in residues_compare}) for segname in segnames_compare}
        for key in results:
            segna, resna, resida, segnb, resnb, residb = _hf.deconst_key(key, self.residuewise)
            if resna not in ['TIP3', 'HOH']+self._ions_list:
                if average: results_per_residue[segna][resida] += results[key]
                else: results_per_residue[segna][resida] +=1
            if resnb not in ['TIP3', 'HOH']+self._ions_list:
                if average: results_per_residue[segnb][residb] += results[key]
                else: results_per_residue[segnb][residb] +=1
        if average:
            for segn in results_per_residue:
                for key in results_per_residue[segn]:
                    results_per_residue[segn][key] = results_per_residue[segn][key].mean()
        if compare_to != None:
            for key in results_compare:
                segna, resna, resida, segnb, resnb, residb = _hf.deconst_key(key, compare_to.residuewise)
                if resna not in ['TIP3', 'HOH']+self._ions_list:
                    if average: results_per_residue_compare[segna][resida] += results_compare[key]
                    else: results_per_residue_compare[segna][resida] +=1
                if resnb not in ['TIP3', 'HOH']+self._ions_list:
                    if average: results_per_residue_compare[segnb][residb] += results_compare[key]
                    else: results_per_residue_compare[segnb][residb] +=1
            if average:
                for segn in results_per_residue_compare:
                    for key in results_per_residue_compare[segn]:
                        results_per_residue_compare[segn][key] = results_per_residue_compare[segn][key].mean()

        if compare_to != None:
            all_residues = _np.unique(list(residues)+list(residues_compare))
            all_segnames = _np.unique(list(segnames)+list(segnames_compare))
            for segn in all_segnames:
                if segn not in results_per_residue_compare: results_per_residue_compare[segn] = {i:0 for i in all_residues}
                if segn not in results_per_residue: results_per_residue[segn] = {i:0 for i in all_residues}
                for key in all_residues:
                    if key not in results_per_residue[segn]: results_per_residue[segn][key]=0
                    if key not in results_per_residue_compare[segn]: results_per_residue_compare[segn][key]=0

            for segn in results_per_residue: results_per_residue[segn] = _odict({key:results_per_residue[segn][key] for key in sorted(list(results_per_residue[segn].keys()))})
            for segn in results_per_residue_compare: results_per_residue_compare[segn] = _odict({key:results_per_residue_compare[segn][key] for key in sorted(list(results_per_residue_compare[segn].keys()))})
            ys_compare = {segn:[] for segn in results_per_residue}

        xs, ys = {segn:[] for segn in results_per_residue}, {segn:[] for segn in results_per_residue}
        for segn in results_per_residue:
            for key in results_per_residue[segn]:
                xs[segn].append(key)
                ys[segn].append(results_per_residue[segn][key])
                if compare_to != None:
                    ys_compare[segn].append(results_per_residue_compare[segn][key])

        for segn in xs: xs[segn] = _np.array(xs[segn])
        for segn in ys: ys[segn] = _np.array(ys[segn])
        if compare_to != None:
            for segn in ys_compare: ys_compare[segn] = _np.array(ys_compare[segn])
        if residues_to_plot != None:
            r_index = {segn:[] for segn in results_per_residue}
            for segn in xs:
                for i,x in enumerate(xs[segn]):
                    if x in residues_to_plot: r_index[segn].append(i)
                xs[segn] = xs[segn][r_index[segn]]
                ys[segn] = ys[segn][r_index[segn]]
                if compare_to != None:
                    ys_compare[segn] = ys_compare[segn][r_index[segn]]

        ys_bottom, ys_compare_bottom, ss = {segn:[] for segn in results_per_residue}, {segn:[] for segn in results_per_residue}, []
        for segname in xs:
            ss.append(segname)
        ys_bottom[ss[0]] = _np.zeros(len(ys[ss[0]]))
        if compare_to != None: ys_compare_bottom[ss[0]] = _np.zeros(len(ys_compare[ss[0]]))
        for i in range(1,len(ss)):
            ys_bottom[ss[i]]=ys[ss[i-1]]
            if compare_to!=None:ys_compare_bottom[ss[i]] = ys_compare[ss[i-1]]
        fig, ax = _plt.subplots()
        for segn in xs:

            if compare_to != None:
                if legend_text != None:
                    if ys[segn].sum()>0:_plt.bar(_np.arange(len(ys[segn])) - 0.1, ys[segn], width=0.2, label=legend_text[0]+' - '+segn, bottom=ys_bottom[segn])
                    if ys_compare[segn].sum()>0:_plt.bar(_np.arange(len(ys_compare[segn])) + 0.1, ys_compare[segn], width=0.2, label=legend_text[1]+' - '+segn, bottom=ys_compare_bottom[segn])
                    _plt.legend()
                else:
                    _plt.bar(_np.arange(len(ys[segn])) - 0.1, ys[segn], width=0.2, label=segn, bottom=ys_bottom[segn])
                    _plt.bar(_np.arange(len(ys_compare[segn])) + 0.1, ys_compare[segn], width=0.2, label=segn, bottom=ys_compare_bottom[segn])
                    _plt.legend()
            else:
                if legend_text != None:
                    if ys[segn].sum()>0:_plt.bar(_np.arange(len(ys[segn])), ys[segn], width=0.4, label = legend_text[0]+' - '+segn, bottom=ys_bottom[segn])
                    _plt.legend()
                else:
                    _plt.bar(_np.arange(len(ys[segn])), ys[segn], width=0.4, label=segn, bottom=ys_bottom[segn])
                    _plt.legend()

        all_xs = _np.unique([xs[key] for key in xs])
        all_labels = [_hf.aa_three2one[resnames[resi]]+str(resi+self.add_missing_residues) if resnames[resi] in _hf.aa_three2one else resnames[resi]+str(resi+self.add_missing_residues) for resi in all_xs]
        return_string = 'labels ' + ' '.join(all_labels) + '\n'
        for segn in xs: return_string += segn + ' ' + ' '.join(_np.array(ys[segn], dtype=_np.str)) + '\n'
        if compare_to != None:
            for segn in xs: return_string += 'compare_' + segn + ' ' + ' '.join(_np.array(ys_compare[segn], dtype=_np.str)) + '\n'
        _plt.xticks(_np.arange(-1, all_xs.size, 1.0), ['']+all_labels, rotation=45)
        _plt.xlim([-0.5, all_xs.size-0.5])
        _plt.xlabel('Residue' , fontsize = 16)
        _plt.ylabel('Count' , fontsize = 16)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        #ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        if return_figure:
            _plt.close()
            return fig, _hf.string_in_columns(return_string)
        self._save_or_draw(filename, return_figure=return_figure)

    def draw_joint_timeseries(self, use_filtered=True, scatter_size=0.5, filename=None, return_figure=False):
        if use_filtered: results = self.filtered_results
        else: results = self.initial_results
        ts = _np.ones(self.nb_frames, dtype=bool)
        for key in results:
            ts &= results[key]
        fig, ax = _plt.subplots()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        _plt.scatter(range(len(ts)), ts, s=scatter_size)
        _plt.title('Joint Occupancy: '+str(_np.round(ts.mean()*100,1)))
        _plt.xlabel('Frame' , fontsize = 16)
        ax.set_xlim(0, self.nb_frames-1)
        _plt.yticks([-0.5,0.01,1.01,1.5], ['', 'false', 'true', ''])
        _plt.tick_params(
                axis='y',
                which='both',
                bottom=False,
                top=False,
                labelbottom=False,
                length=0.0,
                width=0.0)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        return_string = 'Joint_occupancy ' + ' '.join(_np.array(ts.astype(_np.int), dtype=_np.str))
        if return_figure:
            _plt.close()
            return fig, _hf.string_in_columns(return_string)
        self._save_or_draw(filename, return_figure=return_figure)

    def draw_compare_joint_timeseries(self, other_paths, descriptor='path', scatter_size=0.5, use_filtered=True, filename=None, return_figure=False):
        if use_filtered: results = self.filtered_results
        else: results = self.initial_results
        results = [results]
        for other in other_paths:
            if use_filtered: results.append(other.filtered_results)
            else: results.append(other.initial_results)
        ts = _np.ones((len(results), self.nb_frames), dtype=bool)
        for i, result in enumerate(results):
            for key in result:
                ts[i] &= result[key]
        fig, ax = _plt.subplots()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        #ax.spines['left'].set_visible(False)
        occupancies = ts.mean(1)
        sorted_index = _np.argsort(occupancies)
        for i,ind in enumerate(sorted_index):
            _plt.scatter(range(len(ts[ind])), _np.array(ts[ind], dtype=int) * (i+1), s=scatter_size)
        _plt.xlabel('Frame' , fontsize = 16)
        ax.set_xlim(0, self.nb_frames-1)
        ax.set_ylim(0.5, len(results)+0.5)
        #_plt.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False, length=0.0, width=0.0)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax2 = ax.twinx()
        ax2.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax2.set_ylim(0.5, len(results)+0.5)
        ax.set_yticks(range(1, len(ts)+1))
        ax.set_yticklabels(['{}'.format(o*100) for o in occupancies[sorted_index].round(1)])
        ax.set_ylabel('Joint Occupancy [%]', fontsize=16)
        ax2.set_yticks(range(1, len(ts)+1))
        ax2.set_yticklabels([(descriptor + ' {}').format(i+1) for i in range(len(ts))])
        ax2.spines['right'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        return_string = ''
        for i,ind in enumerate(sorted_index):
            return_string += 'Joint_occupancy_' + str(i+1) + ' ' + ' '.join(_np.array(ts[ind].astype(_np.int), dtype=_np.str)) + '\n'
        if return_figure:
            _plt.close()
            fig.add_axes(ax2,  label='paths')
            return fig, _hf.string_in_columns(return_string)

        self._save_or_draw(filename, return_figure=return_figure)

    def draw_occupancy_histogram(self, min_occupancy=0.0, max_occupancy=1.0, occupancy_step=0.1, compare_to=None, legend_text=None, use_filtered=True, filename=None, return_figure=False):
        if use_filtered: results = self.filtered_results
        else: results = self.initial_results
        if compare_to != None:
            if use_filtered: results_compare = compare_to.filtered_results
            else: results_compare = compare_to.initial_results
        occupancies = _np.arange(min_occupancy, max_occupancy, occupancy_step)
        bond_counts = [len(_hf.filter_occupancy(results, occupancy)) for occupancy in occupancies]
        diff_counts = _np.array(bond_counts)
        labels = ['{0:.{1}f}'.format(occupancy, 2) for occupancy in occupancies]
        fig, ax = _plt.subplots()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        return_string = 'occupancies ' + ' '.join(_np.array(occupancies.round(2), dtype=_np.str))  + '\n'
        if compare_to != None:
            diff_counts_c = _np.array([len(_hf.filter_occupancy(results_compare, occupancy)) for occupancy in occupancies])
            if legend_text!=None:
                _plt.bar(_np.array(list(range(len(diff_counts))))-0.2, diff_counts, width=0.4, label=legend_text[0])
                _plt.bar(_np.array(list(range(len(diff_counts))))+0.2, diff_counts_c, width=0.4, label=legend_text[1])
                return_string += legend_text[0] + ' ' + ' '.join(_np.array(diff_counts, dtype=_np.str)) + '\n'
                return_string += legend_text[1] + ' ' + ' '.join(_np.array(diff_counts_c, dtype=_np.str)) + '\n'
                _plt.legend()
            else:
                _plt.bar(_np.array(list(range(len(diff_counts))))-0.2, diff_counts, width=0.4)
                _plt.bar(_np.array(list(range(len(diff_counts))))+0.2, diff_counts_c, width=0.4)
                return_string += 'main' + ' ' + ' '.join(_np.array(diff_counts, dtype=_np.str)) + '\n'
                return_string += 'compare' + ' ' + ' '.join(_np.array(diff_counts_c, dtype=_np.str)) + '\n'
        else:
            if legend_text != None:
                _plt.bar(range(len(diff_counts)), diff_counts, label=legend_text[0])
                return_string +=  legend_text[0] + ' ' + ' '.join(_np.array(diff_counts, dtype=_np.str)) + '\n'
                _plt.legend()
            else:
                _plt.bar(range(len(diff_counts)), diff_counts)
                return_string += 'main' + ' ' + ' '.join(_np.array(diff_counts, dtype=_np.str)) + '\n'
        _plt.xticks(range(len(diff_counts)), labels)
        _plt.xlabel('Occupancy' , fontsize = 16)
        _plt.ylabel('Count' , fontsize = 16)
        if return_figure:
            _plt.close()
            return fig, _hf.string_in_columns(return_string)
        self._save_or_draw(filename, return_figure=return_figure)

    def draw_compare_occupancy_histogram(self, compare_to, xaxis_names=None, use_filtered=True, filename=None, return_figure=False):
        if use_filtered: results = self.filtered_results
        else: results = self.initial_results
        results = [results]
        lengths = [self.nb_frames]
        for other in compare_to:
            lengths.append(other.nb_frames)
            if use_filtered: results.append(other.filtered_results)
            else: results.append(other.initial_results)
        if xaxis_names != None: labels = xaxis_names
        else: labels = _np.arange(1,len(results)+1)
        fig, ax = _plt.subplots()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        occupancies = _np.zeros(len(results))
        for i,result in enumerate(results):
            temp = _np.ones(lengths[i], dtype=bool)
            for key in result:
                temp &= result[key]
            occupancies[i] = temp.mean()
        _plt.bar(_np.arange(len(occupancies)), occupancies)
        _plt.xticks(range(occupancies.size), labels)
        _plt.ylabel('Occupancy' , fontsize = 16)
        if return_figure:
            _plt.close()
            return fig, 'NoData'
        self._save_or_draw(filename, return_figure=return_figure)

    def draw_time_histogram(self, nb_blocks=10, first_frame=0, last_frame=-1, compare_to=None, legend_text=None, use_filtered=True, filename=None, return_figure=False):
        if use_filtered: results = self.filtered_results
        else: results = self.initial_results
        if compare_to != None:
            if use_filtered: results_compare = compare_to.filtered_results
            else: results_compare = compare_to.initial_results
        values = _np.array([results[key] for key in results])
        if (last_frame < 0) or (last_frame > self.nb_frames): last = self.nb_frames
        else: last = last_frame
        block_values = _np.linspace(first_frame, last, nb_blocks, dtype=int)
        y_hist = []
        for i, ii in _hf.pairwise(block_values): y_hist.append(values[:,i:ii].sum())
        if compare_to != None:
            values = _np.array([results_compare[key] for key in results_compare])
            block_values = _np.linspace(0, self.nb_frames, nb_blocks, dtype=int)
            y_hist_c = []
            for i, ii in _hf.pairwise(block_values): y_hist_c.append(values[:,i:ii].sum())
        labels = ['{}'.format(block_value) for block_value in block_values]
        fig, ax = _plt.subplots()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        return_string = 'frame_from ' + ' '.join([labels[i] for i in range(len(labels)-1)]) + '\n'
        return_string += 'frame_to ' + ' '.join([labels[i+1] for i in range(len(labels)-1)]) + '\n'
        if compare_to != None:
            if legend_text!=None:
                _plt.bar(_np.array(list(range(len(y_hist))))+0.3, y_hist, width=0.4, label=legend_text[0])
                _plt.bar(_np.array(list(range(len(y_hist_c))))+0.7, y_hist_c, width=0.4, label=legend_text[1])
                return_string += legend_text[0] + ' ' + ' '.join(_np.array(y_hist, dtype=_np.str)) + '\n'
                return_string += legend_text[1] + ' ' + ' '.join(_np.array(y_hist_c, dtype=_np.str)) + '\n'
                _plt.legend()
            else:
                _plt.bar(_np.array(list(range(len(y_hist))))+0.3, y_hist, width=0.4)
                _plt.bar(_np.array(list(range(len(y_hist_c))))+0.7, y_hist_c, width=0.4)
                return_string += 'main' + ' ' + ' '.join(_np.array(y_hist, dtype=_np.str)) + '\n'
                return_string += 'compare' + ' ' + ' '.join(_np.array(y_hist_c, dtype=_np.str)) + '\n'
        else:
            if legend_text != None:
                _plt.bar(_np.array(list(range(len(y_hist))))+0.5, y_hist, label=legend_text[0])
                return_string += legend_text[0] + ' ' + ' '.join(_np.array(y_hist, dtype=_np.str)) + '\n'
                _plt.legend()
            else:
                _plt.bar(_np.array(list(range(len(y_hist))))+0.5, y_hist)
                return_string += 'main' + ' ' + ' '.join(_np.array(y_hist, dtype=_np.str)) + '\n'
        #_plt.bar(range(len(y_hist)), y_hist)
        _plt.xticks(range(len(y_hist)), labels)
        _plt.xlabel('Frame' , fontsize = 16)
        _plt.ylabel('Count' , fontsize = 16)
        if return_figure:
            _plt.close()
            return fig, _hf.string_in_columns(return_string)
        self._save_or_draw(filename, return_figure=return_figure)

    def draw_residue_range_heatmap(self, ranges, names, average=False, label='Region', use_filtered=True, filename=None, return_figure=False):
        if use_filtered: results = self.filtered_results
        else: results = self.initial_results
        residues = []
        for r in ranges: residues+=list(range(r[0],r[1]+1))
        residues = _np.array(residues)
        residue_trans = _np.empty(_np.max(ranges)+1, dtype=int)
        for i in range(len(ranges)): residue_trans[ranges[i][0]:ranges[i][1]+1]=i
        if average: results_dict = {i:{j:_np.zeros(self.nb_frames) for j in range(len(ranges))} for i in range(len(ranges))}
        results_matrix = _np.zeros((len(ranges), len(ranges)), dtype=float)
        for key in results:
            _, resna, resida, _, resnb, residb = _hf.deconst_key(key, self.residuewise)
            if resna not in ['TIP3', 'HOH'] and resnb not in ['TIP3', 'HOH'] and resida in residues and residb in residues:
                if average: results_dict[residue_trans[resida]][residue_trans[residb]]+=results[key]
                else: results_matrix[residue_trans[resida], residue_trans[residb]] +=1
        if average:
            for i in range(len(ranges)):
                for j in range(len(ranges)):
                    results_matrix[i, j] = results_dict[i][j].mean()
        results_matrix = results_matrix + _np.tril(results_matrix.T, -1)
        #results_matrix = _np.rot90(results_matrix, 3)
        fig, ax = _plt.subplots()
        _plt.imshow(results_matrix, cmap='Reds')
        #_plt.colorbar(ticks=range(int(_np.max(results_matrix))+1))
        _plt.colorbar()
        ax.set_ylim(bottom=-0.5, top=len(names)-0.5)
        ax.set_xticks(_np.arange(len(names)))
        ax.set_yticks(_np.arange(len(names)))
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_xticklabels(names)
        #names.reverse()
        ax.set_yticklabels(names)
        result_string = 'labels ' + ' '.join(names) + '\n'
        result_string += 'resid_from_to ' + ' '.join([str(r[0])+'-'+str(r[1]) for r in ranges]) + '\n'
        result_string = _hf.string_in_columns(result_string) + ' \n\nblock_1 block_2 count\n'
        for i,row in enumerate(results_matrix):
            for j in range(i,len(results_matrix)):
                result_string += names[i] + ' ' + names[j] + ' ' + str(results_matrix[i,j].round(2)) + '\n'
        _plt.setp(ax.get_yticklabels(), rotation=90, ha="center", va="center", rotation_mode="anchor")
        _plt.xlabel(label , fontsize = 16)
        _plt.ylabel(label , fontsize = 16)
        if return_figure:
            _plt.close()
            return fig, result_string
        self._save_or_draw(filename, return_figure=return_figure)

    def draw_residue_residue_heatmap(self, average=False, use_filtered=True, filename=None, return_figure=False):
        if use_filtered: results = self.filtered_results
        else: results = self.initial_results
        residues = _np.sort(_np.unique([_hf.deconst_key(key, self.residuewise)[2] for key in results]+[_hf.deconst_key(key, self.residuewise)[5] for key in results]))
        res_trans = {r:i for i,r in enumerate(residues)}
        results_matrix = _np.zeros((len(residues), len(residues)), dtype=float)
        if average: results_average = _np.zeros((len(residues), len(residues), self.nb_frames))
        for key in results:
            _, resna, resida, _, resnb, residb = _hf.deconst_key(key, self.residuewise)
            if resida in residues and residb in residues:
                resmin, resmax = min(res_trans[resida],res_trans[residb]), max(res_trans[resida],res_trans[residb])
                if average: results_average[resmin, resmax]+=results[key]
                else: results_matrix[resmin, resmax] += 1
        if average:
            results_matrix = results_average.mean(2).round(2)
        results_matrix = results_matrix + _np.tril(results_matrix.T, -1)
        if not average: results_matrix = results_matrix.astype(_np.int)
        fig, ax = _plt.subplots()
        _plt.imshow(results_matrix, cmap='Reds', origin='lower')
        _plt.colorbar(ticks=range(int(_np.max(results_matrix))+1))
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        _plt.title('Save data for actual residue ids.')
        _plt.xlabel('Internal Residue Numbering' , fontsize = 16)
        _plt.ylabel('Internal Residue Numbering' , fontsize = 16)
        names = _np.array(residues+self.add_missing_residues, dtype=_np.str)
        result_string = 'bonded_residue_1 bonded_residue_2 count\n'
        ida, idb = _np.nonzero(results_matrix)
        for i in range(len(ida)):
            if ida[i]>idb[i]:continue
            result_string += names[ida[i]]+ ' ' + names[idb[i]] + ' ' + str(results_matrix[ida[i],idb[i]]) + '\n'
        if return_figure:
            _plt.close()
            return fig, result_string
        self._save_or_draw(filename, return_figure=return_figure)

    def draw_graph(self, draw_edge_occupancy=True, compare_to=None, mutations=None, highlight_interregion=False, in_frame=None, centrality=None, max_centrality=None, default_color='seagreen', node_factor=1.0, draw_labels=True, color_dict={}, filename=None, return_figure=False, draw_base=False, use_filtered=True):
        if compare_to is not None:
            if use_filtered:
                graph1 = self.filtered_graph
                graph2 = compare_to.filtered_graph
                results = self.filtered_results
                compare_results = compare_to.filtered_results
            else:
                graph1 = self.initial_graph
                graph2 = compare_to.initial_graph
                results = self.initial_results
                compare_results = compare_to.initial_results
            if mutations is not None:
                for mutation in mutations:
                    old_node, new_node = mutation.split('_')
                    temp_edges = [(new_node, connected_node) for connected_node in graph1[old_node]]
                    graph1.remove_node(old_node)
                    graph1.add_edges_from(temp_edges)
            composed_graph = _nx.compose(graph1, graph2)
            e_new, e_gone = [], []
            for a,b in graph1.edges():
                if not (a,b) in graph2.edges():
                    e_gone.append((a,b))
            for a,b in graph2.edges():
                if not (a,b) in graph1.edges():
                    e_new.append((a,b))
            edges_gone, edges_new = _nx.Graph(e_gone), _nx.Graph(e_new)
        else:
            if use_filtered:
                composed_graph = self.filtered_graph
                results = self.filtered_results
            else:
                composed_graph = self.initial_graph
                results = self.initial_results
        pos = []
        """
        for node in composed_graph.nodes():
            segid, _, resid = node.split('-')[:3]
            if hasattr(self, 'first_frame_dict'):
                if len(self.first_frame_dict) != 0: self._universe.trajectory[self.first_frame_dict[node]]
            atoms = self._universe.select_atoms("segid {} and resid {}".format(segid, resid))
            pos.append(_np.mean(atoms.positions, axis=0))
            nodes.append(node)
        """
        if in_frame is None:
            self.compute_joint_occupancy()
            if self.joint_occupancy_frames.size > 0: in_frame = self.joint_occupancy_frames[0]
            else: in_frame=0
        self._universe.trajectory[in_frame]
        all_coordinates = _np.vstack((self._da_selection.positions, self._water.positions, _np.array(self._ions.positions).reshape((-1,3))))
        nodes = _np.array(composed_graph.nodes())
        if self.residuewise: all_id = _np.array(self._all_ids+self._ions_ids)
        else: all_id = _np.array(self._all_ids_atomwise+self._ions_ids_atomwise)
        old_nodes, new_nodes = [], []
        if (mutations is not None) and (compare_to is not None):
            for mutation in mutations:
                old_node, new_node = mutation.split('_')
                old_nodes.append(old_node)
                new_nodes.append(new_node)
        for node in nodes:
            if node in new_nodes: node = old_nodes[new_nodes.index(node)]
            pos.append(all_coordinates[all_id == node].mean(0))
        pos = _np.array(pos)
        pos2d, base = _hf.pca_2d_projection(pos)

        _plt.figure(figsize=(_np.sqrt(len(pos))*2,_np.sqrt(len(pos))*1.5))
        if draw_base:
            x_min, y_max = _np.min(pos2d[:,0]), _np.max(pos2d[:,1])
            arrow_pos = base.dot(_np.eye(3).T).T * 1.4
            ax_label_pos = arrow_pos + arrow_pos / _np.linalg.norm(arrow_pos, axis=1).reshape(-1,1) * 0.5
            ax_label_pos[:,0] += x_min
            ax_label_pos[:,1] += y_max
            s = ['x', 'y', 'z']
            for i in range(3):
                _plt.arrow(x_min, y_max, arrow_pos[i,0], arrow_pos[i,1])
                _plt.text(ax_label_pos[i,0], ax_label_pos[i,1], s[i], horizontalalignment='center', verticalalignment='center')

        pos={}
        for p, node in zip(pos2d,nodes):
            pos[node] = p

        _nx.draw_networkx_edges(composed_graph, pos, width=1.5, alpha=0.5)

        if highlight_interregion:
            interregion_graph = _nx.Graph()
            for edge in composed_graph.edges():
                if edge[0].split('-')[0]!=edge[1].split('-')[0]:
                    interregion_graph.add_edge(*edge)
            _nx.draw_networkx_edges(interregion_graph, pos, width=10., alpha=0.7)
        if compare_to != None:
            _nx.draw_networkx_edges(edges_gone, pos, width=3.5, alpha=0.7, edge_color='red')
            _nx.draw_networkx_edges(edges_new, pos, width=3.5, alpha=0.7, edge_color='green')
        labels = {}
        color_graph = {default_color:_nx.Graph()}

        if color_dict:
            for color in color_dict.values():
                color_graph[color]=_nx.Graph()
        else:
            default_colors = ['seagreen', 'red', 'blue', 'yellow', 'grey', 'lightblue', 'orange', 'lightred']
            segnames = _np.unique(_np.array([[bond.split(':')[0].split('-')[0],bond.split(':')[1].split('-')[0]]for bond in results]))
            for i,resnm in enumerate(segnames):
                color_dict[resnm] = default_colors[i]
                color_graph[default_colors[i]] = _nx.Graph()

        for j, node in enumerate(composed_graph.nodes()):
            if self.residuewise: segname, resname, resid = node.split('-')
            else: segname, resname, resid, atomname = node.split('-')
            try:
                if self.residuewise:
                    if resname == 'TIP3': labels[node] = segname+str(int(resid)+self.add_missing_residues)
                    else: labels[node] = _hf.aa_three2one[resname]+str(int(resid)+self.add_missing_residues)
                else:
                    if resname == 'TIP3': labels[node] = segname+str(int(resid)+self.add_missing_residues)+'\n'+atomname
                    else: labels[node] = _hf.aa_three2one[resname]+str(int(resid)+self.add_missing_residues)+'\n'+atomname
            except KeyError:
                labels[node] = resname+resid
            try:color = color_dict[segname]
            except KeyError: color = default_color
            color_graph[color].add_node(node)
        #node_size=node_factor*_np.sqrt(_np.sqrt(_np.sqrt(len(pos))))*600,
        if centrality is not None:
            color_value_list = [centrality[node] if node in centrality else 0 for node in composed_graph.nodes()]
            if max_centrality is not None: mc = max_centrality
            else: mc = max(color_value_list)
            cmap = _plt.get_cmap('jet')
            _nx.draw_networkx_nodes(composed_graph, pos, node_size=900, alpha=0.5, node_color=color_value_list, cmap=cmap, vmin=0.0, vmax=mc)
            sm = _plt.cm.ScalarMappable(cmap=cmap)
            sm.set_array([0.0, mc])
            _plt.colorbar(sm, ax=_plt.axes())
        else:
            for color in color_graph:
                _nx.draw_networkx_nodes(color_graph[color], pos, node_size=900, alpha=0.5, node_color=color)
        #f = _np.sqrt(_np.sqrt(_np.sqrt(_np.sqrt(len(pos)))))*0.5*(1 - 0.5 * (1 - node_factor))
        f=0.6
        len2fontsize = _hf.defaultdict(lambda: 6*f, {2:18*f, 3:18*f, 4:16*f, 5:13*f, 6:11*f, 7:10*f, 8:9*f, 9:8*f, 10:7*f})
        for j, node in enumerate(composed_graph.nodes()):
            tempG = _nx.Graph()
            tempG.add_node(node)
            spl = labels[node].split('\n')
            if len(spl)>1: label_length = len(spl[0])+1
            else: label_length = len(spl[0])
            if draw_labels: _nx.draw_networkx_labels(tempG, pos, {node:labels[node]}, font_weight='bold' , font_size=len2fontsize[label_length])
        if draw_edge_occupancy:
            edge_labels = {}
            for conn in results:
                nodea, nodeb = conn.split(':')
                edge_labels[(nodea,nodeb)] = _np.round(results[conn].mean()*100, 1)

            if compare_to != None:
                compare_edge_labels = {}
                edge_labels_new = {}
                edge_labels_gone = {}
                for conn in compare_results:
                    nodea, nodeb = conn.split(':')
                    compare_edge_labels[(nodea,nodeb)] = _np.round(compare_results[conn].mean()*100, 1)

                for edge in list(edge_labels.keys()):
                    if edge in e_gone:
                        edge_labels_gone[edge] = edge_labels[edge]
                        del edge_labels[edge]

                for edge in list(compare_edge_labels.keys()):
                    if edge in e_new:
                        edge_labels_new[edge] = compare_edge_labels[edge]
                        del compare_edge_labels[edge]

                _nx.draw_networkx_edge_labels(composed_graph, pos, edge_labels=edge_labels, font_size=len2fontsize[4], label_pos=0.33, font_color='red')
                _nx.draw_networkx_edge_labels(composed_graph, pos, edge_labels=edge_labels_gone, font_size=len2fontsize[4], label_pos=0.5, font_color='red')
                _nx.draw_networkx_edge_labels(composed_graph, pos, edge_labels=compare_edge_labels, font_size=len2fontsize[4], label_pos=0.66, font_color='green')
                _nx.draw_networkx_edge_labels(composed_graph, pos, edge_labels=edge_labels_new, font_size=len2fontsize[4], label_pos=0.5, font_color='green')
            else:
                _nx.draw_networkx_edge_labels(composed_graph, pos, edge_labels=edge_labels, font_size=len2fontsize[4], label_pos=0.5)
        _plt.axis('off')
        if return_figure:
            fig = _plt.gcf()
            _plt.close()
            return fig, 'NoData'
        self._save_or_draw(filename, return_figure=return_figure)

    def _reload_universe(self):
        super(NetworkAnalysis, self)._reload_universe()
        sorted_selection = _hf.Selection(self._mda_selection, self.donor_names, self.acceptor_names, self._add_donors_without_hydrogen)
        if not sorted_selection.donors: da_selection = sorted_selection.acceptors
        elif not sorted_selection.acceptors: da_selection = sorted_selection.donors
        else: da_selection = _MDAnalysis.core.groups.AtomGroup(sorted_selection.donors + sorted_selection.acceptors)
        self._da_selection = da_selection
        if sorted_selection.donors: self._donors = _MDAnalysis.core.groups.AtomGroup(sorted_selection.donors)
        else: self._donors = _hf.EmptyGroup()
        if sorted_selection.acceptors: self._acceptors = _MDAnalysis.core.groups.AtomGroup(sorted_selection.acceptors)
        else: self._acceptors = _hf.EmptyGroup()
        water_hydrogen = [h for l in self._water for h in l.residue.atoms[1:]]
        if not sorted_selection.hydrogens and not water_hydrogen:
            if self.check_angle: raise AssertionError('There are no possible hbond donors in the selection and no water. Since check_angle is True, hydrogen is needed for the calculations!')
            else: hydrogen = _hf.EmptyGroup()
        elif not sorted_selection.hydrogens: hydrogen = _MDAnalysis.core.groups.AtomGroup(water_hydrogen)
        elif not water_hydrogen: hydrogen = sorted_selection.hydrogens
        else: hydrogen = sorted_selection.hydrogens + _MDAnalysis.core.groups.AtomGroup(water_hydrogen)
        self._hydrogen = hydrogen

    def _compute_graph_in_frame(self, frame, set_as_filtered_results=False, use_filtered=True):
        if use_filtered: results = self.filtered_results
        else: results = self.initial_results
        keep_edges=[]
        for key in results:
            if results[key][frame]:
                keep_edges.append((key.split(':')[0], key.split(':')[1]))
        graph = _nx.Graph()
        graph.add_edges_from(keep_edges)
        if set_as_filtered_results:
            self.filtered_graph = graph
            self._generate_filtered_results_from_filtered_graph()
        return graph

    def _generate_current_results_from_graph(self):
        temp_res = {}
        for resa, resb in self.initial_graph.edges():
            key, key_check = ':'.join((resa, resb)), ':'.join((resb, resa))
            try: temp_res[key] = self.initial_results[key]
            except KeyError:  temp_res[key_check] = self.initial_results[key_check]
        self.initial_results = temp_res

    def _generate_filtered_results_from_filtered_graph(self):
        temp_res = {}
        for resa, resb in self.filtered_graph.edges():
            key, key_check = ':'.join((resa, resb)), ':'.join((resb, resa))
            try: temp_res[key] = self.initial_results[key]
            except KeyError:  temp_res[key_check] = self.initial_results[key_check]
        self.filtered_results = temp_res

    def _generate_graph_from_current_results(self):
        self.initial_graph = _hf.dict2graph(self.initial_results, self.residuewise)

    def _generate_filtered_graph_from_filtered_results(self):
        self.filtered_graph = _hf.dict2graph(self.filtered_results, self.residuewise)
