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
from .network import NetworkAnalysis
import numpy as _np
import networkx as _nx
import MDAnalysis as _MDAnalysis
from scipy import spatial as _sp
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra
from collections import OrderedDict as _odict
from matplotlib.ticker import MaxNLocator
#import matplotlib
#matplotlib.use('TKAgg', warn=False)
import matplotlib.pyplot as _plt

class WireAnalysis(NetworkAnalysis):

    def __init__(self, selection=None, structure=None, trajectories=None, distance=3.5, cut_angle=60.,
                 start=None, stop=None, step=1, additional_donors=[],
                 additional_acceptors=[], exclude_donors=[], exclude_acceptors=[],
                 ions=[], check_angle=True, residuewise=True, add_donors_without_hydrogen=False, restore_filename=None):

        super(WireAnalysis, self).__init__(selection=selection, structure=structure, trajectories=trajectories,
             distance=distance, cut_angle=cut_angle, start=start, stop=stop, step=step,
             additional_donors=additional_donors, additional_acceptors=additional_acceptors,
             exclude_donors=exclude_donors, exclude_acceptors=exclude_acceptors,
             ions=ions, check_angle=check_angle, residuewise=residuewise,
             add_donors_without_hydrogen=add_donors_without_hydrogen, restore_filename=restore_filename)

        if restore_filename != None: return
        if not self._mda_selection:  raise AssertionError('No atoms match the selection')
        sorted_selection = _hf.Selection(self._mda_selection, self.donor_names, self.acceptor_names, add_donors_without_hydrogen=add_donors_without_hydrogen)
        if not sorted_selection.donors: da_selection = sorted_selection.acceptors
        elif not sorted_selection.acceptors: da_selection = sorted_selection.donors
        else: da_selection = _MDAnalysis.core.groups.AtomGroup(sorted_selection.donors + sorted_selection.acceptors)
        da_ids = _hf.MDA_info_list(da_selection)
        self.hashs = {}
        self.hash_table = {}
        self.wire_lengths = {}
        da_u, da_ind, da_inv = _np.unique(da_ids, return_index=True, return_inverse=True)
        self.da_trans = da_ind[da_inv]

    def set_water_wires_csr(self, max_water=5, allow_direct_bonds = True, water_in_convex_hull=False):

        distances = {}
        path_hashs = {}
        frame_count = 0
        hash_table = {}
        no_direct_bonds = False
        self._allow_direct_bonds = allow_direct_bonds
        for ts in self._universe.trajectory[self._trajectory_slice]:

            water_coordinates = self._water.positions
            selection_coordinates = self._da_selection.positions
            selection_tree = _sp.cKDTree(selection_coordinates)
            try:
                d_tree = _sp.cKDTree(self._donors.positions)
                a_tree = _sp.cKDTree(self._acceptors.positions)
            except:
                no_direct_bonds = True
            hydrogen_coordinates = self._hydrogen.positions

            water_tree = _sp.cKDTree(water_coordinates, leafsize=32)
            local_water_index = []
            [local_water_index.extend(l) for l in water_tree.query_ball_point(selection_coordinates, float(max_water+1)*self.distance/2.)]
            local_water_index = _np.unique(local_water_index)
            local_water_coordinates = water_coordinates[local_water_index]

            if water_in_convex_hull:
                hull = _sp.Delaunay(selection_coordinates)
                local_water_index_hull = (hull.find_simplex(local_water_coordinates) != -1).nonzero()[0]
                local_water_coordinates = water_coordinates[local_water_index[local_water_index_hull]]

            local_water_tree = _sp.cKDTree(local_water_coordinates)

            local_water_index += self._first_water_id
            local_pairs = [(i, local_water_index[j]) for i, bla in enumerate(selection_tree.query_ball_tree(local_water_tree, self.distance)) for j in bla]
            local_water_index -= self._first_water_id
            water_pairs = [(local_water_index[p[0]], local_water_index[p[1]]) for p in local_water_tree.query_pairs(self.distance)]
            if not no_direct_bonds:
                da_pairs = _np.array([[i, j] for i,donors in enumerate(a_tree.query_ball_tree(d_tree, self.distance)) for j in donors])
                da_pairs[:,0] += self._nb_donors
            else:
                da_pairs = _np.array([])
            if self.check_angle:
                all_coordinates = _np.vstack((selection_coordinates, water_coordinates))
                da_hbonds = _hf.check_angle(da_pairs, self.heavy2hydrogen, all_coordinates, hydrogen_coordinates, self.cut_angle)
                water_hbonds = _hf.check_angle_water(water_pairs, water_coordinates, hydrogen_coordinates[self._first_water_hydrogen_id:], self.cut_angle)
                local_hbonds = _hf.check_angle(local_pairs, self.heavy2hydrogen, all_coordinates, hydrogen_coordinates, self.cut_angle)
            else:
                da_hbonds = da_pairs
                water_hbonds = water_pairs
                local_hbonds = local_pairs
            da_hbonds = _np.sort(da_hbonds)
            water_hbonds = _np.sort(water_hbonds) + self._first_water_id
            local_hbonds = _np.sort(_np.array(local_hbonds))
            local_hbonds[:,0]=self.da_trans[local_hbonds[:,0]]
            if no_direct_bonds: hbonds = _np.vstack((local_hbonds, water_hbonds))
            else: hbonds = _np.vstack((_np.vstack((da_hbonds,local_hbonds)), water_hbonds))
            no_direct_bonds = False
            water_da = _np.zeros(len(hbonds), dtype=bool)
            water_da[len(da_hbonds):len(da_hbonds)+len(local_hbonds)]=True
            water_water = _np.zeros(len(hbonds), dtype=bool)
            water_water[len(da_hbonds)+len(local_hbonds):] = True
            uniques, rowsncols = _np.unique(hbonds, return_inverse=True)
            rows, cols = rowsncols.reshape(hbonds.shape).T
            nb_nodes = uniques.size
            residues = (uniques < self._first_water_id).nonzero()[0]
            data = _np.ones(len(hbonds), dtype=_np.float)
            g = csr_matrix((data, (rows, cols)), shape=(nb_nodes, nb_nodes))
            local_rows, local_cols = rows[water_da], cols[water_da]
            g[local_rows, local_cols] = _np.inf
            already_checked = []
            for source in residues:
                source_index = local_rows == source
                if source_index.sum()==0: continue
                g[local_rows[source_index], local_cols[source_index]] = 1.0
                lengths, predecessors = dijkstra(g, directed=False, indices=source, unweighted=False, limit=max_water, return_predecessors=True)
                g[local_rows[source_index], local_cols[source_index]] = _np.inf
                in_range = _np.nonzero(lengths <= max_water)[0]
                target_index = _np.in1d(local_cols,in_range)
                targets = _np.unique(local_rows[target_index])
                for target in targets:
                    if target in already_checked or target==source: continue
                    target_water = local_cols[(local_rows==target) & target_index]
                    length_index = _np.argmin(lengths[target_water])
                    length = lengths[target_water][length_index]
                    wire = uniques[[_hf.predecessor_recursive_1d(ii, predecessors, target_water[length_index]) for ii in range(int(length))[::-1]]+[target]]
                    ai, bi = _np.sort(uniques[[source, target]])
                    wire_hash = hash(wire.tobytes())
                    hash_table[wire_hash] = wire
                    wire_info = self._all_ids[ai]+':'+self._all_ids[bi]
                    check = self._all_ids[bi]+':'+self._all_ids[ai]
                    if check in distances: wire_info=check
                    try:
                        distances[wire_info][frame_count] = length
                        path_hashs[wire_info][frame_count] = wire_hash
                    except KeyError:
                        distances[wire_info] = _np.ones(self.nb_frames, dtype=int) * _np.inf
                        path_hashs[wire_info] = _np.arange(self.nb_frames, dtype=int)
                        distances[wire_info][frame_count] = length
                        path_hashs[wire_info][frame_count] = wire_hash
                already_checked.append(source)

            if allow_direct_bonds:
                for a, b in da_hbonds:
                    wire_info = self._all_ids[a]+':'+self._all_ids[b]
                    check = self._all_ids[b]+':'+self._all_ids[a]
                    if check in distances: wire_info=check
                    try:
                        path_hashs[wire_info][frame_count] = -1
                        distances[wire_info][frame_count] = 0
                    except:
                        distances[wire_info] = _np.ones(self.nb_frames, dtype=_np.int)*_np.inf
                        distances[wire_info][frame_count] = 0
                        path_hashs[wire_info] = _np.arange(self.nb_frames, dtype=_np.int)
                        path_hashs[wire_info][frame_count] = -1
            frame_count += 1
        self._set_results({key:distances[key]!=_np.inf for key in distances})
        self.wire_lengths = distances
        self.hashs = path_hashs
        self.hash_table = hash_table
        self.occupancy_dict = {}
        self.first_frame_dict = {}


    def set_water_wires(self, max_water=5, allow_direct_bonds=True, water_in_convex_hull=False):

        intervals_results = {}
        results = {}
        frame_count = 0
        frames = self.nb_frames
        this_frame_table = {}
        no_direct_bonds = False
        self._allow_direct_bonds = allow_direct_bonds
        for ts in self._universe.trajectory[self._trajectory_slice]:

            water_coordinates = self._water.positions
            selection_coordinates = self._da_selection.positions
            water_tree = _sp.cKDTree(water_coordinates, leafsize=32)
            selection_tree = _sp.cKDTree(selection_coordinates)
            try:
                d_tree = _sp.cKDTree(self._donors.positions)
                a_tree = _sp.cKDTree(self._acceptors.positions)
            except:
                no_direct_bonds = True
            hydrogen_coordinates = self._hydrogen.positions

            local_water_index = []
            [local_water_index.extend(l) for l in water_tree.query_ball_point(selection_coordinates, float(max_water+1)*self.distance/2.)]
            local_water_index = _np.unique(local_water_index)
            local_water_coordinates = water_coordinates[local_water_index]

            if water_in_convex_hull:
                hull = _sp.Delaunay(selection_coordinates)
                local_water_index_hull = (hull.find_simplex(local_water_coordinates) != -1).nonzero()[0]
                local_water_coordinates = water_coordinates[local_water_index[local_water_index_hull]]

            local_water_tree = _sp.cKDTree(local_water_coordinates)

            local_water_index += self._first_water_id
            local_pairs = [(i, local_water_index[j]) for i, bla in enumerate(selection_tree.query_ball_tree(local_water_tree, self.distance)) for j in bla]
            local_water_index -= self._first_water_id
            water_pairs = local_water_index[_np.array(list(local_water_tree.query_pairs(self.distance)))]
            if not no_direct_bonds:
                da_pairs = _np.array([[i, j] for i,donors in enumerate(a_tree.query_ball_tree(d_tree, self.distance)) for j in donors])
                da_pairs[:,0] += self._nb_donors
            else:
                da_pairs = _np.array([])
                no_direct_bonds = False
            if self.check_angle:
                all_coordinates = _np.vstack((selection_coordinates, water_coordinates))
                da_hbonds = _hf.check_angle(da_pairs, self.heavy2hydrogen, all_coordinates, hydrogen_coordinates, self.cut_angle)
                water_hbonds = _hf.check_angle_water(water_pairs, water_coordinates, hydrogen_coordinates[self._first_water_hydrogen_id:], self.cut_angle)
                local_hbonds = _hf.check_angle(local_pairs, self.heavy2hydrogen, all_coordinates, hydrogen_coordinates, self.cut_angle)
            else:
                da_hbonds = da_pairs
                water_hbonds = water_pairs
                local_hbonds = local_pairs

            da_hbonds = _np.sort(da_hbonds)
            water_hbonds = _np.sort(water_hbonds) + self._first_water_id
            local_hbonds = _np.sort(_np.array(local_hbonds))
            local_hbonds[:,0]=self.da_trans[local_hbonds[:,0]]

            g = _nx.Graph()
            g.add_edges_from(water_hbonds)

            residues = _np.unique(local_hbonds[:,0])
            already_checked=[]
            for source in residues:
                already_checked_targets = []
                source_water_index = local_hbonds[:,0]==source

                if not source_water_index.any(): continue

                g.add_edges_from(local_hbonds[source_water_index])
                paths = _nx.single_source_shortest_path(g,source,max_water)
                g.remove_node(source)

                idx = _np.array([self._all_ids[bl]!=self._all_ids[source] for bl in local_hbonds[:,0]])
                target_water_set = set(paths) & set(local_hbonds[:,1][idx])
                twlist = list(target_water_set)

                all_targets_index = _np.in1d(local_hbonds[:,1], _np.array(twlist))

                for target, last_water in local_hbonds[all_targets_index]:

                    if target in already_checked or target in already_checked_targets or self._all_ids[source]==self._all_ids[target]: continue
                    wire = paths[last_water] + [target]
                    wire_hash = hash(str(wire))
                    this_frame_table[wire_hash] = wire
                    wire_info = self._all_ids[source]+':'+self._all_ids[target]
                    check = self._all_ids[target]+':'+self._all_ids[source]
                    if check in results: wire_info=check
                    water_in_wire = len(wire)-2

                    try: water_already_found = results[wire_info][frame_count]
                    except: water_already_found = _np.inf
                    if (water_in_wire >= water_already_found): continue

                    try:
                        results[wire_info][frame_count] = water_in_wire
                        intervals_results[wire_info][frame_count] = wire_hash
                    except:
                        results[wire_info] = _np.ones(frames)*_np.inf
                        results[wire_info][frame_count] = water_in_wire
                        intervals_results[wire_info] = _np.arange(frames, dtype=_np.int)
                        intervals_results[wire_info][frame_count] = wire_hash

                already_checked.append(source)

            if allow_direct_bonds:
                for a, b in da_hbonds:
                    wire_info = self._all_ids[a]+':'+self._all_ids[b]
                    check = self._all_ids[b]+':'+self._all_ids[a]
                    if check in results: wire_info=check
                    try:
                        intervals_results[wire_info][frame_count] = -1
                        results[wire_info][frame_count] = 0
                    except:
                        results[wire_info] = _np.ones(frames, dtype=_np.int)*_np.inf
                        results[wire_info][frame_count] = 0
                        intervals_results[wire_info] = _np.arange(frames, dtype=_np.int)
                        intervals_results[wire_info][frame_count] = -1

            frame_count += 1
        self._set_results({key:results[key]!=_np.inf for key in results})
        self.wire_lengths = results
        self.hashs = intervals_results
        self.hash_table = this_frame_table
        self.occupancy_dict = {}
        self.first_frame_dict = {}

    def set_explicit_water_wires(self, max_water=5, allow_direct_bonds=True, water_in_convex_hull=False):

        self.set_water_wires(max_water, allow_direct_bonds, water_in_convex_hull)

        filtered_results = {}
        for wire_key in self.hashs:
            hashs, hash_index, hash_count = _np.unique(self.hashs[wire_key], return_index=True, return_counts=True)
            filter_occupancy = _np.in1d(hashs, _np.arange(-1, self.nb_frames, dtype=int), invert=True)
            occupancy = (hash_count / self.nb_frames)[filter_occupancy]
            for i,h in enumerate(hashs[filter_occupancy]):
                try:
                    w_string = ':'.join([self._all_ids[w_id] for w_id in self.hash_table[h]])
                except KeyError:
                    if allow_direct_bonds: w_string = wire_key
                    else: break
                try:
                    filtered_results[occupancy[i]].append((w_string, hash_index[filter_occupancy][i], self.hashs[wire_key] == h))
                except:
                    filtered_results[occupancy[i]] = []
                    filtered_results[occupancy[i]].append((w_string, hash_index[filter_occupancy][i], self.hashs[wire_key] == h))
        occupancy_dict = _odict(sorted(filtered_results.items()))

        edges = []
        first_frame_dict = {}
        temp_occupancy_dict = {}
        for occupancy in occupancy_dict:
            for wire, first_appearance_frame, occupancy_series in occupancy_dict[occupancy]:
                nodes_in_wire = wire.split(':')
                for node in nodes_in_wire: first_frame_dict[node] = first_appearance_frame
                for nodea, nodeb in _hf.pairwise(nodes_in_wire):
                    edges.append(':'.join((nodea,nodeb)))
                    temp_occupancy_dict[':'.join((nodea,nodeb))] = occupancy_series

        self.current_results = temp_occupancy_dict
        self.filtered_results = temp_occupancy_dict
        self._generate_graph_from_current_results()
        self._generate_filtered_graph_from_filtered_results()
        self.occupancy_dict = occupancy_dict
        self.first_frame_dict = first_frame_dict

    def compute_average_water_per_wire(self, use_filtered=True):
        if use_filtered: results = self.filtered_results
        else: results = self.initial_results
        return {key:self.wire_lengths[key][self.wire_lengths[key]<_np.inf].mean() for key in results}

    def filter_minmal_bonds_path(self, start, goal, use_filtered=True):
        if use_filtered: graph = self.filtered_graph
        else: graph = self.initial_graph
        if len(graph.nodes()) == 0: raise AssertionError('nothing to filter!')
        if start not in graph.nodes(): raise AssertionError('The start node is not in the graph')
        if goal not in graph.nodes(): raise AssertionError('The goal node is not in the graph')
        for component in _nx.connected_component_subgraphs(graph):
            if start in component.nodes(): break
        if goal in component:
            weighted_component = _nx.Graph()
            avg = self.compute_average_water_per_wire(use_filtered=use_filtered)
            triplets = []
            for u, v in component.edges():
                try:
                    triplets.append((u,v,avg[u+':'+v]))
                except KeyError:
                    triplets.append((u,v,avg[v+':'+u]))
            weighted_component.add_weighted_edges_from(triplets)
            node_set = _nx.all_shortest_paths(weighted_component, start, goal, weight='weight')
        else: raise AssertionError('start and goal nodes are not connected')
        path_graph = _nx.Graph()
        for path in node_set:
            path_graph.add_edges_from(_hf.pairwise(path))
        self.applied_filters['avg_least_bonds']=(start, goal)
        self.filtered_graph = path_graph
        self._generate_filtered_results_from_filtered_graph()

    def compute_all_shortest_paths_info(self, start, goal, use_filtered=True):
        if use_filtered:
            graph = self.filtered_graph
            results = self.filtered_results
        else:
            graph = self.initial_graph
            results = self.initial_results
        if len(graph.nodes()) == 0: raise AssertionError('nothing to filter!')
        if start not in graph.nodes(): raise AssertionError('The start node is not in the graph')
        if goal not in graph.nodes(): raise AssertionError('The goal node is not in the graph')
        try:
            paths = _nx.all_shortest_paths(graph, start, goal)
        except:
            raise AssertionError('start and goal nodes are not connected')
        res = {}
        res_jo = {}
        jos = {}
        avg = self.compute_average_water_per_wire(use_filtered=use_filtered)
        for path in paths:
            key = start+'-'+'-'.join([node.split('-')[2] for node in path[1:-1]])+'-'+goal
            val = _np.sum([avg[a+':'+b] if a+':'+b in avg else avg[b+':'+a] for a,b in _hf.pairwise(path)])
            res[key] = val
            jo = _np.ones(self.nb_frames, dtype=bool)
            jos[key] = jo
            for a,b in _hf.pairwise(path): jo &= results[a+':'+b] if a+':'+b in results else results[b+':'+a]
            res_jo[key] = jo.mean()
        return res, res_jo, jos

    def compute_wire_projection(self):
        frame_count = -1
        projection_results = {}
        for ts in self._universe.trajectory[self._trajectory_slice]:
            frame_count += 1
            for rmsd_key in self.hashs:
                rmsd_hash = self.hashs[rmsd_key][frame_count]
                if rmsd_hash == frame_count or rmsd_hash == -1: continue
                rmsd_wire = self.hash_table[rmsd_hash]
                rmsd_water = rmsd_wire[1:-1]
                rmsd_residues = [rmsd_wire[0], rmsd_wire[-1]]
                nb_water = len(rmsd_water)
                water_coords = self._water.positions
                da_coords = self._da_selection.positions
                all_coordinates = _np.vstack((da_coords, water_coords))
                c_a, c_b, c_c = _np.repeat(all_coordinates[rmsd_residues[0]].reshape((-1,3)), nb_water, axis=0), _np.repeat(all_coordinates[rmsd_residues[1]].reshape((-1,3)), nb_water, axis=0), all_coordinates[rmsd_water].reshape((-1,3))
                angles = _hf.angle(c_a, c_b, c_c)
                dists = _np.sqrt(((c_b - c_c)**2).sum(axis=1))
                proj_sum = (_np.sin(_np.deg2rad(angles)) * dists).mean()
                try:
                    projection_results[rmsd_key][frame_count] = proj_sum
                except:
                    projection_results[rmsd_key] = _np.ones(self.nb_frames)*-1
                    projection_results[rmsd_key][frame_count] = proj_sum
        return projection_results

    def draw_water_timeseries(self, resa, resb, scatter_size=0.5, filename=None, return_figure=False):
        try:
            timeseries = self.wire_lengths[resa+':'+resb]
        except KeyError:
            timeseries = self.wire_lengths[resb+':'+resa]
        fig, ax = _plt.subplots()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        _plt.scatter(range(len(timeseries)), timeseries, s=scatter_size)
        _plt.xlabel('frame' , fontsize = 16)
        _plt.ylabel('# water' , fontsize = 16)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        timeseries[timeseries==_np.inf]=-1
        return_str = resa+':'+resb +' '+ ' '.join(timeseries.astype(_np.int).astype(_np.str))
        if return_figure:
            _plt.close()
            return fig, _hf.string_in_columns(return_str)
        self._save_or_draw(filename, return_figure=return_figure)

    def _save_average_and_occupancy_to_one_file(self, filename):
        avg_water = self.compute_average_water_per_wire(use_filtered=True)
        string = ''
        for key in self.filtered_results:
            string += key + ' ' + str(self.filtered_results[key].mean()) + ' ' + str(avg_water[key]) +'\n'
        with open(filename, 'w') as file:
            file.write(string)
