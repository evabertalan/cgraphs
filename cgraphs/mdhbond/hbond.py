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
from itertools import combinations
from matplotlib.ticker import MaxNLocator
#import matplotlib
#matplotlib.use('TKAgg', warn=False)
import matplotlib.pyplot as _plt


_np.int = int

class HbondAnalysis(NetworkAnalysis):

    def __init__(self, selection=None, structure=None, trajectories=None, distance=3.5, cut_angle=60.,
                 start=None, stop=None, step=1, additional_donors=[],
                 additional_acceptors=[], exclude_donors=[], exclude_acceptors=[],
                 ions=[], check_angle=True, residuewise=True, add_donors_without_hydrogen=False, restore_filename=None):

        super(HbondAnalysis, self).__init__(selection=selection, structure=structure, trajectories=trajectories,
             distance=distance, cut_angle=cut_angle, start=start, stop=stop, step=step,
             additional_donors=additional_donors, additional_acceptors=additional_acceptors,
             exclude_donors=exclude_donors, exclude_acceptors=exclude_acceptors,
             ions=ions, check_angle=check_angle, residuewise=residuewise,
             add_donors_without_hydrogen=add_donors_without_hydrogen, restore_filename=restore_filename)

        if restore_filename is not None: return
        self._i4_distribution = None

    def set_add_ion_interactions(self, include_water=True):
        if not self._ions: raise AssertionError('No ions were specified during initialization.')
        frame_count = 0
        frames = self.nb_frames
        result = {}
        if self.residuewise:
            ion_ids = self._ions_ids
            water_ids = self._water_ids
        else:
            ion_ids = self._ions_ids_atomwise
            water_ids = self._water_ids_atomwise

        for ts in self._universe.trajectory[self._trajectory_slice]:
            selection_coordinates = self._da_selection.positions
            ion_coordinates = self._ions.positions
            selection_tree = _sp.cKDTree(selection_coordinates)
            ion_tree = _sp.cKDTree(ion_coordinates)
            da_ion_pairs = _np.array([[i, j] for i,ions in enumerate(selection_tree.query_ball_tree(ion_tree, self.distance)) for j in ions])
            frame_res = [self._all_ids[i]+':'+ion_ids[j] for i,j in da_ion_pairs]

            if include_water:
                water_coordinates = self._water.positions
                water_tree = _sp.cKDTree(water_coordinates, leafsize=32)
                water_ion_pairs = _np.array([[i, j] for i,ions in enumerate(water_tree.query_ball_tree(ion_tree, self.distance)) for j in ions])
                water_frame_res = [water_ids[i]+':'+ion_ids[j] for i,j in water_ion_pairs]
                frame_res += water_frame_res

            for bond in frame_res:
                try:
                    result[bond][frame_count] = True
                except:
                    result[bond] = _np.zeros(frames, dtype=bool)
                    result[bond][frame_count] = True

            frame_count+=1

        if not self.initial_results: self._set_results(result)
        else: self._add_overwrite_results(result)

    def set_hbonds_in_selection(self, exclude_backbone_backbone=True):
        if exclude_backbone_backbone: backbone_filter = _np.array([(ids.split('-')[3] in ['O', 'N']) for ids in self._all_ids_atomwise])
        frame_count = 0
        frames = self.nb_frames
        result = {}
        if self._nb_acceptors == 0: raise AssertionError('No acceptors in the selection')
        if self._nb_donors == 0: raise AssertionError('No donors in the selection')
        for ts in self._universe.trajectory[self._trajectory_slice]:
            selection_coordinates = self._da_selection.positions
            d_tree = _sp.cKDTree(self._donors.positions)
            a_tree = _sp.cKDTree(self._acceptors.positions)
            hydrogen_coordinates = self._hydrogen.positions

            da_pairs = _np.array([[i, j] for i,donors in enumerate(a_tree.query_ball_tree(d_tree, self.distance)) for j in donors])
            da_pairs[:,0] += self._nb_donors
            if exclude_backbone_backbone: da_pairs = da_pairs[_np.logical_not(_np.all(backbone_filter[da_pairs], axis=1))]

            if self.check_angle:
                all_coordinates = selection_coordinates
                local_hbonds = _hf.check_angle(da_pairs, self.heavy2hydrogen, all_coordinates, hydrogen_coordinates, self.cut_angle)
            else:
                local_hbonds = da_pairs

            sorted_bonds = _np.sort(local_hbonds)
            check = self._resids[sorted_bonds]
            check = check[:,0] < check[:,1]
            frame_res = [self._all_ids[i] + ':' + self._all_ids[j] if check[ii] else self._all_ids[j] + ':' + self._all_ids[i] for ii, (i, j) in enumerate(sorted_bonds)]

            for bond in frame_res:
                if not self.check_angle:
                    a, b = bond.split(':')
                    if a.split('-')[:3] == b.split('-')[:3]: continue
                try:
                    result[bond][frame_count] = True
                except:
                    result[bond] = _np.zeros(frames, dtype=bool)
                    result[bond][frame_count] = True
            frame_count+=1
        self._set_results(result)
        if self._ions_list: self.set_add_ion_interactions(False)

    def set_hbonds_only_water_in_convex_hull(self):
        result = {}
        frame_count = 0
        frames = self.nb_frames

        for ts in self._universe.trajectory[self._trajectory_slice]:
            water_coordinates = self._water.positions
            select_coordinates = self._da_selection.positions
            hydrogen_coordinates = self._hydrogen.positions[self._first_water_hydrogen_id:]

            hull = _sp.Delaunay(select_coordinates)
            local_index = (hull.find_simplex(water_coordinates) != -1).nonzero()[0]

            local_water_coordinates = water_coordinates[local_index]
            local_water_tree = _sp.cKDTree(local_water_coordinates)
            water_pairs = local_index[_np.array([pair for pair in local_water_tree.query_pairs(self.distance)])]

            if self.check_angle:
                local_hbonds = _hf.check_angle_water(water_pairs, water_coordinates, hydrogen_coordinates, self.cut_angle)
            else:
                local_hbonds = water_pairs

            sorted_bonds = _np.sort(local_hbonds)
            check = self._resids[sorted_bonds]
            check = check[:,0] < check[:,1]
            if self.residuewise: frame_res = [self._all_ids[i] + ':' + self._all_ids[j] if check[ii] else self._all_ids[j] + ':' + self._all_ids[i] for ii, (i, j) in enumerate(sorted_bonds)]
            else: frame_res = [self._all_ids_atomwise[i] + ':' + self._all_ids_atomwise[j] if check[ii] else self._all_ids_atomwise[j] + ':' + self._all_ids_atomwise[i] for ii, (i, j) in enumerate(sorted_bonds)]

            for bond in frame_res:
                if not self.check_angle:
                    a, b = bond.split(':')
                    if a.split('-')[:3] == b.split('-')[:3]: continue
                try:
                    result[bond][frame_count] = True
                except:
                    result[bond] = _np.zeros(frames, dtype=bool)
                    result[bond][frame_count] = True
            frame_count+=1
        self._set_results(result)

    def set_hbonds_in_selection_and_water_in_convex_hull(self, exclude_backbone_backbone=True):
        if exclude_backbone_backbone: backbone_filter = _np.array([(ids.split('-')[3] in ['O', 'N']) for ids in self._all_ids_atomwise])
        result = {}
        frame_count = 0
        frames = self.nb_frames
        for ts in self._universe.trajectory[self._trajectory_slice]:
            water_coordinates = self._water.positions
            select_coordinates = self._da_selection.positions
            hydrogen_coordinates = self._hydrogen.positions

            selection_tree = _sp.cKDTree(select_coordinates)
            if self._nb_acceptors > 0 and self._nb_donors > 0:
                d_tree = _sp.cKDTree(self._donors.positions)
                a_tree = _sp.cKDTree(self._acceptors.positions)

            hull = _sp.Delaunay(select_coordinates)
            local_water_index = (hull.find_simplex(water_coordinates) != -1).nonzero()[0]

            local_water_coordinates = water_coordinates[local_water_index]
            local_water_tree = _sp.cKDTree(local_water_coordinates)

            local_pairs = [(i, local_water_index[j]+self._first_water_id) for i, bla in enumerate(selection_tree.query_ball_tree(local_water_tree, self.distance)) for j in bla]
            water_pairs = [(local_water_index[p[0]]+self._first_water_id, local_water_index[p[1]]+self._first_water_id) for p in local_water_tree.query_pairs(self.distance)]
            if self._nb_acceptors > 0 and self._nb_donors > 0:
                da_pairs = _np.array([[i, j] for i,donors in enumerate(a_tree.query_ball_tree(d_tree, self.distance)) for j in donors])
                da_pairs[:,0] += self._nb_donors
                if exclude_backbone_backbone: da_pairs = da_pairs[_np.logical_not(_np.all(backbone_filter[da_pairs], axis=1))]
            else: da_pairs = []

            if self.check_angle:
                all_coordinates = _np.vstack((select_coordinates, water_coordinates))
                hbonds = _hf.check_angle(list(da_pairs)+water_pairs+local_pairs, self.heavy2hydrogen, all_coordinates, hydrogen_coordinates, self.cut_angle)
            else:
                hbonds = list(da_pairs) + water_pairs + local_pairs

            sorted_bonds = _np.sort(hbonds)
            check = self._resids[sorted_bonds]
            check = check[:,0] < check[:,1]
            if self.residuewise: frame_res = [self._all_ids[i] + ':' + self._all_ids[j] if check[ii] else self._all_ids[j] + ':' + self._all_ids[i] for ii, (i, j) in enumerate(sorted_bonds)]
            else: frame_res = [self._all_ids_atomwise[i] + ':' + self._all_ids_atomwise[j] if check[ii] else self._all_ids_atomwise[j] + ':' + self._all_ids_atomwise[i] for ii, (i, j) in enumerate(sorted_bonds)]

            for bond in frame_res:
                if not self.check_angle:
                    a, b = bond.split(':')
                    if a.split('-')[:3] == b.split('-')[:3]: continue
                try:
                    result[bond][frame_count] = True
                except:
                    result[bond] = _np.zeros(frames, dtype=bool)
                    result[bond][frame_count] = True
            frame_count+=1
        self._set_results(result)
        if self._ions_list: self.set_add_ion_interactions(True)


    def set_hbonds_in_selection_and_water_around(self, around_radius, exclude_backbone_backbone=True):
        if exclude_backbone_backbone: backbone_filter = _np.array([(ids.split('-')[3] in ['O', 'N']) for ids in self._all_ids_atomwise])
        result = {}
        frame_count = 0
        frames = self.nb_frames

        for ts in self._universe.trajectory[self._trajectory_slice]:

            water_coordinates = self._water.positions
            selection_coordinates = self._da_selection.positions

            selection_tree = _sp.cKDTree(selection_coordinates)
            if self._nb_acceptors > 0 and self._nb_donors > 0:
                d_tree = _sp.cKDTree(self._donors.positions)
                a_tree = _sp.cKDTree(self._acceptors.positions)
            hydrogen_coordinates = self._hydrogen.positions

            water_tree = _sp.cKDTree(water_coordinates, leafsize=32)
            local_water_index = []
            [local_water_index.extend(l) for l in water_tree.query_ball_point(selection_coordinates, float(around_radius))]
            local_water_index = _np.unique(local_water_index)

            local_water_coordinates = water_coordinates[local_water_index]
            local_water_tree = _sp.cKDTree(local_water_coordinates)

            local_pairs = [(i, local_water_index[j]+self._first_water_id) for i, bla in enumerate(selection_tree.query_ball_tree(local_water_tree, self.distance)) for j in bla]
            water_pairs = [(local_water_index[p[0]]+self._first_water_id, local_water_index[p[1]]+self._first_water_id) for p in local_water_tree.query_pairs(self.distance)]
            if self._nb_acceptors > 0 and self._nb_donors > 0:
                da_pairs = _np.array([[i, j] for i,donors in enumerate(a_tree.query_ball_tree(d_tree, self.distance)) for j in donors])
                da_pairs[:,0] += self._nb_donors
            else: da_pairs = []
            if exclude_backbone_backbone and da_pairs.size != 0: da_pairs = da_pairs[_np.logical_not(_np.all(backbone_filter[da_pairs], axis=1))]

            if self.check_angle:
                all_coordinates = _np.vstack((selection_coordinates, water_coordinates))
                hbonds = _hf.check_angle(list(da_pairs)+water_pairs+local_pairs, self.heavy2hydrogen, all_coordinates, hydrogen_coordinates, self.cut_angle)
            else:
                hbonds = list(da_pairs) + water_pairs + local_pairs

            hbonds = _np.array(hbonds)

            sorted_bonds = _np.sort(hbonds)
            check = self._resids[sorted_bonds]
            check = check[:,0] < check[:,1]
            if self.residuewise: frame_res = [self._all_ids[i] + ':' + self._all_ids[j] if check[ii] else self._all_ids[j] + ':' + self._all_ids[i] for ii, (i, j) in enumerate(sorted_bonds)]
            else: frame_res = [self._all_ids_atomwise[i] + ':' + self._all_ids_atomwise[j] if check[ii] else self._all_ids_atomwise[j] + ':' + self._all_ids_atomwise[i] for ii, (i, j) in enumerate(sorted_bonds)]
            for bond in frame_res:
                if not self.check_angle:
                    a, b = bond.split(':')
                    if a.split('-')[:3] == b.split('-')[:3]: continue
                try:
                    result[bond][frame_count] = True
                except:
                    result[bond] = _np.zeros(frames, dtype=bool)
                    result[bond][frame_count] = True
            frame_count+=1
        self._set_results(result)
        if self._ions_list: self.set_add_ion_interactions(True)

    def filter_only_backbone_bonds(self, use_filtered=True):
        assert self.residuewise == False, 'residuewise has to be False for this filter to work.'
        if use_filtered: results = self.filtered_results
        else: results = self.initial_results
        if len(results) == 0: raise AssertionError('nothing to filter!')
        atomname_pairs = _np.array([[bond.split(':')[0].split('-')[3],bond.split(':')[1].split('-')[3]]for bond in results])
        backbone_backbone_filter = _np.array([(name[0] in ['O','N'] and name[1] in ['O','N']) for name in atomname_pairs])
        self.filtered_results = {key:results[key] for i,key in enumerate(results) if backbone_backbone_filter[i]}
        self._generate_filtered_graph_from_filtered_results()
        self.applied_filters['backbone'] = True

    def filter_water_in_hydration_shells(self, hydration_shells=3, use_filtered=True):
        if use_filtered:
            results = self.filtered_results
            graph = self.filtered_graph
        else:
            results = self.initial_results
            graph = self.initial_graph
        if len(results) == 0: raise AssertionError('nothing to filter!')
        filtered_result = {}
        hysh_lengths = []
        res_nodes = [node for node in graph.nodes() if node.split('-')[1] not in ['TIP3', 'HOH']]
        for res_node in res_nodes:
            hysh_lengths += _nx.single_source_shortest_path_length(graph, res_node, hydration_shells)

        graph.remove_nodes_from(set(graph.nodes())-set(hysh_lengths))

        filtered_result = {}
        for edge in graph.edges():
            bond_name = edge[0]+':'+edge[1]
            check = edge[1]+':'+edge[0]
            try : filtered_result[bond_name] = results[bond_name]
            except KeyError:
                filtered_result[check] = results[check]
        if self.applied_filters['shells'] != None: self.applied_filters['shells'] = _np.min(self.applied_filters['shells'], hydration_shells)
        else: self.applied_filters['shells'] = hydration_shells
        self.filtered_results = filtered_result
        self._generate_filtered_graph_from_filtered_results()

    def compute_i4_bond_lengths(self, filter_above_factor = None, average_over_segments=True, print_table=True, use_filtered=True):
        if self.residuewise: raise AssertionError('need to initialize with residuewise=False!')
        if use_filtered: results = self.filtered_results
        else: results = self.initial_results
        if len(results) == 0: raise AssertionError('nothing to filter!')
        ids = _np.array(self._all_ids_atomwise)
        atomname_pairs = _np.array([[bond.split(':')[0].split('-')[3],bond.split(':')[1].split('-')[3]]for bond in results])
        resid_pairs = _np.array([[bond.split(':')[0].split('-')[2],bond.split(':')[1].split('-')[2]]for bond in results], dtype=_np.int)
        backbone_backbone_filter = _np.array([(name[0] in ['O','N'] and name[1] in ['O','N']) for name in atomname_pairs])
        i4_filter = (_np.in1d(resid_pairs[:, 0] - resid_pairs[:,1], [3,4,5]) | _np.in1d(resid_pairs[:, 0] - resid_pairs[:,1], [-3,-4,-5])) & backbone_backbone_filter
        values = _np.array([results[bond] for bond in results])[i4_filter].T
        keypairs = _np.array([[bond.split(':')[0], bond.split(':')[1]] for bond in results])[i4_filter]
        keyindex = _np.array([[_np.where(ids == keypair[0])[0][0], _np.where(ids == keypair[1])[0][0]] for keypair in keypairs])
        l = _np.array([bond for bond in results])[i4_filter]
        assert keypairs.size != 0, 'No i4-bonds found'
        dists = _np.zeros(values.shape)
        for i,ts in enumerate(self._universe.trajectory[self._trajectory_slice]):
            water_coordinates = self._water.positions
            select_coordinates = self._da_selection.positions
            all_coordinates = _np.vstack((select_coordinates,water_coordinates))
            coords = all_coordinates[keyindex]
            d = coords[:,0]-coords[:,1]
            dists[i][values[i]] = _np.linalg.norm(d, axis=1)[values[i]]
        dists = dists.T
        dists = _np.array([dists[i][dists[i] != 0].mean() for i in range(len(dists))])
        resid_sort_index = _np.argsort(resid_pairs[:,0][i4_filter])
        l = l[resid_sort_index]
        dists = dists[resid_sort_index]
        if average_over_segments:
            resids_average=resid_pairs[i4_filter][resid_sort_index]
            resnames_average = _np.array([[bond.split(':')[0].split('-')[1],bond.split(':')[1].split('-')[1]]for bond in results])[i4_filter][resid_sort_index]
            atomname_average = atomname_pairs[i4_filter][resid_sort_index]
            l_average = []
            dists_average = []
            already_checked = []
            for i in range(len(l)):
                print(resids_average[i])
                if list(resids_average[i]) not in already_checked:
                    l_average.append(resnames_average[i][0]+'-'+str(resids_average[i][0])+'-'+atomname_average[i][0]+':'+resnames_average[i][1]+'-'+str(resids_average[i][1])+'-'+atomname_average[i][1])
                    dists_average.append(dists[(resids_average==resids_average[i]).all(1)].mean())
                    already_checked.append(list(resids_average[i]))
            l = _np.array(l_average)
            dists = _np.array(dists_average)
        if filter_above_factor != None:
            mean, std = dists.mean(), dists.std()
            above_index = dists >= mean + filter_above_factor * std
            l = l[above_index]
            dists = dists[above_index]
        if print_table:
            t = ''
            for i in range(len(l)):
                t+=l[i]+'.'*(45-len(l[i]))+str(_np.round(dists[i],2))+'\n'
            print(t)
        return {l[i]:dists[i] for i in range(len(l))}


    def compute_i4_bonds(self, use_filtered=True, return_print_table=False, print_table=False):
        if self.residuewise: raise AssertionError('need to initialize with residuewise=False!')
        if use_filtered: results = self.filtered_results
        else: results = self.initial_results
        if len(results) == 0: raise AssertionError('nothing to filter!')
        l = _np.array([bond for bond in results])
        segname_pairs = _np.array([[bond.split(':')[0].split('-')[0],bond.split(':')[1].split('-')[0]]for bond in results])
        atomname_pairs = _np.array([[bond.split(':')[0].split('-')[3],bond.split(':')[1].split('-')[3]]for bond in results])
        resnames_pairs = _np.array([[bond.split(':')[0].split('-')[1],bond.split(':')[1].split('-')[1]]for bond in results])
        resid_pairs = _np.array([[bond.split(':')[0].split('-')[2],bond.split(':')[1].split('-')[2]]for bond in results], dtype=_np.int)
        ser_thr_filter = _np.array([(ids[0] in ['SER', 'THR']) or (ids[1] in ['SER', 'THR']) for i, ids in enumerate(resnames_pairs)])
        glu_asp_filter = _np.array([(ids[0] in ['GLU', 'ASP']) or (ids[1] in ['GLU', 'ASP']) for i, ids in enumerate(resnames_pairs)])
        backbone_filter = _np.array([not (name[0] in ['O','N'] or name[1] in ['O','N']) for name in atomname_pairs])
        backbone_backbone_filter = _np.array([not (name[0] in ['O','N'] and name[1] in ['O','N']) for name in atomname_pairs])
        is_backbone = _np.array([[pair[0] in ['O', 'N'], pair[1] in ['O', 'N']] for pair in atomname_pairs], dtype=_np.bool)
        both_filter = ser_thr_filter & glu_asp_filter & backbone_filter
        res_thr_filter = _np.array([[pair[0] in ['SER', 'THR'], pair[1] in ['SER', 'THR']] for pair in resnames_pairs], dtype=_np.bool)
        i4_filter = ((_np.in1d(resid_pairs[:, 0] - resid_pairs[:,1], [3,4,5]) & res_thr_filter[:,0] & is_backbone[:,1]) | \
                    (_np.in1d(resid_pairs[:, 0] - resid_pairs[:,1], [-3,-4,-5]) & res_thr_filter[:,1] & is_backbone[:,0])) & backbone_backbone_filter

        #if i4_filter.sum() == 0: raise AssertionError('No i-4 intrahelical hbonds found for SER/THR and no sidechain hbonds between SER/THR and GLU/ASP found')
        class_dict = {}
        res3={}
        output = ''
        translate_key = {l[i]:segname_pairs[i][0]+'-'+resnames_pairs[i][0]+'-'+str(resid_pairs[i][0]+self.add_missing_residues)+'-'+atomname_pairs[i][0]+':'+segname_pairs[i][1]+'-'+resnames_pairs[i][1]+'-'+str(resid_pairs[i][1]+self.add_missing_residues)+'-'+atomname_pairs[i][1] for i in range(len(results))}
        for i in range(len(l)):
            key = l[i]
            if both_filter[i]:
                resid = resid_pairs[i][res_thr_filter[i]]
                segname = segname_pairs[i][res_thr_filter[i]]
                resname = resnames_pairs[i][res_thr_filter[i]]
                atomname = atomname_pairs[i][res_thr_filter[i]]
                other_bond_filter = ((resid_pairs == resid) & (segname_pairs == segname) & (resnames_pairs == resname) & (atomname_pairs == atomname)).any(1) & i4_filter
                for j in _np.nonzero(other_bond_filter)[0]:
                    sel_i, sel_j = _np.logical_not(res_thr_filter[i]), _np.logical_not(res_thr_filter[j])
                    if not (segname_pairs[i][sel_i] == segname_pairs[j][sel_j] and resnames_pairs[j][sel_j] == resnames_pairs[i][sel_i] and resid_pairs[j][sel_j] == resid_pairs[i][sel_i]):
                        joint_timeseries = results[l[j]] & results[l[i]]
                        if joint_timeseries.mean() != 0.0:
                            res3[l[i]] = l[j]
                class_dict[l[i]] = 2
            if i4_filter[i]:
                class_dict[l[i]] = 1
        segnames = _np.unique(segname_pairs)
        segname_results = {segname:{1:{}, 2:{}} for segname in segnames}
        temp = []
        for segnamea, segnameb in combinations(segname_results, 2): temp.append(segnamea+'-'+segnameb)
        for comb in temp: segname_results[comb]={1:{}, 2:{}}
        for key in class_dict:
            residue1, residue2 = key.split(':')
            sn1, sn2 = residue1.split('-')[0], residue2.split('-')[0]
            comb, comb_check = '-'.join((sn1,sn2)), '-'.join((sn2,sn1))
            try: segname_results[comb][class_dict[key]][key] = results[key]
            except:
                try: segname_results[comb_check][class_dict[key]][key] = results[key]
                except:
                    try: segname_results[sn1][class_dict[key]][key] = results[key]
                    except:
                        try: segname_results[sn2][class_dict[key]][key] = results[key]
                        except: pass
        for segname in segname_results:
            class3_text=''
            for class_id in segname_results[segname]:
                for key in segname_results[segname][class_id]:
                    output += translate_key[key]+ '.'*(45-len(translate_key[key])) + ' | '+ 'class ' + str(class_dict[key]) + ' | '+ str(_np.round(results[key].mean()*100,1))+'\n'
                    if key in res3:
                       class3_text += translate_key[key]+ '.'*(45-len(translate_key[key])) + ' | '+ 'class ' + str(3) + ' | '+ str(_np.round((results[key] & results[res3[key]]).mean()*100,1))+'\n'
                       class3_text += translate_key[res3[key]]+ '.'*(45-len(translate_key[res3[key]])) + ' | '+ 'class ' + str(3) + ' | '+ str(_np.round((results[key] & results[res3[key]]).mean()*100,1))+'\n'
            output+=class3_text+'\n'
        if print_table:
            print(output)
        if return_print_table:
            return output
        return class_dict, res3

    def compute_i4_motif_distribution(self):
        class_dict, res3 = self.compute_i4_bonds(print_table=False)
        z = {i:[] for i in range(1,4)}
        for key in class_dict:
            segida, resna, resida, atoma = key.split(':')[0].split('-')
            segidb, resnb, residb, atomb = key.split(':')[1].split('-')
            atoms = self._mda_selection.select_atoms("(segid {} and resid {} and name {}) or (segid {} and resid {} and name {})".format(segida, resida, atoma, segidb, residb, atomb))
            z[class_dict[key]]+=[atom for atom in atoms]
        for key in res3:
            for i in range(2):
                segida, resna, resida, atoma = key.split(':')[0].split('-')
                segidb, resnb, residb, atomb = key.split(':')[1].split('-')
                atoms = self._mda_selection.select_atoms("(segid {} and resid {} and name {}) or (segid {} and resid {} and name {})".format(segida, resida, atoma, segidb, residb, atomb))
                z[3]+=[atom for atom in atoms]
                if i==0: key = res3[key]
        for key in z:
            if isinstance(z[key], list) and z[key] != []:
                z[key] = _MDAnalysis.core.AtomGroup.AtomGroup(z[key])
        z_plot = {i:[] for i in range(1,4)}
        for ts in self._universe.trajectory[self._trajectory_slice]:
            for key in z:
                if not isinstance(z[key], list):
                    if key==3: z_plot[key].append(z[key].positions[:,2].reshape(-1,4).mean(1))
                    else: z_plot[key].append(z[key].positions[:,2].reshape(-1,2).mean(1))
        z_plot = [_np.array(z_plot[key]).mean(0) if z_plot[key] != [] else _np.array(z_plot[key]) for key in z_plot]
        self._i4_distribution = z_plot
        return z_plot

    def draw_i4_motif_distribution(self, filename=None, return_figure=False):
        class_dict, res3 = self.compute_i4_bonds(print_table=False)
        z = {i:[] for i in range(1,4)}

        for key in class_dict:
            segida, resna, resida, atoma = key.split(':')[0].split('-')
            segidb, resnb, residb, atomb = key.split(':')[1].split('-')
            atoms = self._mda_selection.select_atoms("(segid {} and resid {} and name {}) or (segid {} and resid {} and name {})".format(segida, resida, atoma, segidb, residb, atomb))
            z[class_dict[key]]+=[atom for atom in atoms]
        for key in res3:
            for i in range(2):
                segida, resna, resida, atoma = key.split(':')[0].split('-')
                segidb, resnb, residb, atomb = key.split(':')[1].split('-')
                atoms = self._mda_selection.select_atoms("(segid {} and resid {} and name {}) or (segid {} and resid {} and name {})".format(segida, resida, atoma, segidb, residb, atomb))
                z[3]+=[atom for atom in atoms]
                if i==0: key = res3[key]
        for key in z:
            if isinstance(z[key], list) and z[key] != []:
                z[key] = _MDAnalysis.core.AtomGroup.AtomGroup(z[key])
        z_plot = {i:[] for i in range(1,4)}
        for ts in self._universe.trajectory[self._trajectory_slice]:
            for key in z:
                if not isinstance(z[key], list):
                    if key==3: z_plot[key].append(z[key].positions[:,2].reshape(-1,4).mean(1))
                    else: z_plot[key].append(z[key].positions[:,2].reshape(-1,2).mean(1))
        z_plot = [_np.array(z_plot[key]).mean(0)if z_plot[key] != [] else _np.array(z_plot[key]) for key in z_plot]
        #print(z_plot)
        mi, ma = _np.min([a.min() for a in z_plot if a.size>0]), _np.max([a.max() for a in z_plot if a.size>0])
        fig, ax = _plt.subplots()
        _plt.hist(z_plot,
                 _np.linspace(mi, ma, 10),
                 histtype='bar',
                 orientation=u'horizontal',
                 stacked=False,
                 fill=True,
                 label=['intra-helical','inter-helical','intra-helical + inter-helical'],
                 alpha=0.8, # opacity of the bars
                 edgecolor = "k")
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        _plt.yticks(_np.round(_np.linspace(mi, ma, 10), 0))
        _plt.ylim([mi-1,ma+1])
        _plt.xlabel('Count' , fontsize = 16)
        _plt.ylabel('Z Coordinate' , fontsize = 16)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        _plt.legend()
        return_str = self.compute_i4_bonds(return_print_table=True) + '\n'
        return_str += 'Z-Distribution:\n'
        binedges = _np.linspace(mi, ma, 10)
        interv = _np.round(_np.linspace(mi, ma, 10), 0).astype(_np.str)
        return_str2 = ''
        return_str2 += 'z_intervals ' + ' '.join([interv[i]+'_'+interv[i+1] for i in range(len(interv)-1)]) + '\n'
        return_str2 += 'intra-helical ' + ' '.join(_np.histogram(z_plot[0], binedges)[0].astype(_np.str)) + '\n'
        return_str2 += 'inter-helical ' + ' '.join(_np.histogram(z_plot[1], binedges)[0].astype(_np.str)) + '\n'
        return_str2 += 'intra+inter-helical ' + ' '.join(_np.histogram(z_plot[2], binedges)[0].astype(_np.str)) + '\n'
        return_str += _hf.string_in_columns(return_str2)
        if return_figure:
            _plt.close()
            return fig, return_str
        self._save_or_draw(filename, return_figure=return_figure)

    def draw_hbond_layer_occupancy_histogram(self, min_occupancy, occupancy_step, use_filtered=True, filename=None, return_figure=False):
        if use_filtered: results = self.filtered_results
        else: results = self.initial_results
        if len(results) == 0: raise AssertionError('nothing to filter!')
        occupancies = _np.arange(min_occupancy, 1., occupancy_step)
        bonds_list = [_hf.filter_occupancy(results, occupancy) for occupancy in occupancies]
        graphs = [_hf.dict2graph(bonds, self.residuewise) for bonds in bonds_list]
        hysh_counts = _np.zeros((occupancies.size, 3))

        for i, occupancy in enumerate(occupancies):
            hysh_shells = [[],[],[]]
            res_nodes = [node for node in graphs[i].nodes() if node.split('-')[1] not in ['TIP3', 'HOH']]
            for res_node in res_nodes:
                hysh_lengths = _nx.single_source_shortest_path_length(graphs[i],res_node, 3)
                del hysh_lengths[res_node]
                for hydration_water in hysh_lengths:
                    if hydration_water.split('-')[1]!='TIP3': continue
                    hysh_shells[hysh_lengths[hydration_water]-1].append(hash(hydration_water))

            hysh_counts[i][2] = _np.unique(_np.array(hysh_shells[2])[_np.in1d(hysh_shells[2],hysh_shells[1]+hysh_shells[0], invert=True)]).size
            hysh_counts[i][1] = _np.unique(_np.array(hysh_shells[1])[_np.in1d(hysh_shells[1],hysh_shells[0], invert=True)]).size
            hysh_counts[i][0] = _np.unique(hysh_shells[0]).size

        labels = ['{0:.{1}f}'.format(occupancy, 2) for occupancy in occupancies]
        diff_counts = hysh_counts.T
        fig, ax = _plt.subplots()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        second = _plt.bar(_np.arange(len(diff_counts[1]), dtype=int)*4, diff_counts[1], label='Second Layer')
        _plt.xticks(_np.arange(len(diff_counts[1]), dtype=int)*4, labels)
        first = _plt.bar(_np.arange(len(diff_counts[0]), dtype=int)*4 - 1, diff_counts[0], label='First Layer')
        third = _plt.bar(_np.arange(len(diff_counts[2]), dtype=int)*4 + 1, diff_counts[2], label='Third Layer')
        _plt.legend(handles=[first, second, third])
        _plt.xlabel('Occupancy rate' , fontsize = 16)
        _plt.ylabel('Count' , fontsize = 16)
        if return_figure:
            _plt.close()
            return fig, 'NoData'
        self._save_or_draw(filename, return_figure=return_figure)
