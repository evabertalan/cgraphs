from . import helperfunctions as _hf
import shutil
import numpy as np
import copy
# import pdb
import MDAnalysis as _mda
from .proteingraphanalyser import ProteinGraphAnalyser
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm


class CompareTwo(ProteinGraphAnalyser):
    def __init__(self, type_option, pdb1=None, pdb2=None, psf1=None, psf2=None, dcd1=None, dcd2=None, target_folder='', name1=None, name2=None, plot_parameters={}):
        self.type_option = type_option
        self.target_folder = target_folder
        if self.type_option == 'pdb' and pdb1 and pdb2:
            self.name1 = _hf.retrieve_pdb_code(pdb1, '.pdb')
            self.name2 = _hf.retrieve_pdb_code(pdb2, '.pdb')
            self.compare_folder = _hf.create_directory(target_folder+'/workfolder/compare_'+self.name1+'_'+self.name2)+'/'
            shutil.copy(pdb1, self.compare_folder+self.name1+'.pdb')
            shutil.copy(pdb2, self.compare_folder+self.name2+'.pdb')

            ProteinGraphAnalyser.__init__(self, pdb_root_folder=self.compare_folder, target_folder=target_folder, reference_pdb=pdb1, plot_parameters=plot_parameters)
            self.logger.info('COMPARE STRUCTURES '+ self.name1 + ' WITH ' + self.name2)
            ProteinGraphAnalyser.align_structures(self, superimposition_threshold=30)

        elif self.type_option == 'dcd' and psf1 and psf2 and dcd1 and dcd2:
            self.name1, self.name2 = name1, name2
            self.compare_folder = _hf.create_directory(target_folder+'/workfolder/compare_'+self.name1+'_'+self.name2)+'/'

            ProteinGraphAnalyser.__init__(self, pdb_root_folder=self.compare_folder, target_folder=target_folder, type_option='dcd', psf_files=[psf1, psf2], dcd_files=[dcd1, dcd2], sim_names=[name1, name2], plot_parameters=plot_parameters)
            self.logger.info('COMPARE SIMULATIONS '+ self.name1 + ' WITH ' + self.name2)
        else: self.logger.warning('Required files are missing for the calculation')

    def construct_comparison_objects(self, occupancy=None):
        self.occupancy = occupancy
        if self.type_option == 'dcd':
            if self.occupancy:
                self.logger.info('Filter graphs for '+self.name1 + ' and ' + self.name2+' on occupancy '+str(self.occupancy*100)+'%')
                graphs = []
                for i, (name, objects) in enumerate(self.graph_coord_objects.items()):
                    u = _mda.Universe(objects['psf'], objects['dcd'])
                    self.selection = objects['selection'] if 'selection' in objects.keys() else selection
                    sel = u.select_atoms(self.selection)
                    mda = sel.select_atoms(self.selection)
                    wba = copy.deepcopy(objects['wba'])
                    wba.filter_occupancy(self.occupancy)
                    g = wba.filtered_graph
                    graphs.append(g)
                    self.graph_coord_objects[name].update({'mda': mda})
                    self.reference_coordinates={}
                    self.add_reference_from_structure(mda, g)

                self.graph_coord_objects[self.name1]['graph'] = graphs[0]
                self.graph_coord_objects[self.name2]['graph'] = graphs[1]
            else: self.logger.info('occupancy has to be specified for trajectory comparison!')
            self.pca_positions = _hf.calculate_pca_positions(self.reference_coordinates)


    def plot_graph_comparison(self, color1='blue', color2='green', label_nodes=True, label_edges=True, xlabel='PCA projected xy plane', ylabel='Z coordinates ($\AA$)', color_propka=False, color_data=False, node_color_selection='protein', node_color_map='bwr', calcualte_distance=False):

        if len(self.graph_coord_objects.items()) != 2: self.logger.warning('There are '+str(len(self.graph_coord_objects.items()))+' structures selected. Graph comparison is possible for exactly two structures.')
        else:
            self.logger.info('Plot comparison '+self.graph_type+' graph for '+ self.name1 + ' with ' + self.name2+str(' colored by pKa' if color_propka else '')+str(' colored by user defined values' if color_data else '')+str(' with labels' if label_nodes else ''))
            if 'graph' in self.graph_coord_objects[self.name1].keys() and 'graph' in self.graph_coord_objects[self.name2].keys():
                if hasattr(self, 'occupancy'):
                    occupancy = self.occupancy
                    wba1 = copy.deepcopy(self.graph_coord_objects[self.name1]['wba'])
                    wba1.filter_occupancy(occupancy)
                    graph1 = wba1.filtered_graph

                    wba2 = copy.deepcopy(self.graph_coord_objects[self.name2]['wba']) #TEST THIS
                    wba2.filter_occupancy(occupancy)
                    graph2 = wba2.filtered_graph
                else:
                    occupancy = None
                    graph1 = self.graph_coord_objects[self.name1]['graph']
                    graph2 = self.graph_coord_objects[self.name2]['graph']

                pos1 = self._get_node_positions(self.graph_coord_objects[self.name1], pca=False)
                pos2 = self._get_node_positions(self.graph_coord_objects[self.name2], pca=False)

                if calcualte_distance:

                    bond_distances1 = self._get_edge_distance(self.graph_coord_objects[self.name1], self.name1, self.compare_folder, denumber_waters=True)
                    bond_distances2 = self._get_edge_distance(self.graph_coord_objects[self.name2], self.name2, self.compare_folder, denumber_waters=True)

                conserved_edges=[]
                conserved_nodes=[]
                conserved_waters = {}

                all_pos = {}
                for i, pos in enumerate([pos1, pos2]):
                    for key, value in pos.items():
                        if key.split('-')[1] in _hf.water_types:
                            key=f"{i+1}-{key.split('-')[1]}-{key.split('-')[2]}"
                            if i == 0:
                                conserved_waters.update({key: {
                                    'conserved': False,
                                    'pair': None,
                                    'pos': value
                                    }})
                            if i == 1: #get waters from pos2 since condensed water are initialized as pos1
                                for conserved_key, conseved_vals in conserved_waters.items():
                                    if ((value[0]-conseved_vals['pos'][0])**2 + (value[1]-conseved_vals['pos'][1])**2 + (value[2]-conseved_vals['pos'][2])**2 <= float(1)**2):
                                        # breakpoint()
                                        conserved_waters[conserved_key]['conserved'] = True
                                        conserved_waters[conserved_key]['pair'] = key
                        if key not in all_pos.keys(): all_pos.update({ key:value })

                conserved_waters = {k: v for k,v in conserved_waters.items() if v['conserved'] == True}
                node_pca_pos = _hf.calculate_pca_positions(all_pos)


                plot_name = 'H-bond' if self.graph_type == 'hbond' else 'Water wire'
                fig, ax = _hf.create_plot(title=f'{plot_name} graph comparison of {self.name1} with {self.name2} \n Selection: {self.selection[1:-16]}', xlabel=xlabel, ylabel=ylabel, plot_parameters=self.plot_parameters)

                fig_cons, ax_cons = _hf.create_plot(title=f'{plot_name} conserved graph of {self.name1} and {self.name2} \n Selection: {self.selection[1:-16]}', xlabel=xlabel, ylabel=ylabel, plot_parameters=self.plot_parameters)

                if calcualte_distance:
                    fig_dist, ax_dist = _hf.create_plot(title=f'{plot_name} conserved graph of {self.name1} and {self.name2} \n Selection: {self.selection[1:-16]}', xlabel=xlabel, ylabel=ylabel, plot_parameters=self.plot_parameters)

                dist_plot_data = []
                dist_info = {}
                for e in graph1.edges:
                    e0 = e[0] if e[0].split('-')[1] not in _hf.water_types else f"1-{e[0].split('-')[1]}-{e[0].split('-')[2]}"
                    e1 = e[1] if e[1].split('-')[1] not in _hf.water_types else f"1-{e[1].split('-')[1]}-{e[1].split('-')[2]}"
                    if _hf.is_conserved_edge(np.array([[e2[0], e2[1]] for e2 in graph2.edges]), e0, e1) or e0 in conserved_waters.keys() or e1 in conserved_waters.keys():
                        color = self.plot_parameters['graph_color']
                        conserved_edges.append(e)
                    else: color = color1

                    if e0 in node_pca_pos.keys() and e1 in node_pca_pos.keys():
                        edge_line = [node_pca_pos[e0], node_pca_pos[e1]]
                        x=[edge_line[0][0], edge_line[1][0]]
                        y=[edge_line[0][1], edge_line[1][1]]
                        ax.plot(x, y, color=color, marker='o', linewidth=self.plot_parameters['edge_width'], markersize=self.plot_parameters['node_size']*0.01, markerfacecolor=color, markeredgecolor=color)
                        if e in conserved_edges:
                            ax_cons.plot(x, y, color= self.plot_parameters['graph_color'], marker='o', linewidth=self.plot_parameters['edge_width'], markersize=self.plot_parameters['node_size']*0.01, markerfacecolor= self.plot_parameters['graph_color'], markeredgecolor= self.plot_parameters['graph_color'])

                            if calcualte_distance:
                                e0_chain_id, e0_res_name, e0_res_id  = _hf.get_node_name_pats(e[0])
                                e_0 = f"{e0_chain_id}-{e0_res_name}-{e0_res_id}"

                                e1_chain_id, e1_res_name, e1_res_id  = _hf.get_node_name_pats(e[1])
                                e_1 = f"{e1_chain_id}-{e1_res_name}-{e1_res_id}"
                                if e[0].split('-')[1] in _hf.water_types:
                                    e_0 = 'HOH'
                                if e[1].split('-')[1] in _hf.water_types:
                                    e_1 = 'HOH'

                                key = f"{e_0}_{e_1}"
                                if key in bond_distances1:
                                    dist_1 = bond_distances1[f"{e_0}_{e_1}"]
                                if key in bond_distances2:
                                    dist_2 = bond_distances2[f"{e_0}_{e_1}"]
                                key = f"{e_1}_{e_0}"
                                if key in bond_distances1:
                                    dist_1 = bond_distances1[f"{e_1}_{e_0}"]
                                if key in bond_distances2:
                                    dist_2 = bond_distances2[f"{e_1}_{e_0}"]
                                dist = dist_1 - dist_2

                                dist_plot_data.append([x,y,dist])
                                dist_info.update({key: dist})

                                # if label_edges:
                                #     ax.annotate(round(dist, 1), (x[0] + (x[1]-x[0])/2, y[0] + (y[1]-1.0-y[0])/2), color='blue',  fontsize=self.plot_parameters['edge_label_size'])

                if calcualte_distance:
                    ax_dist, fig_dist = self._plot_dist_plot(dist_plot_data, ax_dist, fig_dist, node_color_map, label_edges)
                    if label_nodes:
                        is_backbone = '_backbone' if self.include_backbone_sidechain else ''
                        _hf.write_text_file(f'{self.compare_folder}H-bond_graph_distances_changes{is_backbone}.txt',[f"{edge} {round(distance_change,3)}\n" for edge, distance_change in dist_info.items()])

                for e in graph2.edges:
                    e0 = e[0] if e[0].split('-')[1] not in _hf.water_types else '2-'+e[0].split('-')[1]+'-'+e[0].split('-')[2]
                    e1 = e[1] if e[1].split('-')[1] not in _hf.water_types else '2-'+e[1].split('-')[1]+'-'+e[1].split('-')[2]
                    conserved_pair = [v['pair'] for v in conserved_waters.values()]
                    if e0 in node_pca_pos.keys() and e1 in node_pca_pos.keys():
                        if not _hf.is_conserved_edge(np.array([[e2[0], e2[1]] for e2 in graph1.edges]), e0, e1) and e0 not in conserved_pair and e1 not in conserved_pair:
                            edge_line = [node_pca_pos[e0], node_pca_pos[e1]]
                            x=[edge_line[0][0], edge_line[1][0]]
                            y=[edge_line[0][1], edge_line[1][1]]
                            ax.plot(x, y, color=color2, marker='o', linewidth=self.plot_parameters['edge_width'], markersize=self.plot_parameters['node_size']*0.01, markerfacecolor=color2, markeredgecolor=color2)

                color_info = {}
                if color_propka or color_data:
                    color_info1 = self._custom_node_coloring(self.graph_coord_objects[self.name1],self.name1, color_propka, color_data, node_color_selection)
                    color_info2 = self._custom_node_coloring(self.graph_coord_objects[self.name2], self.name2, color_propka, color_data, node_color_selection)
                    if len(color_info1) and len(color_info2):
                        for key1, value1 in color_info1.items():
                            if key1 in color_info2.keys():
                                value_diff = float(color_info1[key1]) - float(color_info2[key1])

                                color_info.update({key1: value_diff})
                        value_colors, cmap, norm = _hf.get_color_map(color_info, color_map=node_color_map, center=True)
                    color_bar_label = 'Amino acid data value' if color_data else 'pKa value'
                    if label_nodes:
                        is_backbone = '_backbone' if self.include_backbone_sidechain else ''
                        _hf.write_text_file(f'{self.compare_folder}H-bond_graph_pKa_changes{is_backbone}.txt',[f"{edge} {round(distance_change,3)}\n" for edge, distance_change in color_info.items()])

                for n in graph1.nodes:
                    n = n if n.split('-')[1] not in _hf.water_types else '1-'+n.split('-')[1]+'-'+n.split('-')[2]
                    if n in node_pca_pos.keys():
                        if n in graph2.nodes:
                            conserved_nodes.append(n)
                            color = self.plot_parameters['graph_color']
                        elif n in conserved_waters.keys():
                            color = self.plot_parameters['water_node_color']
                        else: color = color1
                        if n.split('-')[1] in _hf.water_types and n not in conserved_waters.keys():
                            ax.scatter(node_pca_pos[n][0], node_pca_pos[n][1],color=self.plot_parameters['water_node_color'], s=self.plot_parameters['node_size']*0.8, zorder=5, edgecolors=color)
                        else: ax.scatter(node_pca_pos[n][0], node_pca_pos[n][1], s=self.plot_parameters['node_size'], color=color, zorder=5)
                        if n in conserved_nodes or n in conserved_waters.keys():
                            if (color_propka or color_data) and n in color_info.keys():
                                if round(color_info[n],2) == 0: color = self.plot_parameters['graph_color']
                                else: color = value_colors[n]
                            if n.split('-')[1] in _hf.water_types:
                                ax_cons.scatter(node_pca_pos[n][0], node_pca_pos[n][1],color=self.plot_parameters['water_node_color'], s=self.plot_parameters['node_size']*0.8, zorder=5, edgecolors=color)
                                if calcualte_distance:
                                    ax_dist.scatter(node_pca_pos[n][0], node_pca_pos[n][1],color=self.plot_parameters['water_node_color'], s=self.plot_parameters['node_size']*0.8, zorder=5, edgecolors=color)
                            else:
                                ax_cons.scatter(node_pca_pos[n][0], node_pca_pos[n][1], s=self.plot_parameters['node_size'], color=color, zorder=5, edgecolors=self.plot_parameters['graph_color'])
                                if calcualte_distance:
                                    ax_dist.scatter(node_pca_pos[n][0], node_pca_pos[n][1], s=self.plot_parameters['node_size'], color=self.plot_parameters['graph_color'], zorder=5, edgecolors=self.plot_parameters['graph_color'])


                for n in graph2.nodes:
                    n = n if n.split('-')[1] not in _hf.water_types else '2-'+n.split('-')[1]+'-'+n.split('-')[2]
                    if n in node_pca_pos.keys():
                        if n not in graph1.nodes and n not in conserved_pair:
                            if n.split('-')[1] in _hf.water_types:
                                ax.scatter(node_pca_pos[n][0], node_pca_pos[n][1], color=self.plot_parameters['water_node_color'], s=self.plot_parameters['node_size']*0.8, zorder=5, edgecolors=color2)
                            else: ax.scatter(node_pca_pos[n][0], node_pca_pos[n][1], s=self.plot_parameters['node_size'], color=color2, zorder=5)

                if label_nodes:
                    for n in all_pos.keys():
                        values = node_pca_pos[n]
                        chain_id, res_name, res_id = _hf.get_node_name_pats(n)
                        if res_name in _hf.water_types:
                            if n not in conserved_waters.keys() and n not in conserved_pair:
                                ax.annotate(f'W{res_id}', (values[0]+0.2, values[1]-0.25), fontsize=self.plot_parameters['node_label_size'])
                            # if n in conserved_waters.keys():
                                # ax_cons.annotate('W', (values[0]+0.2, values[1]-0.25), fontsize=self.plot_parameters['node_label_size'])
                        elif res_name in _hf.amino_d.keys():
                            l = f'{chain_id}-{_hf.amino_d[res_name]}{res_id}' if self.plot_parameters['show_chain_label'] else f'{_hf.amino_d[res_name]}{res_id}'
                            ax.annotate(l, (values[0]+0.2, values[1]-0.25), fontsize=self.plot_parameters['node_label_size'])
                            if n in conserved_nodes:
                                l = f'{chain_id}-{_hf.amino_d[res_name]}{res_id}' if self.plot_parameters['show_chain_label'] else f'{_hf.amino_d[res_name]}{res_id}'
                                ax_cons.annotate(l, (values[0]+0.2, values[1]-0.25), fontsize=self.plot_parameters['node_label_size'])
                                if calcualte_distance:
                                    ax_dist.annotate(f'{chain_id}-{_hf.amino_d[res_name]}{res_id}', (values[0]+0.2, values[1]-0.25), fontsize=self.plot_parameters['node_label_size'])
                        else:
                            l = f'{chain_id}-{res_name}{res_id}' if self.plot_parameters['show_chain_label'] else f'{res_name}{res_id}'
                            ax.annotate(l, (values[0]+0.2, values[1]-0.25), fontsize=self.plot_parameters['node_label_size'], color=self.plot_parameters['non_prot_color'])
                            if n in conserved_nodes:
                                ax_cons.annotate(l, (values[0]+0.2, values[1]-0.25), fontsize=self.plot_parameters['node_label_size'], color=self.plot_parameters['non_prot_color'])
                                if calcualte_distance:
                                    ax_dist.annotate(l, (values[0]+0.2, values[1]-0.25), fontsize=self.plot_parameters['node_label_size'], color=self.plot_parameters['non_prot_color'])

                if color_info:
                    cbar = fig_cons.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax_cons)
                    cbar.ax.tick_params(labelsize=self.plot_parameters['plot_tick_fontsize'])
                    cbar.set_label(label=color_bar_label, size=self.plot_parameters['plot_label_fontsize'])

                ax.text(0.97, 0.99, self.name1, color=color1, fontsize=20, transform=ax.transAxes, ha='center', va='center')
                ax.text(0.97, 0.97, self.name2, color=color2, fontsize=20, transform=ax.transAxes, ha='center', va='center')

                plt.tight_layout()
                is_label = '_labeled' if label_nodes else ''
                is_propka = '_pKa_color' if color_info and color_propka else ''
                is_conservation = '_data_color' if color_info and color_data else ''
                is_backbone = '_backbone' if self.include_backbone_sidechain else ''
                if self.graph_type == 'hbond':
                    for form in self.plot_parameters['formats']:
                        fig.savefig(f'{self.compare_folder}compare_H-bond_graph_{self.name1}_with_{self.name2}{is_backbone}{is_label}.{form}', format=form, dpi=self.plot_parameters['plot_resolution'])
                        fig_cons.savefig(f'{self.compare_folder}conserved_H-bond_graph_{self.name1}_with_{self.name2}{is_propka}{is_conservation}{is_backbone}{is_label}.{form}', format=form, dpi=self.plot_parameters['plot_resolution'])
                        if calcualte_distance:
                            fig_dist.savefig(f'{self.compare_folder}conserved_H-bond_graph_{self.name1}_with_{self.name2}_distances{is_backbone}{is_label}.{form}', format=form, dpi=self.plot_parameters['plot_resolution'])
                    if is_label:
                        _hf.write_text_file(self.compare_folder+'compare_H-bond_graph_'+self.name1+'_with_'+self.name2+is_backbone+'_info.txt',
                            ['H-bond graph comparison of '+self.name1+' with '+self.name2,
                            '\nSelection string: '+str(self.selection[0:-15]),
                            '\n',
                            '\nNumber of nodes in '+self.name1+': '+str(len(graph1.nodes)),
                            '\nNumber of edges in '+self.name1+': '+str(len(graph1.edges)),
                            '\n',
                            '\nNumber of nodes in '+self.name2+': '+str(len(graph2.nodes)),
                            '\nNumber of edges in '+self.name2+': '+str(len(graph2.edges)),
                            '\n',
                            '\nNumber of nodes in both structures: '+str(len(conserved_nodes)),
                            '\nNumber of edges in both structures: '+str(len(conserved_edges)),
                            '\n',
                            '\nList of nodes in both structures: '+str(conserved_nodes),
                            '\n',
                            '\nList of edges in both structures: '+str(conserved_edges),
                            ])
                elif self.graph_type == 'water_wire':
                    waters = '_max_'+str(self.max_water)+'_water_bridges' if self.max_water > 0 else ''
                    occ = '_min_occupancy_'+str(occupancy) if occupancy  else ''
                    for form in self.plot_parameters['formats']:
                        fig.savefig(f'{self.compare_folder}compare{waters}{occ}_graph_{self.name1}_with_{self.name2}{is_label}.{form}', format=form, dpi=self.plot_parameters['plot_resolution'])
                        fig_cons.savefig(f'{self.compare_folder}conserved{waters}{occ}_graph_{self.name1}_with_{self.name2}{is_propka}{is_conservation}{is_label}.{form}', format=form, dpi=self.plot_parameters['plot_resolution'])
                    if is_label:
                        _hf.write_text_file(self.compare_folder+'compare'+waters+occ+'_graph_'+self.name1+'_with_'+self.name2+'_info.txt',
                            ['Water wire graph comparison of '+self.name1+' with '+self.name2,
                            '\nSelection string: '+str(self.selection[0:-15]),
                            '\nNumber of maximum water molecules allowed in the bridge: '+str(self.max_water),
                            '\nMinimum H-bond occupancy: '+str(occupancy) if occupancy  else '',
                            '\n',
                            '\nNumber of nodes in '+self.name1+': '+str(len(graph1.nodes)),
                            '\nNumber of edges in '+self.name1+': '+str(len(graph1.edges)),
                            '\n',
                            '\nNumber of nodes in '+self.name2+': '+str(len(graph2.nodes)),
                            '\nNumber of edges in '+self.name2+': '+str(len(graph2.edges)),
                            '\n',
                            '\nNumber of nodes in both structures: '+str(len(conserved_nodes)),
                            '\nNumber of edges in both structures: '+str(len(conserved_edges)),
                            '\n',
                            '\nList of nodes in both structures: '+str(conserved_nodes),
                            '\n',
                            '\nList of edges in both structures: '+str(conserved_edges),
                            ])
                plt.close()
            else: self.logger.warning('No '+self.graph_type+' graph for both structures. Comparison can not be performed.')


    def _custom_node_coloring(self, objects, name, color_propka, color_data, node_color_selection):
        color_info = {}
        if  color_propka and color_data:
            self.logger.info(f'Can not color plot by propka and external data values at the same time. Please select just one coloring option!')
        elif color_propka:
            struct_object = objects['structure'] if not hasattr(self, 'occupancy') else objects['mda']
            selected_nodes = struct_object.select_atoms(str(node_color_selection))
            try:
                color_info = _hf.read_propka_file(f'{self.target_folder}/{name}.pka', selected_nodes)
            except:
                self.logger.info(f"{name}.pka not found. To color residues by pKa values, place the propka file in the PDB folder, next to the PDB file.")
            try:
                len(color_info)
            except:
                self.logger.info(f'{name}.propka does not contain the selected residues. Please update the Residues to color!' )
        elif color_data:
            struct_object = objects['structure'] if not hasattr(self, 'occupancy') else objects['mda']
            selected_nodes = struct_object.select_atoms(str(node_color_selection))
            color_info = _hf.read_color_data_file(name, self.target_folder, selected_nodes)
            try:
                len(color_info)
            except:
                self.logger.error(f"No {name}_data .txt file was found in {self.target_folder} or content is invalid. To enable coloring by data values please add a corresponding file.")
        return color_info

    def _plot_dist_plot(self, dist_plot_data, ax, fig, color_map, label_edges):
        import matplotlib.patheffects as pe

        distances = np.array([point[2] for point in dist_plot_data])
        cmap = cm.get_cmap(color_map, len(distances))

        scaled_values = (distances - distances.min()) / (distances.max() - distances.min())
        value_colors = [self.plot_parameters['graph_color'] if round(distances[i], 1) == 0 else cmap(scaled_values[i]) for i in range(len(distances))]

        for i, point in enumerate(dist_plot_data):
            x, y = point[0], point[1]

            ax.plot(x, y, color=value_colors[i], marker='o', linewidth=self.plot_parameters['edge_width'], markersize=self.plot_parameters['node_size']*0.01, markerfacecolor= self.plot_parameters['graph_color'], markeredgecolor=value_colors[i],
                path_effects=[pe.Stroke(linewidth=2, foreground=self.plot_parameters['graph_color']), pe.Normal()])

            if label_edges:
                ax.annotate(round(point[2], 1), (x[0] + (x[1]-x[0])/2, y[0] + (y[1]-1.0-y[0])/2), color='blue',  fontsize=self.plot_parameters['edge_label_size'])

        max_val = np.max([np.abs( distances.min()), np.abs(distances.max())])
        norm = mpl.colors.Normalize(vmin=-1*max_val, vmax=max_val)

        cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)
        cbar.ax.tick_params(labelsize=self.plot_parameters['plot_tick_fontsize'])
        cbar.set_label(label='H-bond distance change', size=self.plot_parameters['plot_label_fontsize'])

        return ax, fig
