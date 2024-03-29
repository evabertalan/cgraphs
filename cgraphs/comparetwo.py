from . import helperfunctions as _hf
import shutil
import numpy as np
import copy
import MDAnalysis as _mda
from .proteingraphanalyser import ProteinGraphAnalyser
import matplotlib.pyplot as plt


class CompareTwo(ProteinGraphAnalyser):
    def __init__(self, type_option, pdb1=None, pdb2=None, psf1=None, psf2=None, dcd1=None, dcd2=None, target_folder='', name1=None, name2=None):
        self.type_option = type_option
        if self.type_option == 'pdb' and pdb1 and pdb2:
            self.name1 = _hf.retrieve_pdb_code(pdb1, '.pdb')
            self.name2 = _hf.retrieve_pdb_code(pdb2, '.pdb')
            self.compare_folder = _hf.create_directory(target_folder+'/workfolder/compare_'+self.name1+'_'+self.name2)+'/'
            shutil.copy(pdb1, self.compare_folder+self.name1+'.pdb')
            shutil.copy(pdb2, self.compare_folder+self.name2+'.pdb')

            ProteinGraphAnalyser.__init__(self, pdb_root_folder=self.compare_folder, target_folder=target_folder, reference_pdb=pdb1)
            self.logger.info('COMPARE STRUCTURES '+ self.name1 + ' WITH ' + self.name2)
            ProteinGraphAnalyser.align_structures(self, superimposition_threshold=30)

        elif self.type_option == 'dcd' and psf1 and psf2 and dcd1 and dcd2:
            self.name1, self.name2 = name1, name2
            self.compare_folder = _hf.create_directory(target_folder+'/workfolder/compare_'+self.name1+'_'+self.name2)+'/'

            ProteinGraphAnalyser.__init__(self, pdb_root_folder=self.compare_folder, target_folder=target_folder, type_option='dcd', psf_files=[psf1, psf2], dcd_files=[dcd1, dcd2], sim_names=[name1, name2])
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


    def plot_graph_comparison(self, color1='blue', color2='green', label_nodes=True, label_edges=True, xlabel='PCA projected xy plane', ylabel='Z coordinates ($\AA$)'):
        if len(self.graph_coord_objects.items()) != 2: self.logger.warning('There are '+str(len(self.graph_coord_objects.items()))+' structures selected. Graph comparison is possible for exactly two structures.')
        else:
            self.logger.info('Plot comparison '+self.graph_type+' graph for '+ self.name1 + ' with ' + self.name2+str(' with labels' if label_nodes else ''))
            if 'graph' in self.graph_coord_objects[self.name1].keys() and 'graph' in self.graph_coord_objects[self.name2].keys():
                if hasattr(self, 'occupancy'):
                    occupancy = self.occupancy
                    wba1 = copy.deepcopy(self.graph_coord_objects[self.name1]['wba'])
                    wba1.filter_occupancy(occupancy)
                    graph1 = wba1.filtered_graph

                    wba2 = copy.deepcopy(self.graph_coord_objects[self.name1]['wba'])
                    wba2.filter_occupancy(occupancy)
                    graph2 = wba2.filtered_graph
                else:
                    occupancy = None
                    graph1 = self.graph_coord_objects[self.name1]['graph']
                    graph2 = self.graph_coord_objects[self.name2]['graph']

                pos1 = self._get_node_positions(self.graph_coord_objects[self.name1], pca=False)
                pos2 = self._get_node_positions(self.graph_coord_objects[self.name2], pca=False)
                conserved_edges=[]
                conserved_nodes=[]

                all_pos = {}
                for i, pos in enumerate([pos1, pos2]):
                    for key, value in pos.items():
                        if key.split('-')[1] in ['HOH', 'TIP3']: key=str(i+1)+'-'+key.split('-')[1]+'-'+key.split('-')[2]
                        if key not in all_pos.keys(): all_pos.update({ key:value })

                node_pca_pos = _hf.calculate_pca_positions(all_pos)

                plot_name = 'H-bond' if self.graph_type == 'hbond' else 'Water wire'
                fig, ax = _hf.create_plot(title=plot_name+' graph comparison of '+self.name1+' with '+self.name2, xlabel=xlabel, ylabel=ylabel)

                for e in graph1.edges:
                    e0 = e[0] if e[0].split('-')[1] not in ['HOH', 'TIP3'] else '1-'+e[0].split('-')[1]+'-'+e[0].split('-')[2]
                    e1 = e[1] if e[1].split('-')[1] not in ['HOH', 'TIP3'] else '1-'+e[1].split('-')[1]+'-'+e[1].split('-')[2]
                    if _hf.is_conserved_edge(np.array([[e2[0], e2[1]] for e2 in graph2.edges]), e0, e1):
                        color = 'gray'
                        conserved_edges.append(e)
                    else: color = color1
                    if e0 in node_pca_pos.keys() and e1 in node_pca_pos.keys():
                        edge_line = [node_pca_pos[e0], node_pca_pos[e1]]
                        x=[edge_line[0][0], edge_line[1][0]]
                        y=[edge_line[0][1], edge_line[1][1]]
                        ax.plot(x, y, color=color, marker='o', linewidth=2, markersize=15, markerfacecolor=color, markeredgecolor=color)

                for e in graph2.edges:
                    e0 = e[0] if e[0].split('-')[1] not in ['HOH', 'TIP3'] else '2-'+e[0].split('-')[1]+'-'+e[0].split('-')[2]
                    e1 = e[1] if e[1].split('-')[1] not in ['HOH', 'TIP3'] else '2-'+e[1].split('-')[1]+'-'+e[1].split('-')[2]

                    if e0 in node_pca_pos.keys() and e1 in node_pca_pos.keys():
                        if not _hf.is_conserved_edge(np.array([[e2[0], e2[1]] for e2 in graph1.edges]), e0, e1):
                            edge_line = [node_pca_pos[e0], node_pca_pos[e1]]
                            x=[edge_line[0][0], edge_line[1][0]]
                            y=[edge_line[0][1], edge_line[1][1]]
                            ax.plot(x, y, color=color2, marker='o', linewidth=2, markersize=15, markerfacecolor=color2, markeredgecolor=color2)

                for n in graph1.nodes:
                    n = n if n.split('-')[1] not in ['HOH', 'TIP3'] else '1-'+n.split('-')[1]+'-'+n.split('-')[2]
                    if n in node_pca_pos.keys():
                        if n in graph2.nodes:
                            color = 'gray'
                            conserved_nodes.append(n)
                        else: color = color1
                        if n.split('-')[1] in ['HOH', 'TIP3']:
                            ax.scatter(node_pca_pos[n][0], node_pca_pos[n][1],color='#db5c5c', s=150, zorder=5, edgecolors=color)
                        else: ax.scatter(node_pca_pos[n][0], node_pca_pos[n][1], s=200, color=color, zorder=5)

                for n in graph2.nodes:
                    n = n if n.split('-')[1] not in ['HOH', 'TIP3'] else '2-'+n.split('-')[1]+'-'+n.split('-')[2]
                    if n in node_pca_pos.keys():
                        if n not in graph1.nodes:
                            if n.split('-')[1] in ['HOH', 'TIP3']:
                                ax.scatter(node_pca_pos[n][0], node_pca_pos[n][1], color='#db5c5c', s=150, zorder=5, edgecolors=color2)
                            else: ax.scatter(node_pca_pos[n][0], node_pca_pos[n][1], s=200, color=color2, zorder=5)

                if label_nodes:
                    for n in all_pos.keys():
                        values = node_pca_pos[n]
                        chain_id, res_name, res_id = _hf.get_node_name_pats(n)
                        if res_name in ['HOH', 'TIP3']:
                            ax.annotate(f'W{res_id}', (values[0]+0.2, values[1]-0.25), fontsize=12)
                        elif res_name in _hf.amino_d.keys():
                            ax.annotate(f'{chain_id}-{_hf.amino_d[res_name]}{res_id}', (values[0]+0.2, values[1]-0.25), fontsize=12)
                        else:
                            ax.annotate(f'{chain_id}-{res_name}{res_id}', (values[0]+0.2, values[1]-0.25), fontsize=12, color='blue')

                ax.text(0.97, 0.99, self.name1, color=color1, fontsize=20, transform=ax.transAxes, ha='center', va='center')
                ax.text(0.97, 0.97, self.name2, color=color2, fontsize=20, transform=ax.transAxes, ha='center', va='center')

                plt.tight_layout()
                is_label = '_labeled' if label_nodes else ''
                if self.graph_type == 'hbond':
                    plt.savefig(self.compare_folder+'compare_H-bond_graph_'+self.name1+'_with_'+self.name2+is_label+'.png')
                    plt.savefig(self.compare_folder+'compare_H-bond_graph_'+self.name1+'_with_'+self.name2+is_label+'.eps', format='eps')
                    if is_label:
                        _hf.write_text_file(self.compare_folder+'compare_H-bond_graph_'+self.name1+'_with_'+self.name2+'_info.txt',
                            ['H-bond graph comparison of '+self.name1+' with '+self.name2,
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
                    plt.savefig(self.compare_folder+'compare'+waters+occ+'_graph_'+self.name1+'_with_'+self.name2+is_label+'.png')
                    plt.savefig(self.compare_folder+'compare'+waters+occ+'_graph_'+self.name1+'_with_'+self.name2+is_label+'.eps', format='eps')
                    if is_label:
                        _hf.write_text_file(self.compare_folder+'compare'+waters+occ+'_graph_'+self.name1+'_with_'+self.name2+'_info.txt',
                            ['Water wire graph comparison of '+self.name1+' with '+self.name2,
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
