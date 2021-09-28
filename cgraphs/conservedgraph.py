from . import helperfunctions as _hf
import numpy as np
import copy
from .proteingraphanalyser import ProteinGraphAnalyser
import matplotlib.pyplot as plt


class ConservedGraph(ProteinGraphAnalyser):
    def __init__(self, pdb_root_folder='',  type_option='pdb', target_folder='', reference_pdb='', reference_coordinates=None, sequance_identity_threshold=0.75, psf_files=[], dcd_files=[[]], sim_names=[]):
        if type_option == 'pdb':
            ProteinGraphAnalyser.__init__(self, pdb_root_folder, target_folder, reference_pdb)
            self.logger.info('CONSERVED NETWORK ANALYSIS')
            ProteinGraphAnalyser.align_structures(self, sequance_identity_threshold=sequance_identity_threshold, superimposition_threshold=30)
            if reference_coordinates is not None:
                self.reference_coordinates = reference_coordinates
                self.logger.info('Using water cluster coordinates as conserved water molecules.')
            # self.pca_positions = _hf.calculate_pca_positions(self.reference_coordinates)

        elif type_option == 'dcd':
            ProteinGraphAnalyser.__init__(self, target_folder=target_folder, type_option='dcd', psf_files=psf_files, dcd_files=dcd_files, sim_names=sim_names)

        else: raise ValueError('Given type_option should be "pdb" or "dcd"')

    def get_conserved_graph(self, conservation_threshold=0.9, occupancy=None, eps=1.5):
        self.logger.info('Conservation threshold across structures is set to: '+str(conservation_threshold*100)+'%')
        if occupancy: self.logger.info('H-bond occupancy is set to: '+str(occupancy*100)+'%' )
        self.occupancy = occupancy
        nodes = []
        edges = []
        self.avg_water_per_conserved_edges = None
        avg_water_per_edge = {}
        if self.graph_type == 'water_wire':
            for objects in self.graph_coord_objects.values():
                if 'graph' in objects.keys():
                    avg_waters = objects['wba'].compute_average_water_per_wire()
                    if occupancy:
                        wba = copy.deepcopy(objects['wba'])
                        wba.filter_occupancy(occupancy)
                        graph = wba.filtered_graph
                    else: graph = objects['graph']
                    for node in graph.nodes:
                        node = _hf.get_node_name(node)
                        nodes.append(node)

                    for edge in graph.edges:
                        e0 =  _hf.get_node_name(edge[0])
                        e1 =  _hf.get_node_name(edge[1])
                        if ([e1, e0]) in edges:
                            edges.append([e1, e0])
                        else: edges.append([e0, e1])

                        key = e0+':'+e1
                        key2 = e1+':'+e0
                        if key in avg_waters:
                            if key in avg_water_per_edge: avg_water_per_edge[key].append(avg_waters[key])
                            elif key2 in avg_water_per_edge: avg_water_per_edge[key2].append(avg_waters[key])
                            else:
                                avg_water_per_edge.update( {key: [avg_waters[key]]} )
                        elif key2 in avg_waters:
                            if key in avg_water_per_edge: avg_water_per_edge[key].append(avg_waters[key2])
                            elif key2 in avg_water_per_edge: avg_water_per_edge[key2].append(avg_waters[key2])
                            else:
                                avg_water_per_edge.update( {key2: [avg_waters[key2]]} )

        elif self.graph_type == 'hbond':
            for objects in self.graph_coord_objects.values():
                if 'graph' in objects.keys():
                    graph = objects['graph']
                    waters = {}
                    for res in objects['structure'].select_atoms('(protein and name CA) or'+_hf.water_def):
                        if res.resname in ['HOH', 'TIP3']: waters[res.resid]=res.position
    #             for graph in self.graphs:
                    #here select conserved water by the superimposed ones
    #                 if useEstimatedConservedWaters:
    #                     check whether water coordinates belong to any clustrer thatn give name
                    for node in graph.nodes:
                        nodes.append(_hf.get_node_name(node))
                    for edge in graph.edges:
                        edge = [edge[0], edge[1]]
                        _i = [0, 1]
                        for i in _i:
                            if edge[i].split('-')[1] in ['HOH', 'TIP3']:
                                #TODO: fix issue  with water ID in the graph
                                # if int(edge[i].split('-')[2]) >= 10000: n = int(edge[i].split('-')[2])-10000
                                # else:
                                n =  edge[i].split('-')[2]
                                for key, cc in self.reference_coordinates.items():
                                    if key.startswith('X-w'):
                                        w = waters.get(int(n))
                                        #TODO change radius regarding EPS
                                        if w is not None and ((w[0]-cc[0])**2 + (w[1]-cc[1])**2 + (w[2]-cc[2])**2 <= float(eps)**2):
                                            edge[i] = 'X-w-'+key.split('-')[-1]
                        e0 = _hf.get_node_name(edge[0])
                        e1 = _hf.get_node_name(edge[1])
                        if ([e1, e0]) in edges:
                            edges.append([e1, e0])
                        else: edges.append([e0, e1])

        th = np.round(len(self.graph_coord_objects) * conservation_threshold)
        u_nodes, c_nodes = np.unique(nodes, return_counts=True)
        self.conserved_nodes = u_nodes[np.where(c_nodes >= th)[0]]
        u_edges, c_edges = np.unique(edges, return_counts=True, axis=0)
        self.conserved_edges = u_edges[np.where(c_edges >= th)[0]]
        if len(avg_water_per_edge) > 0:
            self.avg_water_per_conserved_edges = {key: np.mean(value) for key, value in avg_water_per_edge.items() if len(value) >= th }


    def plot_conserved_graph(self, label_nodes=True, label_edges=True, xlabel='PCA projected xy plane', ylabel='Z coordinates ($\AA$)'):
        self.logger.info('Plotting conserved '+self.graph_type+' graph'+str(' with labels' if label_nodes else ''))
        self.pca_positions = _hf.calculate_pca_positions(self.reference_coordinates)
        #TODO set back lables
        plot_name = 'H-bond' if self.graph_type == 'hbond' else 'water wire'
        fig, ax = _hf.create_plot(title='Conserved '+plot_name+' graph',
                                  xlabel=xlabel,
                                  ylabel=ylabel)
        for e in self.conserved_edges:
            if e[0] in self.pca_positions.keys() and e[1] in self.pca_positions.keys():
                edge_line = [self.pca_positions[e[0]], self.pca_positions[e[1]]]
                x=[edge_line[0][0], edge_line[1][0]]
                y=[edge_line[0][1], edge_line[1][1]]
                ax.plot(x, y, color='gray', marker='o', linewidth=2, markersize=15, markerfacecolor='gray', markeredgecolor='gray')
            if label_edges and self.avg_water_per_conserved_edges:
                key1 = e[0]+':'+e[1]
                key2 = e[1]+':'+e[0]
                water = [value for key, value in self.avg_water_per_conserved_edges.items() if key == key1 or key == key2][0]
                ax.annotate(np.round(water,1), (x[0] + (x[1]-x[0])/2, y[0] + (y[1]-y[0])/2), color='indianred',  fontsize=10, weight='bold',)

        for n in self.conserved_nodes:
            if n in self.pca_positions.keys():
                color = '#db5c5c' if n.split('-')[1] in ['HOH', 'TIP3'] else 'gray'
                ax.scatter(self.pca_positions[n][0], self.pca_positions[n][1], color=color, s=200, zorder=5)

        if self.graph_type == 'hbond':
            for r in self.reference_coordinates:
                if r.split('-')[1].startswith('w'):
                    ax.scatter(self.pca_positions[r][0], self.pca_positions[r][1], color='#db5c5c', s=80, zorder=5)
                    if label_nodes: ax.annotate('W'+r.split('-')[-1], (self.pca_positions[r][0]+0.2, self.pca_positions[r][1]-0.25), fontsize=13, zorder=6)


        if label_nodes:
            for node in self.conserved_nodes:
                chain_id, res_name, res_id = _hf.get_node_name_pats(node)
                if node in self.pca_positions.keys():
                    if res_name not in ['HOH', 'TIP3'] and res_name in _hf.amino_d.keys():
                        ax.annotate(f'{chain_id}-{_hf.amino_d[res_name]}{res_id}', (self.pca_positions[node][0]+0.2, self.pca_positions[node][1]-0.25), fontsize=13, zorder=6)
                    elif res_name not in ['HOH', 'TIP3'] and res_name not in _hf.amino_d.keys():
                        ax.annotate(f'{chain_id}-{res_name}{res_id}', (self.pca_positions[node][0]+0.2, self.pca_positions[node][1]-0.25), fontsize=13, zorder=6, color='blue')

        plt.tight_layout()
        is_label = '_labeled' if label_nodes else ''
        if self.graph_type == 'hbond':
            plot_folder = _hf.create_directory(self.workfolder+'/H-bond_graphs/')
            plt.savefig(plot_folder+'conserved_H-bond_graph'+is_label+'.png')
            plt.savefig(plot_folder+'conserved_H-bond_graph'+is_label+'.eps', format='eps')
            if is_label:
                _hf.write_text_file(plot_folder+'conserved_H-bond_graph_info.txt',
                    ['Conserved H-bond graph of '+str(len(self.graph_coord_objects.keys()))+' PDB structures',
                    '\n',
                    '\nNumber of conserved nodes : '+str(len(self.conserved_nodes)),
                    '\nNumber of conserved edges : '+str(len(self.conserved_edges)),
                    '\n',
                    '\nList of conserved nodes: '+str(self.conserved_nodes),
                    '\n',
                    '\nList of conserved edges: '+str(self.conserved_edges),
                    ])
        elif self.graph_type == 'water_wire':
            plot_folder = _hf.create_directory(self.workfolder+'/'+str(self.max_water)+'_water_wires/')
            waters = '_max_'+str(self.max_water)+'_water_bridges' if self.max_water > 0 else ''
            occ = '_min_occupancy_'+str(self.occupancy) if self.occupancy  else ''
            plt.savefig(plot_folder+'conserved'+waters+occ+'_graph'+is_label+'.png')
            plt.savefig(plot_folder+'conserved'+waters+occ+'_graph'+is_label+'.eps', format='eps')
            if is_label:
                _hf.write_text_file(plot_folder+'conserved'+waters+occ+'_graph_inof.txt',
                    ['Conserved water wire graph of '+str(len(self.graph_coord_objects.keys()))+
                    str(' PDB structures' if not self.occupancy else ' simulations'),
                    '\nNumber of maximum water molecules allowed in the bridge: '+str(self.max_water),
                    '\nMinimum H-bond occupancy: '+str(self.occupancy) if self.occupancy  else '',
                    '\n',
                    '\nNumber of conserved nodes : '+str(len(self.conserved_nodes)),
                    '\nNumber of conserved edges : '+str(len(self.conserved_edges)),
                    '\n',
                    '\nList of conserved nodes: '+str(self.conserved_nodes),
                    '\n',
                    '\nList of conserved edges: '+str(self.conserved_edges),
                    ])
        plt.close()


    def plot_difference(self, label_nodes=True, label_edges=True, xlabel='PCA projected xy plane', ylabel='Z coordinates ($\AA$)'):
        self.logger.info('Plotting difference '+self.graph_type+' graphs'+str(' with labels' if label_nodes else ''))
        for name, objects in self.graph_coord_objects.items():
            if 'graph' in objects.keys():
                if self.occupancy:
                    wba = copy.deepcopy(objects['wba'])
                    wba.filter_occupancy(self.occupancy)
                    graph = wba.filtered_graph
                else: graph = objects['graph']

                self.logger.debug('Calculating '+self.graph_type+' difference graph for: '+name)
                plot_name = 'H-bond' if self.graph_type == 'hbond' else 'water wire'
                fig, ax = _hf.create_plot(title='Difference '+plot_name+' graph of structure '+name,
                                          xlabel=xlabel,
                                          ylabel=ylabel)
                node_pca_pos = self._get_node_positions(objects)
                node_pca_pos = _hf.check_projection_sign(node_pca_pos, self.pca_positions)

                for e in graph.edges:
                    e0 = _hf.get_node_name(e[0])
                    e1 = _hf.get_node_name(e[1])
                    if e0 in node_pca_pos.keys() and e1 in node_pca_pos.keys():
                        edge_line = [node_pca_pos[e0], node_pca_pos[e1]]
                        x=[edge_line[0][0], edge_line[1][0]]
                        y=[edge_line[0][1], edge_line[1][1]]

                        if _hf.is_conserved_edge(self.conserved_edges, e0, e1):
                            ax.plot(x, y, color='gray', marker='o', linewidth=2, markersize=15, markerfacecolor='gray', markeredgecolor='gray')
                        else:
                            ax.plot(x, y, color='#129fe6', marker='o', linewidth=2, markersize=15, markerfacecolor='#129fe6', markeredgecolor='#129fe6')
                        if label_edges and self.graph_type == 'water_wire':
                            waters, occ_per_wire, _ = _hf.get_edge_params(objects['wba'], graph.edges)
                            ax.annotate(np.round(waters[list(graph.edges).index(e)],1), (x[0] + (x[1]-x[0])/2, y[0] + (y[1]-y[0])/2), color='indianred',  fontsize=10, weight='bold',)
                            ax.annotate(int(occ_per_wire[list(graph.edges).index(e)]*100), (x[0] + (x[1]-x[0])/2, y[0] + (y[1]-1.0-y[0])/2), color='green',  fontsize=10)

                for node in graph.nodes:
                    n = _hf.get_node_name(node)
                    if n in node_pca_pos.keys():
                        if n in self.conserved_nodes:
                            ax.scatter(node_pca_pos[n][0], node_pca_pos[n][1], s=200, color='gray', zorder=5)
                        else: ax.scatter(node_pca_pos[n][0], node_pca_pos[n][1], s=200, color='orange')

                if self.graph_type == 'hbond':
                    for n, values in node_pca_pos.items():
                        if n.split('-')[1] in ['HOH', 'TIP3']:
                            ax.scatter(values[0],values[1], color='#db5c5c', s=110, zorder=5)

                if label_nodes:
                    for n in graph.nodes:
                        n = _hf.get_node_name(n)
                        if n in node_pca_pos.keys():
                            values = node_pca_pos[n]
                            chain_id, res_name, res_id = _hf.get_node_name_pats(n)
                            if res_name in ['HOH', 'TIP3']: ax.annotate(f'W{res_id}', (values[0]+0.2, values[1]-0.25), fontsize=12)
                            elif res_name in _hf.amino_d.keys():
                                ax.annotate(f'{chain_id}-{_hf.amino_d[res_name]}{res_id}', (values[0]+0.2, values[1]-0.25), fontsize=12)
                            else: ax.annotate(f'{chain_id}-{res_name}{res_id}', (values[0]+0.2, values[1]-0.25), fontsize=12, color='blue')

                plt.tight_layout()
                is_label = '_labeled' if label_nodes else ''
                if self.graph_type == 'hbond':
                    plot_folder = _hf.create_directory(self.workfolder+'/H-bond_graphs/'+name+'/')
                    plt.savefig(plot_folder+name+'_H-bond_difference_graph'+is_label+'.png')
                    plt.savefig(plot_folder+name+'_H-bond_difference_graph'+is_label+'.eps', format='eps')
                elif self.graph_type == 'water_wire':
                    plot_folder = _hf.create_directory(self.workfolder+'/'+str(self.max_water)+'_water_wires/'+name+'/')
                    waters = '_max_'+str(self.max_water)+'_water_bridges' if self.max_water > 0 else ''
                    occ = '_min_occupancy_'+str(self.occupancy) if self.occupancy  else ''
                    plt.savefig(plot_folder+name+waters+occ+'_difference_graph'+is_label+'.png')
                    plt.savefig(plot_folder+name+waters+occ+'_difference_graph'+is_label+'.eps', format='eps')
                plt.close()


    def plot_unique_grpahs(self):
        pass

    def plot_unique_linear_lenght(self):
        pass
