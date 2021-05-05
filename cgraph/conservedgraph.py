from . import helperfunctions as _hf
import numpy as np
import copy
from .proteingraphanalyser import ProteinGraphAnalyser
import matplotlib.pyplot as plt


class ConservedGraph(ProteinGraphAnalyser):
    def __init__(self, pdb_root_folder='',  type_option='pdb', target_folder='', reference_pdb='', reference_coordinates=None, sequance_identity_threshold=0.75):
        if type_option == 'pdb':
            ProteinGraphAnalyser.__init__(self, pdb_root_folder, target_folder, reference_pdb)
            self.logger.info('CONSERVED NETWORK ANALYSIS')
            ProteinGraphAnalyser.align_structures(self, sequance_identity_threshold=sequance_identity_threshold)
            if reference_coordinates is not None:
                self.reference_coordinates = reference_coordinates
                self.logger.info('Using water cluster coordinates as conserved water molecules.')
            self.pca_positions = _hf.calculate_pca_positions(self.reference_coordinates)

        elif type_option == 'dcd':
            ProteinGraphAnalyser.__init__(self, target_folder=target_folder, type_option='dcd')

        else: raise ValueError('Given type_option should be "pdb" or "dcd"')

    @staticmethod
    def get_conserved_graph(self, conservation_threshold=0.9, occupancy=None):
        self.logger.info('Conservation threshold across structures is set to: '+str(conservation_threshold*100)+'%')
        if occupancy: self.logger.info('H-bond occupancy is set to: '+str(occupancy*100)+'%' )
        self.occupancy = occupancy
        nodes = []
        edges = []
        if self.graph_type == 'water_wire':
            for objects in self.graph_coord_objects.values():
                if 'graph' in objects.keys():
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

        elif self.graph_type == 'hbond':
            for objects in self.graph_coord_objects.values():
                if 'graph' in objects.keys():
                    graph = objects['graph']
                    waters = {}
                    for res in list(objects['structure'][0].get_residues()):
                        if res.get_id()[0] == 'W': waters[res.get_id()[1]]=res['O'].get_coord()
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
                            if edge[i].split('-')[1] == 'HOH':
                                #TODO: fix issue  with water ID in the graph
                                if int(edge[i].split('-')[2]) >= 10000: n = int(edge[i].split('-')[2])-10000
                                else: n =  edge[i].split('-')[2]
                                for key, cc in self.reference_coordinates.items():
                                    if key.startswith('w'):
                                        w = waters.get(int(n))
                                        #TODO change radius regarding EPS
                                        if ((w[0]-cc[0])**2 + (w[1]-cc[1])**2 + (w[2]-cc[2])**2 <= 1.5**2):
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

    def plot_conserved_graph(self, label_nodes=True, label_edges=True, xlabel='PCA projected xy plane', ylabel='Z coordinates ($\AA$)'):
        self.logger.info('Plotting conserved '+self.graph_type+' graph'+str(' with labels' if label_nodes else ''))
        #TODO set back lables
        plot_name = 'H-bond' if self.graph_type == 'hbond' else 'water wire'
        fig, ax = _hf.create_plot(title='Conserved '+plot_name+' graph',
                                  xlabel=xlabel,
                                  ylabel=ylabel)
        for e in self.conserved_edges:
            edge_line = [self.pca_positions[e[0]], self.pca_positions[e[1]]]
            x=[edge_line[0][0], edge_line[1][0]]
            y=[edge_line[0][1], edge_line[1][1]]
            ax.plot(x, y, color='gray', marker='o', linewidth=2, markersize=18, markerfacecolor='gray', markeredgecolor='gray')

        for n in self.conserved_nodes:
            ax.scatter(self.pca_positions[n][0], self.pca_positions[n][1], color='gray', s=180, zorder=5)

        if self.graph_type == 'hbond':
            for r in self.reference_coordinates:
                if r.split('-')[1].startswith('w'):
                    ax.scatter(self.pca_positions[r][0], self.pca_positions[r][1], color='#db5c5c', s=80, zorder=5)
                    if label_nodes: ax.annotate('W'+r.split('-')[-1], (self.pca_positions[r][0]+0.2, self.pca_positions[r][1]-0.25), fontsize=17, zorder=6)


        if label_nodes:
            for node in self.conserved_nodes:
                if node.split('-')[1] != 'HOH':
                    ax.annotate(str(node.split('-')[0])+'-'+str(_hf.amino_d[node.split('-')[1]])+str(int(node.split('-')[2])), (self.pca_positions[node][0]+0.2, self.pca_positions[node][1]-0.25), fontsize=17, zorder=6)
        plt.tight_layout()
        is_label = '_labeled' if label_nodes else ''
        if self.graph_type == 'hbond':
            plot_folder = _hf.create_directory(self.workfolder+'/H-bond_graphs/')
            plt.savefig(plot_folder+'conserved_H-bond_graph'+is_label+'.png')
            plt.savefig(plot_folder+'conserved_H-bond_graph'+is_label+'.eps', format='eps')
        elif self.graph_type == 'water_wire':
            plot_folder = _hf.create_directory(self.workfolder+'/'+str(self.max_water)+'_water_wires/')
            waters = '_max_'+str(self.max_water)+'_water_bridges' if self.max_water > 0 else ''
            occ = '_min_occupancy_'+str(self.occupancy) if self.occupancy  else ''
            plt.savefig(plot_folder+'conserved'+waters+occ+'_graph'+is_label+'.png')
            plt.savefig(plot_folder+'conserved'+waters+occ+'_graph'+is_label+'.eps', format='eps')
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
                    edge_line = [node_pca_pos[e0], node_pca_pos[e1]]
                    x=[edge_line[0][0], edge_line[1][0]]
                    y=[edge_line[0][1], edge_line[1][1]]

                    if _hf.is_conserved_edge(self.conserved_edges, e0, e1):
                        ax.plot(x, y, color='gray', marker='o', linewidth=2, markersize=18, markerfacecolor='gray', markeredgecolor='gray')
                    else:
                        ax.plot(x, y, color='#129fe6', marker='o', linewidth=2, markersize=18, markerfacecolor='#129fe6', markeredgecolor='#129fe6')

                for node in graph.nodes:
                    n = _hf.get_node_name(node)
                    if n in self.conserved_nodes:
                        ax.scatter(node_pca_pos[n][0], node_pca_pos[n][1], s=180, color='gray')
                    else: ax.scatter(node_pca_pos[n][0], node_pca_pos[n][1], s=180, color='orange')

                if self.graph_type == 'hbond':
                    for n, values in node_pca_pos.items():
                        if n.split('-')[0] == 'HOH':
                            ax.scatter(values[0],values[1], color='#db5c5c', s=110, zorder=5)

                if label_nodes:
                    for n in graph.nodes:
                        n = _hf.get_node_name(n)
                        values = node_pca_pos[n]
                        if n.split('-')[1] == 'HOH': ax.annotate('W'+str(int(n.split('-')[2])), (values[0]+0.2, values[1]-0.25), fontsize=12)
                        else: ax.annotate(str(n.split('-')[0])+'-'+str(_hf.amino_d[n.split('-')[1]])+str(int(n.split('-')[2])), (values[0]+0.2, values[1]-0.25), fontsize=12)

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
