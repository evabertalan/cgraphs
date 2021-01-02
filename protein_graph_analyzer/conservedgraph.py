import helperfunctions as _hf
import numpy as np
from proteingraphanalyser import ProteinGraphAnalyser
import matplotlib.pyplot as plt
import waterclusters as wc


class ConservedGraph(ProteinGraphAnalyser):
    def __init__(self, pdb_root_folder, target_folder='', reference_pdb='', reference_coordinates=None):
        ProteinGraphAnalyser.__init__(self, pdb_root_folder, target_folder, reference_pdb)
        ProteinGraphAnalyser.align_structures(self)
        if reference_coordinates is not None: self.reference_coordinates = reference_coordinates
        self.pca_positions = _hf.calculate_pca_positions(self.reference_coordinates)

    
    def get_conserved_graph(self, threshold=0.9):
        nodes = []
        edges = []
        if self.graph_type == 'water_wire':
            for objects in self.graph_coord_objects.values():
                if 'graph' in objects.keys():
                    graph = objects['graph']
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
    #                         print(edge[i])
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
                        
        th = np.round(len(self.file_list) * threshold)
        u_nodes, c_nodes = np.unique(nodes, return_counts=True)
        self.conserved_nodes = u_nodes[np.where(c_nodes >= th)[0]]
        u_edges, c_edges = np.unique(edges, return_counts=True, axis=0)
        self.conserved_edges = u_edges[np.where(c_edges >= th)[0]]
    
    def plot_conserved_graph(self, label_nodes=True, label_edges=True):
        fig, ax = _hf.create_plot(title='Conserved '+self.graph_type+' graph',
                                  xlabel='PCA projected xy plane',
                                  ylabel='Z coordinates')
        for e in self.conserved_edges:
            edge_line = [self.pca_positions[e[0]], self.pca_positions[e[1]]]
            x=[edge_line[0][0], edge_line[1][0]]
            y=[edge_line[0][1], edge_line[1][1]]
            ax.plot(x, y, color='gray', marker='o', linewidth=2, markersize=13, markerfacecolor='gray', markeredgecolor='gray')
        
        if self.graph_type == 'hbond':
            for n in self.conserved_nodes:
                ax.scatter(self.pca_positions[n][0], self.pca_positions[n][1], color='gray', s=100, zorder=5)
            for r in self.reference_coordinates:
                if r.startswith('w'):
                    ax.scatter(self.pca_positions[r][0], self.pca_positions[r][1], color='#db5c5c', s=80, zorder=5)
                    if label_nodes: ax.annotate('W'+r.split('-')[-1], (self.pca_positions[r][0]+0.2, self.pca_positions[r][1]-0.25), fontsize=17, zorder=6)
                
                    
        if label_nodes:
            for node in self.conserved_nodes:
                ax.annotate(str(_hf.amino_d[node.split('-')[0]])+str(int(node.split('-')[1])), (self.pca_positions[node][0]+0.2, self.pca_positions[node][1]-0.25), fontsize=17, zorder=6)
        plt.savefig(self.plot_folder+'conserved_'+str(self.max_water)+self.graph_type+'_lables_'+str(label_nodes)+'_graph.png')
    
    
    
    
    
    
    def plot_unique_grpahs(self):
        pass
    
    def plot_unique_linear_lenght(self):
        pass