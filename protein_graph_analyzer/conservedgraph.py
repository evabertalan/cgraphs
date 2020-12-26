import helperfunctions as _hf
import numpy as np
from proteingraphanalyser import ProteinGraphAnalyser
import matplotlib.pyplot as plt
import waterclusters as wc


class ConservedGraph(ProteinGraphAnalyser):
    def __init__(self, pdb_root_folder, target_folder='', reference_pdb=''):
        ProteinGraphAnalyser.__init__(self, pdb_root_folder, target_folder, reference_pdb)
        self.pca_positions = _hf.calculate_pca_positions(self.reference_coordinates)

    
    def get_conserved_graph(self, threshold=1):
        nodes = []
        edges = []
        if self.graph_type == 'water_wire':
            for objects in self.graph_coord_objects.values():
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
                graph = objects['graph']
#             for graph in self.graphs:
                #here select conserved water by the superimposed ones 
#                 if useEstimatedConservedWaters:
#                     check whether water coordinates belong to any clustrer thatn give name
                for node in graph.nodes:
                    nodes.append(_hf.get_node_name(node))
                for edge in graph.edges:
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
                ax.scatter(self.pca_positions[n][0], self.pca_positions[n][1], color='gray', s=100, zorder=10)
                if n.startswith('w'):
                    ax.scatter(self.pca_positions[n][0], self.pca_positions[n][1], color='red', s=100, zorder=10)
                
                    
        if label_nodes:
            for node in self.conserved_nodes:
                ax.annotate(str(_hf.amino_d[node.split('-')[0]])+str(int(node.split('-')[1])), (self.pca_positions[node][0]+0.2, self.pca_positions[node][1]-0.25), fontsize=17)
        plt.savefig(self.plot_folder+'Conserved_'+str(self.max_water)+self.graph_type+'_graph.png')
    
    
    
    
    
    
    def plot_unique_grpahs(self):
        pass
    
    def plot_unique_linear_lenght(self):
        pass