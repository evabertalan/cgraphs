import helperfunctions as _hf
import numpy as np
from proteingraphanalyser import ProteinGraphAnalyser

class ConservedGraph(ProteinGraphAnalyser):
    def __init__(self, pdb_root_folder, target_folder=''):
        ProteinGraphAnalyser.__init__(self, pdb_root_folder, target_folder='')
    
    def get_conserved_graph(self, threshold=1):
        nodes = []
        edges = []
        if self.graph_type == 'water_wire':
            for graph in self.graphs:
                for node in graph.nodes:
                    node = node.split('-')[1]+'-'+ node.split('-')[2]
                    nodes.append(node)
                for edge in graph.edges:
                    e0 =  edge[0].split('-')[1]+'-'+ edge[0].split('-')[2]
                    e1 =  edge[1].split('-')[1]+'-'+ edge[1].split('-')[2]
                    if ([e1, e0]) in edges:
                        edges.append([e1, e0])
                    else: edges.append([e0, e1])
        
        elif self.graph_type == 'hbond':
            for graph in self.graphs:
                #here select conserved water by the superimposed ones 
#                 if useEstimatedConservedWaters:
#                     check whether water coordinates belong to any clustrer thatn give name
                for node in graph.nodes:
                    if node.split('-')[1] == 'HOH': node = node.split('-')[1]+'-w'
                    else: node = node.split('-')[1]+'-'+ node.split('-')[2]
                    nodes.append(node)
                for edge in graph.edges:
                    if edge[0].split('-')[1] == 'HOH': e0 = edge[0].split('-')[1]+'-w'
                    else: e0 =  edge[0].split('-')[1]+'-'+ edge[0].split('-')[2]

                    if edge[1].split('-')[1] == 'HOH': e1 = edge[1].split('-')[1]+'-w'
                    else: e1 =  edge[1].split('-')[1]+'-'+ edge[1].split('-')[2]

                    if ([e1, e0]) in edges:
                        edges.append([e1, e0])
                    else: edges.append([e0, e1])

        th = np.round(len(self.file_list) * threshold)
        u_nodes, c_nodes = np.unique(nodes, return_counts=True)
        self.conserved_nodes = u_nodes[np.where(c_nodes >= threshold)[0]]
        u_edges, c_edges = np.unique(edges, return_counts=True, axis=0)
        self.conserved_edges = u_edges[np.where(c_edges >= threshold)[0]]
        print(self.conserved_nodes)
        print(self.conserved_edges)
    
    def plot_conserved_graph(self):
        pass
    
    def plot_unique_grpahs(self):
        pass
    
    def plot_unique_linear_lenght(self):
        pass