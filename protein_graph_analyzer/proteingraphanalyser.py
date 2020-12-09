import networkx as nx
from mdhbond import WireAnalysis
import helperfunctions as _hf


class ProteinGraphAnalyser():
    def __init__(self, pdb_root_folder, target_folder=''):
        self.pdb_root_folder = pdb_root_folder+'/'
        if target_folder == '':
            self.target_folder = pdb_root_folder+'/'
        else: self.target_folder = target_folder+'/'

        self.file_list = _hf.get_pdb_files(self.pdb_root_folder)
            
    def align_structures(self, sequance_similarity_threshold=0.75, reference_pdb=''):
        if reference_pdb == '': reference_pdb = self.file_list[0]
        print('Reference strucure: ', reference_pdb)
        
        for pdb_move in self.file_list[1:]:
            ref_aligned, move_aligned = _hf.align_sequence(self.pdb_root_folder+reference_pdb,
                                                       self.pdb_root_folder+pdb_move,
                                                       threshold=sequance_similarity_threshold)
            if (ref_aligned is not None) and (move_aligned is not None):
                _hf.superimpose_aligned_atoms(ref_aligned, self.pdb_root_folder+reference_pdb,
                                          move_aligned, self.pdb_root_folder+pdb_move,
                                          file_name= self.target_folder+pdb_move)
    def number_of_waters_per_structure(self):
        for file in self.file_list:
            waters = _hf.water_in_pdb(self.pdb_root_folder+file)
            number_of_waters = len(waters)
            print(number_of_waters)
    
    def calculate_graphs(self, graph_type='water_wire', selection='protein', max_water=3, write_to_file=True):
        self.graphs = []
#         try: graph_type in ['water_wire', 'hbond']
        if graph_type == 'water_wire':
            for file in self.file_list:
                if file.endswith('.pdb'):
                    pdb_file = self.pdb_root_folder+file
                    print(pdb_file)
                    print(selection)
                    wba = WireAnalysis('protein',
                                       pdb_file,
                                       residuewise=True,
                                       check_angle=False,
                                       add_donors_without_hydrogen=True)
                    wba.set_water_wires(max_water=max_water)
                    wba.compute_average_water_per_wire()
                    g = wba.filtered_graph
                    self.graphs.append(g)
                    if write_to_file: nx.write_gpickle(g, self.target_folder+file.split('.pdb')[0]+'_graphs.pickle')

                else:
                    print('Pleases provide a list of pdb files for analysis')
                    return
                    #later support psf and dcd as well
#             if psf and dcd:
#                 pass
            
        elif graph_type == 'hbond':
#             only for pdb
            for file in self.file_list:
                pass
        else:
            print("graph_type has to be 'water_wire' or 'hbond' ")
            return
    
    def get_clusters(self):
        pass
    
    def plot_clusters(self):
        self.get_clusters()
    
    def get_linear_lenght(self):
        pass
    
    def plot_linear_lenght(self):
        self.get_linear_lenght()