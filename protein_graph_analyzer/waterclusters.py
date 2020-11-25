import helperfunctions as _hf
# import numpy as np

class WaterClusters():
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

    def fit_parameters(self):
        pass

    def evaluate_parameters(self):
        #write to file
        pass

    def calculate_cluster_centers(self):
        pass

    def plot_clusters(self):
        pass

    def draw_cluster_centers(self):
        pass

    def calculate_water_clusters(self):
        self.fit_parameters()
        self.evaluate_parameters()
        self.calculate_cluster_centers()
        self.plot_clusters()
        self.draw_cluster_centers()
