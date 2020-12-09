import helperfunctions as _hf
# import numpy as np

class WaterClusters(ProteinGraphAnalyser):
    def __init__(self, pdb_root_folder, target_folder=''):
        self.pdb_root_folder = pdb_root_folder+'/'
        if target_folder == '':
            self.target_folder = pdb_root_folder+'/'
        else: self.target_folder = target_folder+'/'
    

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
