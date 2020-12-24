import helperfunctions as _hf
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA

from proteingraphanalyser import ProteinGraphAnalyser

class WaterClusters(ProteinGraphAnalyser):
    def __init__(self, pdb_root_folder, target_folder=''):
        ProteinGraphAnalyser.__init__(self, pdb_root_folder, target_folder='')
        ProteinGraphAnalyser.align_structures(self)
        self.superimposed_files = _hf.get_files(self.pdb_root_folder, '_superimposed.pdb')
        self.water_coordinates = self._get_water_coordinates()

    def fit_parameters(self):
        pass

    def evaluate_parameters(self):
        #write to file
        pass

    def calculate_cluster_centers(self):
#         append water center coordinates to reference coordinates with water cluser number
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
        
    def plot_waters_along_membrane_normal(self, file_name=''):
#         if self.membrnaeProtein:
        pass
    
    def plot_all_waters(self,  file_name=''):
        XY = self.water_coordinates[:,0:2]
        pca = PCA(n_components=1)
        xy = pca.fit_transform(XY)
        
        fig, ax = _hf.create_plot(title='Projection of all water molecules',
                                 xlabel='PCA projected xy plane',
                                 ylabel='Z coordinates')
        ax.scatter(xy, self.water_coordinates[:,2], s=11, c='darkblue')
        
        if file_name:
            fig.savefig(file_name)
    
    def _get_water_coordinates(self):
        water_coordinates = []
        for file in self.superimposed_files:
            water_coord = _hf.water_coordinates(self.pdb_root_folder+file)
            water_coordinates.append(water_coord)
            _file_name = file.split('.pdb')[0]
            np.savetxt(self.target_folder+_file_name+'_water_coordinates.txt', water_coord)
        return _hf.concatenate_arrays(water_coordinates)
