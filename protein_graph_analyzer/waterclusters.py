from . import helperfunctions as _hf
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA
from sklearn import metrics
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt

from .proteingraphanalyser import ProteinGraphAnalyser #TODO: make independnet form this module

class WaterClusters(ProteinGraphAnalyser):
    def __init__(self, pdb_root_folder, target_folder='', reference_pdb='', sequance_identity_threshold=0.75):
        ProteinGraphAnalyser.__init__(self, pdb_root_folder=pdb_root_folder, target_folder=target_folder, reference_pdb=reference_pdb)
        waters = 0
        for file in self.file_list:
            waters += len(_hf.water_in_pdb(self.pdb_root_folder+file))
        if waters <= len(self.file_list)*2: #check this number or fiugre out somethinf
            print('There are not enough waters to cluster in the PDB files')
            return
        else:
            ProteinGraphAnalyser.align_structures(self, sequance_identity_threshold=sequance_identity_threshold)
            self.superimposed_files = _hf.get_files(self.target_folder, '_superimposed.pdb')
            self.water_coordinates = self._get_water_coordinates()

    def fit_parameters(self):
        neigh = NearestNeighbors(n_neighbors=5)
        nbrs = neigh.fit(self.water_coordinates[:,0:3])
        distances, indices = nbrs.kneighbors(self.water_coordinates[:,0:3])

        fig, ax = _hf.create_plot(figsize=(7,5),
                                  title='Optimal eps value',
                                  xlabel='Data points',
                                  ylabel='kNN distance')
        distances = np.sort(distances, axis=0)
        distances = distances[:,1]
        plt.plot(distances, color='#29335C')
        plt.yticks(np.arange(min(distances), max(distances)+1, 1.0))
        ax.grid(axis='y')
        plt.savefig(self.plot_folder+'kNN_distance_evaluation.png')
        plt.close()

        u = np.unique(distances)
        slope = []
        for i in range(len(u)-1):
            s = u[i+1] - u[i]
            slope.append(s)
        slope = np.array(slope)
        print(np.argmax(slope))
        print(u[np.argmax(slope)])
        print(distances[np.argmax(slope)+1])


#         fig.savefig(folder+'optimal_eps.png')

    def evaluate_parameters(self):
        #maybe set min_samples with hte conservation threshold
        EPS=1.5
        db = DBSCAN(eps=EPS, min_samples=len(self.superimposed_files)).fit(self.water_coordinates[:,0:3])
        core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True
        self.labels = db.labels_

        # Number of clusters in labels, ignoring noise if present.
        self.n_clusters_ = len(set(self.labels)) - (1 if -1 in self.labels else 0)
        n_noise_ = list(self.labels).count(-1)

        print('Estimated number of clusters: %d' % self.n_clusters_)
        print('Estimated number of noise points: %d' % n_noise_)
        print("Silhouette Coefficient: %0.3f"
              % metrics.silhouette_score(self.water_coordinates, self.labels))

        XY = self.water_coordinates[:,0:2]
        pca = PCA(n_components=1)
        xy = pca.fit_transform(XY)

        fig, ax = _hf.create_plot(title='Projection of the '+str(self.n_clusters_)+' water clusters and '+str(n_noise_)+' outlier points',
                                      xlabel='PCA projected xy plane',
                                      ylabel='Z coordinates')

        colormap = plt.cm.get_cmap('tab20', len(set(self.labels)))
        mycolors = colormap(np.linspace(0, 1, len(set(self.labels))))

        for x, y, l in zip(xy, self.water_coordinates[:,2], self.labels):
            if l == -1:
                ax.scatter(x, y, color='black', s=13)
            else:
                ax.scatter(x, y, color=mycolors[l], s=13)
        plt.savefig(self.plot_folder+'water_clusters.png')
        plt.close()



    def calculate_cluster_centers(self):
#         append water center coordinates to reference coordinates with water cluser number
        self.clusters = {}
        for i, lab in enumerate(self.labels):
            if lab in self.clusters.keys():
                self.clusters[lab].append(self.water_coordinates[i])
            else:
                self.clusters[lab] = [self.water_coordinates[i]]

        self.cluster_centers = {}
        size = []
        for key, value in self.clusters.items():
            coords = np.array(value)[:,:3]
            if key != -1:
                self.cluster_centers[key] = np.mean(coords, axis=0)
                size.append(len(value))

        for key, value in self.cluster_centers.items():
            res = 'w-'+str(key+1)
            self.reference_coordinates.update( {res:value} )

    def plot_clusters(self):
        pass

    def write_cluster_center_coordinates(self):
        np.savetxt(self.target_folder+'DBSCAN_water_cc_coordinates.txt', np.array(list(self.cluster_centers.values())))
        with open(self.target_folder+'DBSCAN_water_cc_coordinates.xyz', 'a') as file:
            file.write(str(len(self.cluster_centers.values()))+ ' \n')
            for i in self.cluster_centers.values():
                file.write('O '+ str(i[0]) + ' ' + str(i[1]) + ' ' + str(i[2]) + ' \n')

    def draw_clusters_centers_chimera(self):
        with open(self.target_folder+'DBSCAN_water_cc_chimera.bild', 'w') as f:
            f.write('.transparency 0.2\n')

            for cc in self.cluster_centers.values():
                f.write('.color  cornflower blue\n')
                f.write('.sphere '+ str(cc[0]) + ' '
                        + str(cc[1]) + ' '
                        + str(cc[2])+' '
                        + str(1)+' \n')


    def calculate_water_clusters(self):
#         self.fit_parameters()
        self.evaluate_parameters()
        self.calculate_cluster_centers()
#         self.plot_clusters()
#         self.draw_cluster_centers()

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
        ax.scatter(xy, self.water_coordinates[:,2], s=18, c='darkblue')
        plt.savefig(self.plot_folder+'all_water_projection.png')
        plt.close()

#         if file_name:
#             fig.savefig(file_name)

    def _get_water_coordinates(self):
        water_coordinates = []
        for file in self.superimposed_files:
            water_coord = _hf.water_coordinates(self.target_folder+file)
            water_coordinates.append(water_coord)
            _file_name = file.split('.pdb')[0]
            np.savetxt(self.target_folder+_file_name+'_water_coordinates.txt', water_coord)
        return _hf.concatenate_arrays(water_coordinates)
