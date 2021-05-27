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
            self.logger.warning('There are not enough waters to cluster in the PDB files. Cluster analysis can not be performed.')
            return
        else:
            self.water_cluster_folder = _hf.create_directory(self.workfolder+'/water_clusters/')
            ProteinGraphAnalyser.align_structures(self, sequance_identity_threshold=sequance_identity_threshold)
            self.superimposed_files = _hf.get_files(self.superimposed_structures_folder, '_superimposed.pdb')
            self.water_coordinates = self._get_water_coordinates()
            self.logger.info('There are '+str(len(self.water_coordinates))+' water molecules in the '+str(len(self.superimposed_files))+' superimposed files')

    def fit_parameters(self):
        self.logger.info('DBSCAN PARAMETER ANALYSIS')
        self.parameter_folder = _hf.create_directory(self.water_cluster_folder+'/parameter_analysis/')

        self.logger.info('Calculating nearest neighbor distances of datat points.')
        neigh = NearestNeighbors(n_neighbors=5)
        nbrs = neigh.fit(self.water_coordinates[:,0:3])
        distances, indices = nbrs.kneighbors(self.water_coordinates[:,0:3])

        fig, ax = _hf.create_plot(figsize=(12,8),
                                  title='Optimal eps value',
                                  xlabel='Data points',
                                  ylabel='kNN distance')
        distances = np.sort(distances, axis=0)
        distances = distances[:,1]
        plt.plot(distances, color='#29335C')
        plt.yticks(np.arange(min(distances), max(distances)+1, 1.0))
        ax.grid(axis='y')
        plt.tight_layout()
        plt.savefig(self.parameter_folder+'kNN_distances.png')
        plt.close()

        fx = np.unique(distances)
        second_deriv = []
        for i in range(1, len(fx)-1):
            ddx = (fx[i+1] - 2*fx[i] + fx[i-1]) / 1
            second_deriv.append(ddx)

        self.logger.debug('Calculating second derivative of nearest neighbor distances.')
        fig, ax = _hf.create_plot(figsize=(12,8),
                                  title='Second derivatives of kNN_distances',
                                  xlabel='Data points',
                                  ylabel='Second derivative')
        plt.tight_layout()
        plt.plot(second_deriv, color='#07a61c')
        plt.savefig(self.parameter_folder+'kNN_distances_second_derivative.png')
        plt.close()

        second_deriv = np.round(np.abs(np.array(second_deriv)), 2)
        indexes = [i for i in range(len(second_deriv)) if 0.02 <= second_deriv[i] <= 0.10]
        eps_candidates = fx[range(indexes[0], indexes[-1])]

        eps_cluster = []
        eps_noise = []
        silhouette_score = []
        for _eps in eps_candidates:
            db = DBSCAN(eps=_eps).fit(self.water_coordinates[:,0:3])
            core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
            core_samples_mask[db.core_sample_indices_] = True
            labels = db.labels_

            # Number of clusters in labels, ignoring noise if present.
            n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
            n_noise_ = list(labels).count(-1)

            eps_cluster.append(n_clusters_)
            eps_noise.append(n_noise_)
            if n_clusters_  < 2:
                self.logger.warning('With eps: '+str(np.round(_eps,2))+' only 1 water cluster was found. Result of water clutering might be unreliable.')
                silhouette_score.append(-1)
            else:
                silhouette_score.append(metrics.silhouette_score(self.water_coordinates, labels))

        self.logger.info('Saving evaluation of possible eps values.')
        np.savetxt(self.parameter_folder+'eps_evaluation.csv', np.c_[eps_candidates, eps_cluster, eps_noise, silhouette_score], delimiter=',', fmt='%.2f', header='eps_candidate, number_of_cluster, number_of_outliers, silhouette_score')

        self.logger.info('Best eps is where the number of clusters are the maximal with low number of outliers and the silhouette score is maximal.')
        # max_sc = np.argmax(silhouette_score)
        if np.max(silhouette_score) < 0:
            self.logger.warning('Silhouette score is negative. Result of water clutering might be unreliable.')
            self.logger.info('Water molecules in the selected structre set do not show reliable clustering pattern.')
        else:
            indexes = [i for i in range(len(silhouette_score)) if silhouette_score[i] >= np.max(silhouette_score)*0.98]
        self.logger.info('Maximum silhouette score is: '+str(np.max(silhouette_score)))
        _eps = np.array(eps_candidates)[indexes]
        max_clusters = np.array(eps_cluster)[indexes]
        min_outlier = np.array(eps_noise)[indexes]
        best_eps_index = np.argmin(max_clusters/min_outlier)
        self.logger.info('The best eps is: '+str(np.round(_eps[best_eps_index],2))+'. With '+ str(max_clusters[best_eps_index])+' clusters and '+str(min_outlier[best_eps_index])+' outliers.')


    def evaluate_parameters(self, eps=1.4, min_samples=None):
        self.logger.info('WATER CLUSTER ANALYSIS')
        min_samples = min_samples if min_samples else len(self.superimposed_files)
        self.logger.info('DBSCAN eps is set to: '+str(eps))
        self.logger.info('Minimum number of water molecules to be considered as a cluster: '+str(min_samples))
        if eps < 1.0 or eps > 1.6:
            self.logger.warning('For water cluster analysis the recommanded eps is 1.4. See parameter study in supporting information of: https://www.sciencedirect.com/science/article/pii/S1047847720302070')

        db = DBSCAN(eps=eps, min_samples=min_samples).fit(self.water_coordinates[:,0:3])
        core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True
        self.labels = db.labels_

        # Number of clusters in labels, ignoring noise if present.
        self.n_clusters_ = len(set(self.labels)) - (1 if -1 in self.labels else 0)
        n_noise_ = list(self.labels).count(-1)

        self.logger.info('Estimated number of clusters: %d' % self.n_clusters_)
        self.logger.info('Estimated number of noise points: %d' % n_noise_)
        self.logger.info("Silhouette Coefficient: %0.3f"
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
        plt.savefig(self.water_cluster_folder+'water_clusters.png')
        plt.savefig(self.water_cluster_folder+'water_clusters.eps', format='eps')
        self.logger.debug('Water cluster plot is saved to :'+self.water_cluster_folder+'water_clusters.png')
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
            res = 'X-w-'+str(key+1)
            self.reference_coordinates.update( {res:value} )

    def plot_clusters(self):
        pass

    def write_cluster_center_coordinates(self):
        np.savetxt(self.water_cluster_folder+'DBSCAN_water_cc_coordinates.txt', np.array(list(self.cluster_centers.values())))
        with open(self.water_cluster_folder+'DBSCAN_water_cc_coordinates.xyz', 'a') as file:
            file.write(str(len(self.cluster_centers.values()))+ ' \n')
            for i in self.cluster_centers.values():
                file.write('O '+ str(i[0]) + ' ' + str(i[1]) + ' ' + str(i[2]) + ' \n')

    def draw_clusters_centers_chimera(self):
        with open(self.water_cluster_folder+'DBSCAN_water_cc_chimera.bild', 'w') as f:
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
        plt.savefig(self.water_cluster_folder+'all_water_projection.png')
        plt.close()

#         if file_name:
#             fig.savefig(file_name)

    def _get_water_coordinates(self):
        water_coordinates = []
        for file in self.superimposed_files:
            water_coord = _hf.water_coordinates(self.superimposed_structures_folder+file)
            water_coordinates.append(water_coord)
            _file_name = file.split('.pdb')[0]
            self.logger.info('Number of water molecules in '+_file_name+' is: '+str(len(water_coord)))
            np.savetxt(self.water_cluster_folder+_file_name+'_water_coordinates.txt', water_coord)
        return _hf.concatenate_arrays(water_coordinates)
