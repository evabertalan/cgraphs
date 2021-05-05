from . import helperfunctions as _hf
import shutil
import numpy as np

from .proteingraphanalyser import ProteinGraphAnalyser
from .conservedgraph import ConservedGraph
import matplotlib.pyplot as plt



class CompareTwo(ProteinGraphAnalyser):
    def __init__(self, pdb1, pdb2, target_folder=''):
        self.pdb1_code = _hf.retrieve_pdb_code(pdb1, '.pdb')
        self.pdb2_code = _hf.retrieve_pdb_code(pdb2, '.pdb')
        self.workfolder = _hf.create_directory(target_folder+'/workfolder/compare_'+self.pdb1_code+'_'+self.pdb2_code)+'/'
        shutil.copy(pdb1, self.workfolder+self.pdb1_code+'.pdb')

        shutil.copy(pdb2, self.workfolder+self.pdb2_code+'.pdb')

        ProteinGraphAnalyser.__init__(self, pdb_root_folder=self.workfolder, target_folder=target_folder, reference_pdb=pdb1)
        self.logger.info('COMPARE STRUCTURES '+ self.pdb1_code + ' WITH ' + self.pdb2_code)
        ProteinGraphAnalyser.align_structures(self)
        self.pca_positions = _hf.calculate_pca_positions(self.reference_coordinates)

    def plot_graph_comparison(self, color1='blue', color2='green', label_nodes=True, label_edges=True, xlabel='PCA projected xy plane', ylabel='Z coordinates ($\AA$)',occupancy=None):
        if len(self.graph_coord_objects.items()) > 2: self.logger.warning('More than 2 structures are selected. Graph comparison is possible for two structures.')
        else:
            ConservedGraph.get_conserved_graph(self)
            plot_name = 'H-bond' if self.graph_type == 'hbond' else 'Water wire'
