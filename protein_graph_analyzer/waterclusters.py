import helperfunctions as _hf
# import numpy as np

class WaterClusters():
  def __init__(self, pdb_root_folder, target_folder=''):
    print(pdb_root_folder)
    self.pdb_root_folder = pdb_root_folder+'/'
    if target_folder == '':
      self.target_folder = pdb_root_folder
    else: self.target_folder = target_folder+'/'

    self.file_list = _hf.get_pdb_files(pdb_root_folder)
#     _hf._superimpose(pdb_root_folder, target_folder)

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
