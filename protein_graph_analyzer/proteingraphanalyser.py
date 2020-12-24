import networkx as nx
import shutil
from mdhbond import WireAnalysis
from mdhbond import HbondAnalysis
import helperfunctions as _hf


class ProteinGraphAnalyser():
    def __init__(self, pdb_root_folder, target_folder='', reference_pdb=''):
        #here set refernce file form the modal
        self.pdb_root_folder = pdb_root_folder+'/'
        if target_folder == '':
            self.target_folder = pdb_root_folder+'/'
        else: self.target_folder = target_folder+'/'

        self.file_list = _hf.get_pdb_files(self.pdb_root_folder)
        if reference_pdb == '': _ref = self.file_list[0]
        else: _ref = reference_pdb
        
        shutil.copy(_ref, self.target_folder+_ref.split('/')[-1].split('.pdb')[0]+'_ref.pdb')
        self.reference_pdb = self.target_folder+_hf.get_files(self.target_folder, '_ref.pdb')[0]
        self.get_reference_coordinates()
        self._load_structures()
    
    def _load_structures(self):
        self.graph_coord_objects = {}
        for file in self.file_list:
            structure = _hf.load_pdb_model(self.pdb_root_folder+file)
            self.graph_coord_objects.update( {file.split('/')[-1].split('.pdb')[0]: {'structure': structure} } )
        print(self.graph_coord_objects)

#     def set_reference_file(self, reference_pdb=''):
#         #         print('1) Please select a reference file in case of membrane protein a correspongin OPM oriented strucurre is recommended. 2) use one of the files as refenecre form the list')

# #         if isMembraneProtein:
# #             print('Could not find OPM representation')
# #         self.reference_pdb =  self.target_folder 
#         self.reference_pdb = reference_pdb
    
    def get_reference_coordinates(self, save=True):
        self.reference_coordinates = {}
        structure = _hf.load_pdb_model(self.reference_pdb)
        for i in range(1, len(list(structure[0].get_residues()))):
            res_name = list(structure[0].get_residues())[i-1].get_resname()
            res_id = list(structure[0].get_residues())[i-1].get_id()[1]
        
            if res_name in _hf.amino_d.keys():
                res = res_name+'-'+str(res_id)
                coord = list(structure[0].get_residues())[i-1]['CA'].get_coord()
                self.reference_coordinates.update( {res:coord} )
        
        if save: _hf.pickle_write_file(self.target_folder+'reference_coordinate_positions.pickle', self.reference_coordinates)
    
                      
    def align_structures(self, sequance_similarity_threshold=0.75, isMembraneProtein=True):
        print('Reference strucure: ', self.reference_pdb)
        
        for pdb_move in self.file_list:
            ref_aligned, move_aligned = _hf.align_sequence(self.reference_pdb,
                                                       self.pdb_root_folder+pdb_move,
                                                       threshold=sequance_similarity_threshold)
            if (ref_aligned is not None) and (move_aligned is not None):
                struct = _hf.superimpose_aligned_atoms(ref_aligned, self.reference_pdb,
                                          move_aligned, self.pdb_root_folder+pdb_move,
                                          file_name= self.target_folder+pdb_move)
                self.graph_coord_objects.update( {pdb_move.split('/')[-1].split('.pdb')[0]: {'structure': struct} } )
    
    def number_of_waters_per_structure(self):
        for file in self.file_list:
            waters = _hf.water_in_pdb(self.pdb_root_folder+file)
            number_of_waters = len(waters)
            print(number_of_waters)
    
    def calculate_graphs(self, graph_type='water_wire', selection='protein', max_water=3, write_to_file=True):
        self.graph_type = graph_type
        self.graphs = []
#         try: graph_type in ['water_wire', 'hbond']
        if self.graph_type == 'water_wire':
            for file in self.file_list:
                if file.endswith('.pdb'):
                    pdb_file = self.pdb_root_folder+file
                    wba = WireAnalysis(selection,
                                       pdb_file,
                                       residuewise=True,
                                       check_angle=False,
                                       add_donors_without_hydrogen=True)
                    wba.set_water_wires(max_water=max_water)
                    wba.compute_average_water_per_wire()
                    g = wba.filtered_graph
                    self.graphs.append(g) # TODO maybe this graph object is not needed and coord obj could be used everywhere
                    self.graph_coord_objects[file.split('/')[-1].split('.pdb')[0]].update( {'graph': g} )
                    if write_to_file: nx.write_gpickle(g, self.target_folder+file.split('.pdb')[0]+'_graphs.pickle')

                else:
                    print('Pleases provide a list of pdb files for analysis')
                    return
                    #later support psf and dcd as well
#             if psf and dcd:
#                 pass
            
        elif self.graph_type == 'hbond':
#             only for pdb
            for file in self.file_list:
        
                if file.endswith('.pdb'):
                    print(file)
                    pdb_file = self.pdb_root_folder+file
                    hba = HbondAnalysis(selection,
                                        pdb_file, 
                                        residuewise=True, 
                                        check_angle=False,
                                        additional_donors=['N'], 
                                        additional_acceptors=['O'], 
                                        add_donors_without_hydrogen=True)
                    hba.set_hbonds_in_selection(exclude_backbone_backbone=True)
                    hba.set_hbonds_in_selection_and_water_around(max_water)
                    g = hba.filtered_graph
                    self.graphs.append(g)
                else:
                    print('For H-bond analysis only pdb files are supported')
                    return
        else:
            print("graph_type has to be 'water_wire' or 'hbond' ")
            return
    
    def get_node_positions(self, objects):
        node_pos = {}
        for n in objects['graph'].nodes:
            n = n.split('-')[1]+'-'+ n.split('-')[2]
            if n in self.reference_coordinates.keys(): node_pos.update( {n:self.reference_coordinates[n]} )
            else: 
                chains = list(objects['structure'][0].get_chains())
                chain_id = chains[0].get_id()
                res_id = n.split('-')[-1]
                node_pos.update( {n:list(objects['structure'][0][chain_id][int(res_id)]['CA'].get_coord())} )
        return _hf.calculate_pca_positions(node_pos)
    
    def plot_graphs(self, label_nodes=True, label_edges=True):
        for name, objects in self.graph_coord_objects.items():
            fig, ax = _hf.create_plot(title='H-bond graph of '+name,
                                      xlabel='PCA projected xy plane',
                                      ylabel='Z coordinates')
            node_pca_pos = self.get_node_positions(objects)
            for values in node_pca_pos.values():
                ax.scatter(values[0], values[1], s=90, c='gray')
            plt.sav


            
        
        
    def get_clusters(self):
        pass
    
    def plot_clusters(self):
        self.get_clusters()
    
    def get_linear_lenghts(self):
        pass
    
    def plot_linear_lenghts(self):
        self.get_linear_lenghts()
    
