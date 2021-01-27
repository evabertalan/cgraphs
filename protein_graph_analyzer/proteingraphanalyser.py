from . import helperfunctions as _hf
import shutil
import networkx as nx
import numpy as np
from . import mdhbond as mdh
import matplotlib.pyplot as plt


class ProteinGraphAnalyser():
    def __init__(self, pdb_root_folder, target_folder='', reference_pdb='', type_option='pdb'):
        #here set refernce file form the modal
        self.type_option = type_option
        self.pdb_root_folder = pdb_root_folder+'/'
        if target_folder == '':
            self.target_folder = _hf.create_directory(pdb_root_folder+'/workfolder')+'/'
        else: self.target_folder = target_folder+'/'
        self.plot_folder = _hf.create_directory(self.target_folder+'/plots/')
        self.max_water = ''

        if self.type_option == 'pdb':
            self.file_list = _hf.get_pdb_files(self.pdb_root_folder)
            shutil.copy(reference_pdb, self.target_folder+reference_pdb.split('/')[-1].split('.pdb')[0]+'_ref.pdb')
            self.reference_pdb = self.target_folder+_hf.get_files(self.target_folder, '_ref.pdb')[0]
            self.get_reference_coordinates()
            self._load_structures()

        elif self.type_option == 'dcd':
            dcd = sorted([folder+file for file in os.listdir(folder) if re.match(r'.*n\d{2}.dcd$', file) ])
            psf = [folder+file for file in os.listdir(folder) if re.match(r'read_protein.*.psf$', file) ][0]
            reference_pdb = [folder+file for file in os.listdir(folder) if re.match(r'read_protein.*.pdb$', file) ][0]
            ref = mda.Universe(pdb) #WHAT

        else: raise ValueError('Given type_option should be "pdb" or "dcd"')

    def _load_structures(self):
        self.graph_coord_objects = {}
        for file in self.file_list:
            print(file)
            print(len(_hf.water_in_pdb(self.pdb_root_folder+file)))
            structure = _hf.load_pdb_structure(self.pdb_root_folder+file)
            self.graph_coord_objects.update( {file.split('/')[-1].split('.pdb')[0]: {'structure': structure} } )

    def _load_exisitng_graphs(self):
        pass

#     def set_reference_file(self, reference_pdb=''):
#         #         print('1) Please select a reference file in case of membrane protein a correspongin OPM oriented strucurre is recommended. 2) use one of the files as refenecre form the list')

# #         if isMembraneProtein:
# #             print('Could not find OPM representation')
# #         self.reference_pdb =  self.target_folder
#         self.reference_pdb = reference_pdb

    def get_reference_coordinates(self, save=True):
        self.reference_coordinates = {}
        structure = _hf.load_pdb_structure(self.reference_pdb)
        for i in range(1, len(list(structure[0].get_residues()))):
            res_name = list(structure[0].get_residues())[i-1].get_resname()
            res_id = list(structure[0].get_residues())[i-1].get_id()[1]
            chain = list(structure[0].get_residues())[i-1].get_parent()
            if res_name in _hf.amino_d.keys():
#                 if res_name == 'HSD' or res_name == 'HSE': res_name='HIS'
                res = res_name+'-'+str(res_id)
                coord = list(structure[0].get_residues())[i-1]['CA'].get_coord()
                self.reference_coordinates.update( {res:coord} )
            elif res_name == 'HOH':
                res = res_name+'-'+str(res_id)
                coord = _hf.get_water_coordinates(chain, res_id)

                self.reference_coordinates.update( {res:coord} )

#                 print(res)

#                 chain = structure[0][node.split('-')[0]]
# #                 print(n)
# #                 print(chain)
#                 res_id = n.split('-')[-1]
#                 if n.split('-')[0] == 'HOH':
#                     #quickfix:
#                     if int(res_id) > 10000:
#                         res_id = int(res_id) - 10000
#                     coords = _hf.get_water_coordinates(chain, res_id)



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

    def calculate_graphs(self, graph_type='water_wire', selection='protein', max_water=3):
        self.graph_type = graph_type
        if self.type_option == 'pdb':
            try:
                self.graph_type in ['water_wire', 'hbond']
            except ValueError:
                raise ValueError('Given graph_type has to be "water_wire" or "hbond"')
            #TODO correct this logic
            self.file_list = [f for f in _hf.get_files(self.target_folder, '_superimposed.pdb')]
            if self.graph_type == 'water_wire':
                self.max_water = max_water
                for file in self.file_list:
                    pdb_file = self.target_folder+file
                    wba = mdh.WireAnalysis(selection,
                                       pdb_file,
                                       residuewise=True,
                                       check_angle=False,
                                       add_donors_without_hydrogen=True)
                    wba.set_water_wires(max_water=max_water)
                    wba.compute_average_water_per_wire()
                    g = wba.filtered_graph
                    nx.write_gpickle(g, self.target_folder+file.split('.pdb')[0]+self.graph_type+'_graphs.pickle')
                    self.graph_coord_objects[file.split('/')[-1].split('_superimposed.pdb')[0]].update( {'graph': g} )

            elif self.graph_type == 'hbond':
                for file in self.file_list:
                    pdb_file = self.target_folder+file
                    hba = mdh.HbondAnalysis(selection,
                                        pdb_file,
                                        residuewise=True,
                                        check_angle=False,
                                        add_donors_without_hydrogen=True,
                                        additional_donors=['N'],
                                        additional_acceptors=['O'])
                    hba.set_hbonds_in_selection(exclude_backbone_backbone=True)
                    hba.set_hbonds_in_selection_and_water_around(max_water)
                    g = hba.filtered_graph
                    nx.write_gpickle(g, self.target_folder+file.split('.pdb')[0]+'_'+self.graph_type+'_graphs.pickle')
                    self.graph_coord_objects[file.split('/')[-1].split('_superimposed.pdb')[0]].update( {'graph': g} )

        elif self.type_option == 'dcd' and self.graph_type == 'water_wire':
            for name, files in self.file_list.items():
                wba =  mdh.WireAnalysis(selection,
                                    files['psf'],
                                    files['dcd'],
                                    residuewise=True,
                                    check_angle=False,
                                    add_donors_without_hydrogen=True)
                wba.set_water_wires(water_in_convex_hull=max_water, max_water=max_water)
                wba.compute_average_water_per_wire()
                g = wba.filtered_graph
                wba.dump_to_file(self.target_folder+name+'_'+self.graph_type+'_wba_object.pickle')

        else: raise ValueError('For dcd analysis only graph_type="water_wire" is supported.')


    def _get_node_positions(self, objects):
        node_pos = {}
        for node in objects['graph'].nodes:
            n = _hf.get_node_name(node)
            if n not in self.reference_coordinates.keys():
                chain = objects['structure'][0][node.split('-')[0]]
                res_id = n.split('-')[-1]
                if n.split('-')[0] == 'HOH': coords = _hf.get_water_coordinates(chain, res_id)
                else: coords = chain[int(res_id)]['CA'].get_coord()
                node_pos.update( {n:list(coords)} )
            else: node_pos.update( {n:self.reference_coordinates[n]} )
        return _hf.calculate_pca_positions(node_pos)

    def plot_graphs(self, label_nodes=True, label_edges=True, xlabel='PCA projected xy plane', ylabel='Z coordinates'):
        for name, objects in self.graph_coord_objects.items():
            if 'graph' in objects.keys():
                print(name)
                fig, ax = _hf.create_plot(title=self.graph_type+' graph of '+name,
                                          xlabel=xlabel,
                                          ylabel=ylabel)
                node_pca_pos = self._get_node_positions(objects)
                node_pca_pos = _hf.check_projection_sign(node_pca_pos, self.pca_positions)

                for e in objects['graph'].edges:
                    e0 = _hf.get_node_name(e[0])
                    e1 = _hf.get_node_name(e[1])
                    edge_line = [node_pca_pos[e0], node_pca_pos[e1]]
                    x=[edge_line[0][0], edge_line[1][0]]
                    y=[edge_line[0][1], edge_line[1][1]]
                    ax.plot(x, y, color='gray', marker='o', linewidth=2, markersize=18, markerfacecolor='gray', markeredgecolor='gray')
                if self.graph_type == 'hbond':
                    for n, values in node_pca_pos.items():
                        if n.split('-')[0] == 'HOH':
                            ax.scatter(values[0],values[1], color='#db5c5c', s=120, zorder=5)

                if label_nodes:
                    for n, values in node_pca_pos.items():
                        if n.split('-')[0] == 'HOH': ax.annotate('W'+str(int(n.split('-')[1])), (values[0]+0.2, values[1]-0.25), fontsize=12)
                        else: ax.annotate(str(_hf.amino_d[n.split('-')[0]])+str(int(n.split('-')[1])), (values[0]+0.2, values[1]-0.25), fontsize=12)
                plt.tight_layout()
                is_label = '_labeled' if label_nodes else ''
                plt.savefig(self.plot_folder+name+'_'+str(self.max_water)+self.graph_type+is_label+'_graph.png')
                plt.close()


    def get_clusters(self):
        pass

    def plot_clusters(self):
        self.get_clusters()

    def get_linear_lenght(self, objects):
        connected_components = _hf.get_connected_components(objects['graph'])
        protein_chain = list(objects['structure'][0])[0]
        return _hf.calculate_connected_compontents_coordinates(connected_components, protein_chain)

    def plot_linear_lenghts(self):
        for name, objects in self.graph_coord_objects.items():
            if 'graph' in objects.keys():
                connected_components_coordinates = self.get_linear_lenght(objects)

                fig, ax = _hf.create_plot(figsize=(9,16),
                                        title=self.graph_type+' chains along the Z-axis of '+name,
                                        xlabel='# of nodes in the chain',
                                        ylabel='Z-axis coordinate')

                for i, g in enumerate(connected_components_coordinates):
                    for j in range(len(g)):
                        if connected_components_coordinates[i][j][0].split('-')[1] == 'HOH': color = '#db5c5c'
                        else: color = 'dimgray'
                        z_coords = connected_components_coordinates[i][j][1][2]
                        ax.scatter(i, z_coords, color=color, s=140)

                ax.set_xticks(np.arange(len(connected_components_coordinates)))
                ax.set_xticklabels([len(c) for c in connected_components_coordinates])

                plt.tight_layout()
                plt.savefig(self.plot_folder+name+'_'+str(self.max_water)+self.graph_type+'_linear_length.png')
                plt.close()
