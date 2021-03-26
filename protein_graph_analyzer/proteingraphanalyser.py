from . import helperfunctions as _hf
import shutil
import pickle
import networkx as nx
import numpy as np
from . import mdhbond as mdh
import matplotlib.pyplot as plt


class ProteinGraphAnalyser():
    def __init__(self, pdb_root_folder='', target_folder='', reference_pdb='', type_option='pdb', psf_file='', dcd_files=[], sim_name=''):
        #here set refernce file form the modal
        self.type_option = type_option
        self.pdb_root_folder = pdb_root_folder+'/'
        if target_folder == '':
            self.workfolder = _hf.create_directory(pdb_root_folder+'/workfolder')+'/'
        else: self.workfolder = _hf.create_directory(target_folder+'/workfolder')+'/'
        self.graph_object_folder = _hf.create_directory(self.workfolder+'/.graph_objects/')
        self.helper_files_folder = _hf.create_directory(self.workfolder+'/.helper_files/')
        self.plot_folder = _hf.create_directory(self.workfolder+'/plots/')
        self.max_water = 0
        self.graph_coord_objects = {}
        self.logger = _hf.create_logger(self.helper_files_folder)

        if self.type_option == 'pdb':
            self.logger.debug('Analysis for PDB crystal structures')
            self.file_list = _hf.get_pdb_files(self.pdb_root_folder)
            shutil.copy(reference_pdb, self.helper_files_folder+reference_pdb.split('/')[-1].split('.pdb')[0]+'_ref.pdb')
            self.reference_pdb = self.helper_files_folder+_hf.get_files(self.helper_files_folder, '_ref.pdb')[0]
            self.get_reference_coordinates(self.reference_pdb)
            self._load_structures()

        elif self.type_option == 'dcd' and len(dcd_files) and len(psf_file):
            self.logger.debug('Analysis for MD trajectories')

            # assert len(psf_file) == len(dcd_files) == len(sim_names) #later use this
            # for i in range(len(psf_files)):
            self.graph_coord_objects.update( { sim_name: {'psf': psf_file, 'dcd': dcd_files} } )
            # psf = [folder+file for file in os.listdir(folder) if re.match(r'read_protein.*.psf$', file) ][0]
            # reference_pdb = [folder+file for file in os.listdir(folder) if re.match(r'read_protein.*.pdb$', file) ][0]
            # ref = mda.Universe(pdb) #WHAT

        # else: raise ValueError('Given type_option should be "pdb" or "dcd"')

    def _load_structures(self):
        self.logger.info('Loading PDB crystal structures')
        for file in self.file_list:
            self.logger.debug('Loading structure: ', file)
            self.logger.debug('Number of water molecules in '+file+' is: '+str(len(_hf.water_in_pdb(self.pdb_root_folder+file))))
            structure = _hf.load_pdb_structure(self.pdb_root_folder+file)
            pdb_name = file.split('/')[-1].split('.pdb')[0]
            self.graph_coord_objects.update( { pdb_name: {'structure': structure} } )
            _hf.create_directory(self.plot_folder+pdb_name+'/')
            _hf.create_directory(self.plot_folder+pdb_name+'/hbond_graphs/')
            _hf.create_directory(self.plot_folder+pdb_name+'/water_wires/')
            self.logger.debug('Plot folders for each pdb structures are created')

    # def _load_trajectories(self, psf, dcd, sim_name):
    #     self.graph_coord_objects.update( { sim_name: {'psf': psf, 'dcd': dcd} } )

    def _load_exisitng_graphs(self, graph_files=None, reference_pdb='', graph_type='water_wire'):
        self.graph_type = graph_type
        self.logger.info('Loading graphs for simulations')
        self.logger.info('This step takes some time.')

        for i, graph_file in enumerate(graph_files):
            name = graph_file.split('/')[-1].split('_water_wire')[0]
            _hf.create_directory(self.plot_folder+name+'/')
            _hf.create_directory(self.plot_folder+name+'/water_wires/')
            self.logger.info('Loading '+name+'...')
            graph_coord_object = _hf.pickle_load_file(self.helper_files_folder+name+'_water_wire_graph_coord_objects.pickle')
            psf = graph_coord_object[name]['psf']
            dcd = graph_coord_object[name]['dcd']
            # wba = self.graph_coord_objects[name]['wba'] # try to use this later, but WBA needs to be saved in graph_coord_objects
            wba = mdh.WireAnalysis('protein', psf, dcd)
            wba.load_from_file(self.graph_object_folder+name+'_water_wire_wba_object.pickle', reload_universe=True)
            g = wba.filtered_graph
            mda = wba._mda_selection.select_atoms('name CA')
            self.graph_coord_objects[name] = {'psf': psf,  'dcd': dcd, 'wba': wba, 'graph': g, 'mda': mda}
            if i == 0:
                self.get_reference_coordinates(mda)
                self.pca_positions = _hf.calculate_pca_positions(self.reference_coordinates)
            self.logger.info(name+' loading completed')


    def get_reference_coordinates(self, reference, save=True):
        self.reference_coordinates = {}
        if self.type_option == 'pdb':
            structure = _hf.load_pdb_structure(reference)
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
        else:
            positions = reference.positions
            for i, resisdue in enumerate(reference):
                res_name, res_id = resisdue.resname, resisdue.resid
                res = res_name+'-'+str(res_id)
                self.reference_coordinates.update( {res:positions[i]} )

        if save: _hf.pickle_write_file(self.helper_files_folder+'reference_coordinate_positions.pickle', self.reference_coordinates)


    def align_structures(self, sequance_identity_threshold=0.75, isMembraneProtein=True):
        self.logger.debug('Reference strucure: ', self.reference_pdb)
        self.logger.info('Sequance identity threshold is set to: '+str(sequance_identity_threshold*100)+'%')
        self.superimposed_structures_folder = _hf.create_directory(self.workfolder+'/superimposed_structures/')


        for pdb_move in self.file_list:
            ref_aligned, move_aligned = _hf.align_sequence(self.logger, self.reference_pdb,
                                                       self.pdb_root_folder+pdb_move,
                                                       threshold=sequance_identity_threshold)
            if (ref_aligned is not None) and (move_aligned is not None):
                struct = _hf.superimpose_aligned_atoms(self.logger, ref_aligned, self.reference_pdb,
                                          move_aligned, self.pdb_root_folder+pdb_move,
                                          save_file_to= self.superimposed_structures_folder+pdb_move)
                self.graph_coord_objects.update( {pdb_move.split('/')[-1].split('.pdb')[0]: {'structure': struct} } )

    def number_of_waters_per_structure(self):
        for file in self.file_list:
            waters = _hf.water_in_pdb(self.pdb_root_folder+file)
            number_of_waters = len(waters)
            self.logger.info('Number of water molecules in '+file+' is: '+str(number_of_waters))


    def calculate_graphs(self, graph_type='water_wire', selection='protein', max_water=3, exclude_backbone_backbone=True, include_backbone_sidechain=False):
        self.graph_type = graph_type
        self.logger.info('Calculating graphs for '+self.graph_type+' analysis.')
        if self.type_option == 'pdb':
            try:
                self.graph_type in ['water_wire', 'hbond']
            except ValueError:
                raise ValueError('Given graph_type has to be "water_wire" or "hbond"')
            #TODO correct this logic
            self.file_list = [f for f in _hf.get_files(self.superimposed_structures_folder, '_superimposed.pdb')]
            if self.graph_type == 'water_wire':
                self.logger.info('Maximum number of water in water bridges is set to : '+str(max_water))
                self.max_water = max_water
                for file in self.file_list:
                    self.logger.debug('Calculating '+self.graph_type+' graph for: '+file)
                    pdb_file = self.superimposed_structures_folder+file
                    wba = mdh.WireAnalysis(selection,
                                       pdb_file,
                                       residuewise=True,
                                       check_angle=False,
                                       add_donors_without_hydrogen=True)
                    wba.set_water_wires(max_water=max_water)
                    wba.compute_average_water_per_wire()
                    g = wba.filtered_graph
                    nx.write_gpickle(g, self.graph_object_folder+file.split('.pdb')[0]+self.graph_type+'_graphs.pickle')
                    self.graph_coord_objects[file.split('/')[-1].split('_superimposed.pdb')[0]].update( {'graph': g} )

            elif self.graph_type == 'hbond':
                donors = []
                acceptors = []
                if include_backbone_sidechain:
                    self.logger.info('Including sidechain-backbone interactions')
                    donors.append('N')
                    acceptors.append('O')
                if not exclude_backbone_backbone:
                    self.logger.info('Including backbone-backbone interactions')
                for file in self.file_list:
                    self.logger.debug('Calculating '+self.graph_type+' graph for: '+file)
                    pdb_file = self.superimposed_structures_folder+file
                    hba = mdh.HbondAnalysis(selection,
                                        pdb_file,
                                        residuewise=True,
                                        check_angle=False,
                                        add_donors_without_hydrogen=True,
                                        additional_donors=donors,
                                        additional_acceptors=acceptors)
                    hba.set_hbonds_in_selection(exclude_backbone_backbone=exclude_backbone_backbone)
                    hba.set_hbonds_in_selection_and_water_around(max_water)
                    g = hba.filtered_graph
                    nx.write_gpickle(g, self.graph_object_folder+file.split('.pdb')[0]+'_'+self.graph_type+'_graphs.pickle')
                    self.graph_coord_objects[file.split('/')[-1].split('_superimposed.pdb')[0]].update( {'graph': g} )

        elif self.type_option == 'dcd' and self.graph_type == 'water_wire':
            self.logger.info('Maximum number of water in water bridges is set to : '+str(max_water))
            for name, files in self.graph_coord_objects.items():
                self.logger.info('Loading '+str(len(files['dcd']))+' trajectory files for '+name)
                self.logger.info('This step may take some time.')
                wba =  mdh.WireAnalysis(selection,
                                    files['psf'],
                                    files['dcd'],
                                    residuewise=True,
                                    check_angle=False,
                                    add_donors_without_hydrogen=True)
                wba.set_water_wires(water_in_convex_hull=max_water, max_water=max_water)
                wba.compute_average_water_per_wire()
                _hf.pickle_write_file(self.helper_files_folder+name+'_'+self.graph_type+'_graph_coord_objects.pickle', self.graph_coord_objects)
                g = wba.filtered_graph
                self.graph_coord_objects[name].update( {'wba': wba})
                self.graph_coord_objects[name].update( {'graph': g})
                wba.dump_to_file(self.graph_object_folder+name+'_'+self.graph_type+'_wba_object.pickle')
                nx.write_gpickle(g, self.helper_files_folder+name+'_'+self.graph_type+'_'+str(max_water)+'_water_graphs.pickle')
                self.logger.info('Graph object is saved as: '+self.graph_object_folder+name+'_'+self.graph_type+'_'+str(max_water)+'_water_graphs.pickle')

        else: raise ValueError('For dcd analysis only graph_type="water_wire" is supported.')


    def _get_node_positions(self, objects):
        node_pos = {}
        for node in objects['graph'].nodes:
            n = _hf.get_node_name(node)
            if n not in self.reference_coordinates.keys() or n.split('-')[0] == 'HOH':
                chain = objects['structure'][0][node.split('-')[0]]
                res_id = n.split('-')[-1]
                if n.split('-')[0] == 'HOH': coords = _hf.get_water_coordinates(chain, res_id)
                else: coords = chain[int(res_id)]['CA'].get_coord()
                node_pos.update( {n:list(coords)} )
            else: node_pos.update( {n:self.reference_coordinates[n]} )
        return _hf.calculate_pca_positions(node_pos)

    def plot_graphs(self, label_nodes=True, label_edges=True, xlabel='PCA projected xy plane', ylabel='Z coordinates', occupancy=0.1):
        for pdb_name, objects in self.graph_coord_objects.items():
            if 'graph' in objects.keys():
                self.logger.debug('Creating '+self.graph_type+' plot for: '+pdb_name)
                fig, ax = _hf.create_plot(title=self.graph_type+' graph of '+pdb_name,
                                          xlabel=xlabel,
                                          ylabel=ylabel)
                if self.type_option == 'dcd': node_pca_pos = self.pca_positions
                else:
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
                if self.graph_type == 'hbond':
                    plt.savefig(self.plot_folder+pdb_name+'/hbond_graphs/'+pdb_name+'_Hbond_graph'+is_label+'.png')
                    plt.savefig(self.plot_folder+pdb_name+'/hbond_graphs/'+pdb_name+'_Hbond_graph'+is_label+'.eps', format='eps')
                elif self.graph_type == 'water_wire':
                    waters = '_max_'+str(self.max_water)+'_water_bridges' if self.max_water > 0 else ''
                    plt.savefig(self.plot_folder+pdb_name+'/water_wires/'+pdb_name+waters+'_graph'+is_label+'.png')
                    plt.savefig(self.plot_folder+pdb_name+'/water_wires/'+pdb_name+waters+'_graph'+is_label+'.eps', format='eps')
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
        for pdb_name, objects in self.graph_coord_objects.items():
            if 'graph' in objects.keys():
                self.logger.debug('Creating '+self.graph_type+' linear length plot for: '+pdb_name)
                connected_components_coordinates = self.get_linear_lenght(objects)

                fig, ax = _hf.create_plot(figsize=(9,16),
                                        title=self.graph_type+' chains along the Z-axis of '+pdb_name,
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
                if self.graph_type == 'hbond':
                    plt.savefig(self.plot_folder+pdb_name+'/hbond_graphs/'+pdb_name+'_Hbond_linear_length.png')
                    plt.savefig(self.plot_folder+pdb_name+'/hbond_graphs/'+pdb_name+'_Hbond_linear_length.eps', format='eps')
                elif self.graph_type == 'water_wire':
                    waters = '_'+str(self.max_water)+'_water_bridges' if self.max_water > 0 else ''
                    plt.savefig(self.plot_folder+pdb_name+'/water_wires/'+pdb_name+waters+'_linear_length.png')
                    plt.savefig(self.plot_folder+pdb_name+'/water_wires/'+pdb_name+waters+'_linear_length.eps', format='eps')
                plt.close()
