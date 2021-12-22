from . import helperfunctions as _hf
import copy
import networkx as nx
import numpy as np
import MDAnalysis as _mda
from . import mdhbond as mdh
import matplotlib.pyplot as plt


class ProteinGraphAnalyser():
    def __init__(self, pdb_root_folder='', target_folder='', reference_pdb='', type_option='pdb', psf_files=[], dcd_files=[[]], sim_names=[], plot_parameters={}):
        #here set refernce file form the modal
        self.plot_parameters = {
            'edge_width': plot_parameters['edge_width'] if 'edge_width' in plot_parameters.keys() else 2,
            'node_label_size': plot_parameters['node_label_size'] if 'node_label_size' in plot_parameters.keys() else 12,
            'edge_label_size': plot_parameters['edge_label_size'] if 'edge_label_size' in plot_parameters.keys() else 10,
            'node_size': plot_parameters['node_size'] if 'node_size' in plot_parameters.keys() else 150,
            'graph_color': plot_parameters['graph_color'] if 'graph_color' in plot_parameters.keys() else 'gray',
            'water_node_color': plot_parameters['water_node_color'] if 'water_node_color' in plot_parameters.keys() else '#db5c5c',
            'difference_graph_color': plot_parameters['difference_graph_color'] if 'difference_graph_color' in plot_parameters.keys() else '#129fe6',
            'non_prot_color': plot_parameters['non_prot_color'] if 'non_prot_color' in plot_parameters.keys() else 'blue',
            'plot_title_fontsize': plot_parameters['plot_title_fontsize'] if 'plot_title_fontsize' in plot_parameters.keys() else 20,
            'plot_label_fontsize': plot_parameters['plot_label_fontsize'] if 'plot_label_fontsize' in plot_parameters.keys() else 36,
            'plot_tick_fontsize': plot_parameters['plot_tick_fontsize'] if 'plot_tick_fontsize' in plot_parameters.keys() else 33,
            'plot_resolution': plot_parameters['plot_resolution'] if 'plot_resolution' in plot_parameters.keys() else 400,
            'figsize': plot_parameters['figsize'] if 'figsize' in plot_parameters.keys() else (15, 16),
            'formats': plot_parameters['formats'] if 'formats' in plot_parameters.keys() else ['png'],
        }
        self.type_option = type_option
        self.pdb_root_folder = pdb_root_folder+'/'
        self.workfolder = _hf.create_directory(f'{pdb_root_folder}/workfolder/') if target_folder == '' else _hf.create_directory(f'{target_folder}/workfolder/')
        self.graph_object_folder = _hf.create_directory(f'{self.workfolder}/graph_objects/')
        self.helper_files_folder = _hf.create_directory(f'{self.workfolder}/.helper_files/')
        self.max_water = 0
        self.graph_coord_objects = {}
        self.logger = _hf.create_logger(self.helper_files_folder)

        if self.type_option == 'pdb':
            self.logger.debug('Analysis for PDB crystal structures')
            self.file_list = _hf.get_pdb_files(self.pdb_root_folder)
            self.reference_pdb = reference_pdb
            self.get_reference_coordinates(self.reference_pdb)
            self._load_structures()

        elif self.type_option == 'dcd' and len(dcd_files) and psf_files:
            self.reference_coordinates = {}
            self.logger.debug('Analysis for MD trajectories')
            if len(psf_files) == len(dcd_files) == len(sim_names):
                for i in range(len(psf_files)):
                    self.graph_coord_objects.update( { sim_names[i]: {'psf': psf_files[i], 'dcd': dcd_files[i]} } )
            else: self.logger.warning('Not equal number of parameters. Each simulation has to have a psf file, a name and a list of dcd files.')

        # else: raise ValueError('Given type_option should be "pdb" or "dcd"')

    def _load_structures(self):
        self.logger.info('Loading '+str(len(self.file_list))+' PDB crystal structures')
        for file in self.file_list:
            self.logger.debug('Loading structure: ', file)
            self.logger.debug(f'Number of water molecules in {file} is: {len(_hf.water_in_pdb(self.pdb_root_folder+file))}')
            structure = _hf.load_pdb_structure(self.pdb_root_folder+file)
            pdb_code = _hf.retrieve_pdb_code(file, '.pdb')
            self.graph_coord_objects.update( { pdb_code: {'structure': structure} } )
            self.logger.debug('Plot folders for each pdb structures are created')


    def _load_exisitng_graphs(self, graph_files=None, reference_pdb='', graph_type='water_wire', selection='protein'):
        self.graph_type = graph_type
        self.logger.info('Loading graphs for simulations')
        self.logger.info('This step takes some time.')
        self.reference_coordinates = {}

        for i, graph_file in enumerate(graph_files):
            _split = _hf.retrieve_pdb_code(graph_file, '_water_wire')
            name = _split[0:-2]
            self.max_water = int(_split[-1:])
            self.logger.info('Loading '+name+'...')
            graph_coord_object = _hf.pickle_load_file(f'{self.helper_files_folder}{name}_{self.max_water}_water_wires_coord_objects.pickle')
            psf = graph_coord_object[name]['psf']
            dcd = graph_coord_object[name]['dcd']

            self.selection = graph_coord_object[name]['selection'] if 'selection' in graph_coord_object[name].keys() else selection
            self.logger.info(f'Loading selection from the graph calculation. Selection: {self.selection}')

            wba = mdh.WireAnalysis(self.selection, psf, dcd)
            wba.load_from_file(graph_file, reload_universe=False)
            g = wba.filtered_graph
            u = _mda.Universe(psf, dcd)
            mda = u.select_atoms(str(self.selection))
            self.graph_coord_objects[name] = {'psf': psf,  'dcd': dcd, 'wba': wba, 'graph': g, 'mda': mda}
            self.add_reference_from_structure(mda, g)
            self.logger.info(f'{name} loading completed')
        self.pca_positions = _hf.calculate_pca_positions(self.reference_coordinates)

    def get_reference_coordinates(self, reference, save=True):
        self.reference_coordinates = {}
        if self.type_option == 'pdb':
            structure = _hf.load_pdb_structure(reference)
            if hasattr(self, 'selection') and self.selection!='protein':
                selection = structure.select_atoms(f'resname BWX or {self.selection} or {_hf.water_def}')
            else: selection = structure.select_atoms(f'resname BWX or (protein and name CA) or {_hf.water_def}')
            positions = selection.positions
            for i, resisdue in enumerate(selection):
                chain, res_name, res_id = resisdue.segid ,resisdue.resname, resisdue.resid
                res = chain+'-'+res_name+'-'+str(res_id)
                self.reference_coordinates.update( {res:positions[i]} )
        else:
            positions = reference.positions
            for i, resisdue in enumerate(reference):
                chain, res_name, res_id = resisdue.segid ,resisdue.resname, resisdue.resid
                res = chain+'-'+res_name+'-'+str(res_id)
                self.reference_coordinates.update( {res:positions[i]} )

        if save: _hf.pickle_write_file(self.helper_files_folder+'reference_coordinate_positions.pickle', self.reference_coordinates)


    def align_structures(self, sequance_identity_threshold=0.75, superimposition_threshold=5):
        self.logger.debug('Reference structure: ', self.reference_pdb)
        self.logger.info(f'Sequence identity threshold is set to: {sequance_identity_threshold*100}%')
        self.logger.info(f'Superimposition RMS threshold is set to: {superimposition_threshold}')
        _hf.delete_directory(self.workfolder+'/.superimposed_structures/')
        self.superimposed_structures_folder = _hf.create_directory(self.workfolder+'/.superimposed_structures/')

        for pdb_move in self.file_list:
            struct = None
            pdb_code = _hf.retrieve_pdb_code(pdb_move, '.pdb')
            ref_aligned, move_aligned = _hf.align_sequence(self.logger, self.reference_pdb,
                                                       self.pdb_root_folder+pdb_move,
                                                       threshold=sequance_identity_threshold)
            if (ref_aligned is not None) and (move_aligned is not None):
                struct = _hf.superimpose_aligned_atoms(self.logger, ref_aligned, self.reference_pdb,
                                          move_aligned, self.pdb_root_folder+pdb_move,
                                          save_file_to=self.superimposed_structures_folder+pdb_move, superimposition_threshold=superimposition_threshold)
                if struct is not None:
                    self.graph_coord_objects.update( {pdb_code: {'structure': struct, 'file': self.superimposed_structures_folder+pdb_code+'_superimposed.pdb'}} )
            if struct is None:
                self.graph_coord_objects.pop(pdb_code)

    def number_of_waters_per_structure(self):
        for file in self.file_list:
            waters = _hf.water_in_pdb(self.pdb_root_folder+file)
            self.logger.info(f'Number of water molecules in {file} is: {len(waters)}')

    def add_reference_from_structure(self, s, g):
        for i, resisdue in enumerate(s):
            chain, res_name, res_id = resisdue.segid ,resisdue.resname, resisdue.resid
            res = chain+'-'+res_name+'-'+str(res_id)
            if res not in self.reference_coordinates.keys() and res in g.nodes:
                self.reference_coordinates.update( {res: resisdue.position} )

    def calculate_graphs(self, graph_type='water_wire', selection='protein', max_water=3, exclude_backbone_backbone=True, include_backbone_sidechain=False, include_waters=True, distance=3.5, cut_angle=60., check_angle=False, additional_donors=[], additional_acceptors=[]):
        assert (type(additional_donors) is list and type(additional_acceptors) is list)
        if additional_donors or additional_acceptors:
            self.logger.info(f'List of additional donors: {additional_donors}\nList of additional acceptors: {additional_acceptors}')
        if selection!='protein': self.logger.info(f'Atom selection string: {selection}')
        self.graph_type = graph_type
        self.logger.info(f'Calculating graphs for {self.graph_type} analysis.')
        if self.type_option == 'pdb':
            self.selection = f'({selection}) or resname BWX'
            # self.get_reference_coordinates(self.reference_pdb)
            try:
                self.graph_type in ['water_wire', 'hbond']
            except ValueError:
                raise ValueError('Given graph_type has to be "water_wire" or "hbond"')
            self.logger.info(f'H-bond criteria cut off distance: {distance} A')
            if check_angle: self.logger.info(f'H-bond criteria cut off angle: {cut_angle} degree')
            if self.graph_type == 'water_wire':
                self.water_graphs_folder = _hf.create_directory(f'{self.graph_object_folder}{self.max_water}_water_wires/')
                self.logger.info(f'Maximum number of water in water bridges is set to: {max_water}')
                self.max_water = max_water
                for name, objects in self.graph_coord_objects.items():
                    pdb_file, pdb_code = objects['file'], name
                    if len(_hf.water_in_pdb(pdb_file)) == 0:
                        self.logger.warning(f'There are no water molecules in {pdb_code}. Water wire can not be calculated.')
                    else:
                        self.logger.debug(f'Calculating {self.graph_type} graph for: {pdb_code}')
                        if len(_hf.water_in_pdb(pdb_file)) == 0:
                            self.logger.warning(f'There are no water molecules in {pdb_code}. Water wire can not be calculated. Please use the H-bond network option.')
                        else:
                            wba = mdh.WireAnalysis(self.selection,
                                               pdb_file,
                                               residuewise=True,
                                               check_angle=check_angle,
                                               add_donors_without_hydrogen=not check_angle,
                                               additional_donors=additional_donors,
                                               additional_acceptors=additional_acceptors,
                                               distance=distance,
                                               cut_angle=cut_angle)
                            wba.set_water_wires(max_water=max_water)
                            wba.compute_average_water_per_wire()
                            g = wba.filtered_graph
                            nx.write_gpickle(g, self.water_graphs_folder+pdb_code+'_'+self.graph_type+'_graphs.pickle')
                            self.graph_coord_objects[pdb_code].update( {'graph': g} )
                            tmp = copy.copy(wba)
                            tmp._universe=None
                            self.graph_coord_objects[pdb_code].update( {'wba': tmp})
                            edge_info = _hf.edge_info(wba, g.edges)
                            _hf.json_write_file(self.helper_files_folder+pdb_code+'_'+self.graph_type+'_graph_edge_info.json', edge_info)
                            s = objects['structure'].select_atoms(self.selection)
                            self.add_reference_from_structure(s, g)


            elif self.graph_type == 'hbond':
                if include_backbone_sidechain:
                    self.logger.info('Including sidechain-backbone interactions')
                    additional_donors.append('N')
                    additional_acceptors.append('O')
                if not exclude_backbone_backbone:
                    self.logger.info('Including backbone-backbone interactions')
                for name, objects in self.graph_coord_objects.items():
                    pdb_file, pdb_code = objects['file'], name
                    self.logger.debug(f'Calculating {self.graph_type} graph for: {pdb_code}')
                    hba = mdh.HbondAnalysis(self.selection,
                                        pdb_file,
                                        residuewise=True,
                                        check_angle=check_angle,
                                        add_donors_without_hydrogen=not check_angle,
                                        additional_donors=additional_donors,
                                        additional_acceptors=additional_acceptors,
                                        distance=distance,
                                        cut_angle=cut_angle)
                    hba.set_hbonds_in_selection(exclude_backbone_backbone=exclude_backbone_backbone)
                    if len(_hf.water_in_pdb(pdb_file)) > 0 and include_waters: hba.set_hbonds_in_selection_and_water_around(max_water)
                    g = hba.filtered_graph
                    nx.write_gpickle(g, self.graph_object_folder+pdb_code+'_'+self.graph_type+'_graphs.pickle')
                    self.graph_coord_objects[pdb_code].update( {'graph': g} )
                    s = objects['structure'].select_atoms(f'resname BWX or {self.selection} or {_hf.water_def}')
                    self.add_reference_from_structure(s, g)

        elif self.type_option == 'dcd' and self.graph_type == 'water_wire':
            self.selection = selection
            self.max_water = max_water
            self.water_graphs_folder = _hf.create_directory(self.graph_object_folder+str(self.max_water)+'_water_wires/')
            self.logger.info(f'Maximum number of water in water bridges is set to : {max_water}')
            self.logger.info(f'H-bond criteria cut off distance: {distance} A')
            if check_angle: self.logger.info(f'H-bond criteria cut off angle: {cut_angle} degree')
            for i, (name, files) in enumerate(self.graph_coord_objects.items()):
                self.logger.info(f'Loading '+str(len(files['dcd']))+' trajectory files for '+name)
                self.logger.info('This step takes some time.')
                wba =  mdh.WireAnalysis(self.selection,
                                    files['psf'],
                                    files['dcd'],
                                    residuewise=True,
                                    check_angle=check_angle,
                                    add_donors_without_hydrogen=not check_angle,
                                    additional_donors=additional_donors,
                                    additional_acceptors=additional_acceptors,
                                    distance=distance,
                                    cut_angle=cut_angle)
                wba.set_water_wires(water_in_convex_hull=max_water, max_water=max_water)
                wba.compute_average_water_per_wire()
                self.graph_coord_objects[name].update( {'selection': self.selection})
                _hf.pickle_write_file(self.helper_files_folder+name+'_'+str(self.max_water)+'_water_wires_coord_objects.pickle', { name: self.graph_coord_objects[name] })
                g = wba.filtered_graph
                self.graph_coord_objects[name].update( {'graph': g})
                tmp = copy.copy(wba)
                tmp._universe=None
                self.graph_coord_objects[name].update( {'wba': tmp})
                u = _mda.Universe(files['psf'], files['dcd'])
                mda = u.select_atoms(self.selection)
                self.graph_coord_objects[name].update( {'mda': mda})
                wba_loc = self.water_graphs_folder+name+'_'+str(self.max_water)+'_water_wires_graph.pickle'
                wba.dump_to_file(wba_loc)
                nx.write_gpickle(g, self.helper_files_folder+name+'_'+self.graph_type+'_'+str(max_water)+'_water_nx_graphs.pickle')
                edge_info = _hf.edge_info(wba, g.edges)
                _hf.json_write_file(self.helper_files_folder+name+'_'+self.graph_type+'_'+str(max_water)+'_water_graph_edge_info.json', edge_info)
                self.add_reference_from_structure(mda, g)
                # if i == 0:
                #     self.get_reference_coordinates(mda)
                #     self.pca_positions = _hf.calculate_pca_positions(self.reference_coordinates)
                self.logger.info('Graph object is saved as: '+wba_loc)

        else: raise ValueError('For dcd analysis only graph_type="water_wire" is supported.')
        self.pca_positions = _hf.calculate_pca_positions(self.reference_coordinates)


    def _get_node_positions(self, objects, pca=True):
        node_pos = {}
        for node in objects['graph'].nodes:
            n = _hf.get_node_name(node)
            if n not in self.reference_coordinates.keys() or n.split('-')[1] in ['HOH', 'TIP3']:
                chain_id, res_name, res_id  = _hf.get_node_name_pats(n)
                if self.type_option == 'pdb':
                    chain = objects['structure'].select_atoms('segid '+ chain_id)
                    if res_name in ['HOH', 'TIP3']:
                        coords = _hf.get_water_coordinates(chain, res_id)
                    elif res_name in _hf.amino_d.keys():
                        coords = chain.select_atoms(f'resname BWX or protein and name CA and resid {res_id}').positions[0]
                    else: coords = chain.select_atoms(f'resname {res_name} and resid {res_id}').positions[0]
                elif self.type_option == 'dcd':
                    coords = objects['mda'].select_atoms('resid '+ res_id).positions[0]
                if coords is not None: node_pos.update({ n : list(coords) })
            else: node_pos.update({ n : self.reference_coordinates[n] })
        if pca: return _hf.calculate_pca_positions(node_pos)
        else: return node_pos

    def plot_graphs(self, label_nodes=True, label_edges=True, xlabel='PCA projected xy plane', ylabel='Z coordinates ($\AA$)', occupancy=None):
        if occupancy is not None or hasattr(self, 'occupancy'):
            if occupancy is not None: occupancy = occupancy
            else: occupancy = self.occupancy
        for name, objects in self.graph_coord_objects.items():
            if 'graph' in objects.keys():
                if occupancy:
                    wba = copy.deepcopy(objects['wba'])
                    wba.filter_occupancy(occupancy)
                    graph = wba.filtered_graph
                else: graph = objects['graph']
                self.logger.debug('Creating '+self.graph_type+' plot for: '+name)
                plot_name = 'H-bond' if self.graph_type == 'hbond' else 'Water wire'
                fig, ax = _hf.create_plot(title=f'{plot_name} graph of structure {name}\nSelection:{self.selection[1:-16]}',
                                          xlabel=xlabel,
                                          ylabel=ylabel, plot_parameters=self.plot_parameters)
                node_pca_pos = self._get_node_positions(objects)
                node_pca_pos = _hf.check_projection_sign(node_pca_pos, self.pca_positions)

                for e in graph.edges:
                    e0 = _hf.get_node_name(e[0])
                    e1 = _hf.get_node_name(e[1])
                    if e0 in node_pca_pos.keys() and e1 in node_pca_pos.keys():
                        edge_line = [node_pca_pos[e0], node_pca_pos[e1]]
                        x=[edge_line[0][0], edge_line[1][0]]
                        y=[edge_line[0][1], edge_line[1][1]]
                        ax.plot(x, y, color=self.plot_parameters['graph_color'], marker='o', linewidth=self.plot_parameters['edge_width'], markersize=self.plot_parameters['node_size']*0.01, markerfacecolor=self.plot_parameters['graph_color'], markeredgecolor=self.plot_parameters['graph_color'])
                        if label_edges and self.graph_type == 'water_wire':
                            waters, occ_per_wire, _ = _hf.get_edge_params(objects['wba'], graph.edges)
                            ax.annotate(np.round(waters[list(graph.edges).index(e)],1), (x[0] + (x[1]-x[0])/2, y[0] + (y[1]-y[0])/2), color='indianred',  fontsize=self.plot_parameters['edge_label_size'], weight='bold',)
                            ax.annotate(int(occ_per_wire[list(graph.edges).index(e)]*100), (x[0] + (x[1]-x[0])/2, y[0] + (y[1]-1.0-y[0])/2), color='green',  fontsize=self.plot_parameters['edge_label_size'])
                    # else:
                    #     self.logger.warning('Edge '+e0+'-'+e1+' is not in node positions. Can be an due to too many atoms in the PDB file.')

                for n, values in node_pca_pos.items():
                    if n.split('-')[1] in ['HOH', 'TIP3']:
                        ax.scatter(values[0],values[1], color=self.plot_parameters['water_node_color'], s=self.plot_parameters['node_size']*0.7, zorder=5)
                    elif n.split('-')[1] in _hf.amino_d.keys():
                        ax.scatter(values[0],values[1], color=self.plot_parameters['graph_color'], s=self.plot_parameters['node_size'], zorder=5)
                    else:
                        ax.scatter(values[0],values[1], color=self.plot_parameters['non_prot_color'], s=self.plot_parameters['node_size'], zorder=5)

                if label_nodes:
                    for n in graph.nodes:
                        n = _hf.get_node_name(n)
                        if n in node_pca_pos.keys():
                            values = node_pca_pos[n]
                            chain_id, res_name, res_id = _hf.get_node_name_pats(n)
                            if res_name in ['HOH', 'TIP3']: ax.annotate(f'W{res_id}', (values[0]+0.2, values[1]-0.25), fontsize=self.plot_parameters['node_label_size'])
                            elif res_name in _hf.amino_d.keys():
                                ax.annotate(f'{chain_id}-{_hf.amino_d[res_name]}{res_id}', (values[0]+0.2, values[1]-0.25), fontsize=self.plot_parameters['node_label_size'])
                            else: ax.annotate(f'{chain_id}-{res_name}{res_id}', (values[0]+0.2, values[1]-0.25), fontsize=self.plot_parameters['node_label_size'], color=self.plot_parameters['non_prot_color'])

                plt.tight_layout()
                is_label = '_labeled' if label_nodes else ''
                if self.graph_type == 'hbond':
                    plot_folder = _hf.create_directory(self.workfolder+'/H-bond_graphs/'+name+'/')
                    for form in self.plot_parameters['formats']:
                        plt.savefig(f'{plot_folder}{name}_H-bond_graph{is_label}.{form}', format=form, dpi=self.plot_parameters['plot_resolution'])
                    if is_label:
                        _hf.write_text_file(plot_folder+name+'_H-bond_graph_info.txt',
                            ['H-bond graph of '+name,
                            '\nSelection string: '+str(self.selection[0:-15]),
                            '\n',
                            '\nNumber of nodes in '+name+': '+str(len(graph.nodes)),
                            '\nNumber of edges in '+name+': '+str(len(graph.edges)),
                            '\n',
                            '\nList of nodes: '+str(graph.nodes),
                            '\n',
                            '\nList of edges: '+str(graph.edges),
                            ])
                elif self.graph_type == 'water_wire':
                    plot_folder = _hf.create_directory(self.workfolder+'/'+str(self.max_water)+'_water_wires/'+name+'/')
                    waters = '_max_'+str(self.max_water)+'_water_bridges' if self.max_water > 0 else ''
                    occ = '_min_occupancy_'+str(occupancy) if occupancy  else ''
                    for form in self.plot_parameters['formats']:
                        plt.savefig(f'{plot_folder}{name}{waters}{occ}_graph{is_label}.{form}', format=form, dpi=self.plot_parameters['plot_resolution'])
                    if is_label:
                        _hf.write_text_file(plot_folder+name+waters+occ+'_water_wire_graph_info.txt',
                            ['Water wire graph of '+name,
                            '\nSelection string: '+str(self.selection[0:-15]),
                            '\nNumber of maximum water molecules allowed in the bridge: '+str(self.max_water),
                            '\nMinimum H-bond occupancy: '+str(occupancy) if occupancy  else '',
                            '\n',
                            '\nNumber of nodes in '+name+': '+str(len(graph.nodes)),
                            '\nNumber of edges in '+name+': '+str(len(graph.edges)),
                            '\n',
                            '\nList of nodes: '+str(graph.nodes),
                            '\n',
                            '\nList of edges: '+str(graph.edges),
                            ])
                plt.close()


    def get_clusters(self):
        pass

    def plot_clusters(self):
        self.get_clusters()

    def get_linear_lenght(self, objects, graph):
        connected_components = _hf.get_connected_components(graph)
        struct_object = objects['structure'] if self.type_option == 'pdb' else objects['mda']
        return _hf.calculate_connected_compontents_coordinates(connected_components, struct_object, option=self.type_option)

    def plot_linear_lenghts(self, occupancy=None, label_nodes=True):
        self.logger.info('Plotting linear lengths for continuous network components')
        for name, objects in self.graph_coord_objects.items():
            self.logger.debug('Creating linear length plot for '+name)
            if 'graph' in objects.keys():
                if occupancy:
                    wba = copy.deepcopy(objects['wba'])
                    wba.filter_occupancy(occupancy)
                    graph = wba.filtered_graph
                else: graph = objects['graph']
                self.logger.debug('Creating '+self.graph_type+' linear length plot for: '+name)
                connected_components_coordinates = self.get_linear_lenght(objects, graph)
                plot_parameters = copy.deepcopy(self.plot_parameters)
                plot_parameters['figsize'] = (1+int(len(connected_components_coordinates)),plot_parameters['figsize'][1])

                plot_name = 'H-bond' if self.graph_type == 'hbond' else 'water wire'
                fig, ax = _hf.create_plot(title=f'Linear length of continuous {plot_name} subnetworks \nalong the Z-axis in structure {name}\nSelection: {self.selection[1:-16]}',
                                        xlabel='# of nodes',
                                        ylabel='Z-axis coordinates ($\AA$)',plot_parameters=plot_parameters)

                for i, g in enumerate(connected_components_coordinates):
                    for j in range(len(g)):
                        node_name = connected_components_coordinates[i][j][0]
                        chain_id, res_name, res_id = node_name.split('-')[0], node_name.split('-')[1], node_name.split('-')[2]
                        if res_name in ['HOH', 'TIP3']: color = self.plot_parameters['water_node_color']
                        elif res_name in  _hf.amino_d.keys(): color = self.plot_parameters['graph_color']
                        else: color = self.plot_parameters['non_prot_color']
                        z_coords = connected_components_coordinates[i][j][1][2]
                        ax.scatter(i, z_coords, color=color, s=self.plot_parameters['node_size']*0.8)
                        if label_nodes:
                            # if res_name in ['HOH', 'TIP3']: ax.annotate(f'W{res_id}', (i, z_coords), fontsize=self.plot_parameters['node_label_size'], zorder=6)
                            if res_name in _hf.amino_d.keys():
                                ax.annotate(f'{chain_id}-{_hf.amino_d[res_name]}{res_id}', (i, z_coords), fontsize=self.plot_parameters['node_label_size'], zorder=6)
                            else:
                                ax.annotate(f'{chain_id}-{res_name}{res_id}', (i, z_coords), fontsize=self.plot_parameters['node_label_size'], zorder=6)

                ax.set_xticks(np.arange(len(connected_components_coordinates)))
                ax.set_xticklabels([len(c) for c in connected_components_coordinates])

                plt.tight_layout()
                if self.graph_type == 'hbond':
                    plot_folder = _hf.create_directory(self.workfolder+'/H-bond_graphs/'+name+'/')
                    for form in self.plot_parameters['formats']:
                        plt.savefig(f'{plot_folder}{name}_H-bond_linear_length.{form}', format=form, dpi=self.plot_parameters['plot_resolution'])
                elif self.graph_type == 'water_wire':
                    plot_folder = _hf.create_directory(self.workfolder+'/'+str(self.max_water)+'_water_wires/'+name+'/')
                    waters = '_'+str(self.max_water)+'_water_bridges' if self.max_water > 0 else ''
                    occ = '_min_occupancy_'+str(occupancy) if occupancy  else ''
                    for form in self.plot_parameters['formats']:
                        plt.savefig(f'{plot_folder}{name}{waters}{occ}_linear_length.{form}', format=form, dpi=self.plot_parameters['plot_resolution'])
                plt.close()
            else: self.logger.warning(f'{name} has no {self.graph_type} graph. Linear length can not be calculated for this structure.')
