from cgraph import WaterClusters
from cgraph import ConservedGraph


# pdb_root_folder = '/Users/evabertalan/Documents/protein_graph_analyzer/workfolder/test_files_GlplG'
# reference_pdb = '/Users/evabertalan/Documents/protein_graph_analyzer/workfolder/2irv_aout.pdb'


# pdb_root_folder = '/Users/evabertalan/Documents/JSR/rhodopsin_crystal_structures/squid'
# reference_pdb = '/Users/evabertalan/Documents/JSR/rhodopsin_crystal_structures/squid/2z73_sup.pdb'

# pdb_root_folder = '/Users/evabertalan/Documents/JSR/rhodopsin_crystal_structures/bovine/'
# reference_pdb = '/Users/evabertalan/Documents/JSR/rhodopsin_crystal_structures/bovine/1u19_sup.pdb'

# pdb_root_folder = '/Users/evabertalan/Documents/JSR/rhodopsin_crystal_structures/spider/'
# reference_pdb = '/Users/evabertalan/Documents/JSR/rhodopsin_crystal_structures/spider/6i9k_sup.pdb'

pdb_root_folder = '/Users/evabertalan/Documents/c_test_files/pdb_from_sim/'
reference_pdb = '/Users/evabertalan/Documents/c_test_files/pdb_from_sim/9cis_m103a.pdb'

# pdb_root_folder = '/Users/evabertalan/Documents/c_test_files/comp_2/squid'
# reference_pdb = '/Users/evabertalan/Documents/c_test_files/comp_2/squid/2z73_sup.pdb'
# pdb_root_folder = '/Users/evabertalan/Documents/c_test_files/test_new_features_GPCR/'
# reference_pdb = '/Users/evabertalan/Documents/c_test_files/4eiy_opm.pdb'

# pdb_root_folder = '/Users/evabertalan/Documents/c_test_files/adenosin_TEST/'
# reference_pdb = '/Users/evabertalan/Documents/c_test_files/4eiy_opm_clean.pdb'

# graph_type = 'water_wire'
graph_type = 'hbond'

#in the case of dcd analysis the ref file can be the initla crystal strucure or the pdb file of any frames
# maximum number of comperable simulations is 4 for conserved graph calculation, but if you want to compare more sim results ,you can export for example the last frame or multiple frames as pdb and than compare as pbd

w = WaterClusters(pdb_root_folder, reference_pdb=reference_pdb)
w.evaluate_parameters()
w.calculate_cluster_centers()
w.write_cluster_center_coordinates()
w.draw_clusters_centers_chimera()
ref_coordinates = w.reference_coordinates
c = ConservedGraph(pdb_root_folder, reference_pdb=reference_pdb, reference_coordinates=ref_coordinates)

# c = ConservedGraph(pdb_root_folder, reference_pdb=reference_pdb)
c.calculate_graphs(graph_type=graph_type)
c.plot_graphs(label_nodes=True)
c.plot_linear_lenghts()
c.get_conserved_graph()
c.plot_conserved_graph(label_nodes=True, xlabel='PCA projected membrane plane')
c.plot_conserved_graph(label_nodes=False, ylabel='Membrane normal ($\AA$)')
c.plot_difference(label_nodes=True, xlabel='PCA projected membrane plane')
c.plot_difference(label_nodes=False, ylabel='Membrane normal ($\AA$)')



# c = ConservedGraph(pdb_root_folder, type_option='pdb', reference_pdb=reference_pdb)
# c.calculate_graphs(graph_type=graph_type, selection='protein', max_water=3)
# c.plot_graphs(label_nodes=True, xlabel='PCA projected membrane plane')
# c.plot_graphs(label_nodes=False, ylabel='Membrane normal ($\AA$)')
# c.plot_linear_lenghts()
# c.get_conserved_graph()
# c.plot_conserved_graph(label_nodes=True, xlabel='PCA projected membrane plane')
# c.plot_conserved_graph(label_nodes=False, ylabel='Membrane normal ($\AA$)')
# c.plot_difference(label_nodes=True, xlabel='PCA projected membrane plane')
# c.plot_difference(label_nodes=False, ylabel='Membrane normal ($\AA$)')



# c = cg.ConservedGraph(pdb_root_folder, type_option='psdsddb', reference_pdb=reference_pdb)
# c.calculate_graphs(graph_type='water_wire', selection='protein or resname LYR', max_water=3)
# c.plot_graphs(label_nodes=True, xlabel='PCA projected membrane plane')
# c.plot_graphs(label_nodes=False, ylabel='Membrane normal ($\AA$)')
# c.plot_linear_lenghts()



# def _plot_conserved_graphs(c,occupancy=None):
#     c.plot_graphs(label_nodes=True, occupancy=0.9)
#     c.plot_graphs(label_nodes=False, occupancy=0.1)
#     c.plot_linear_lenghts(occupancy=0.9)
#     c.get_conserved_graph(conservation_threshold=0.90, occupancy=0.6)
#     c.plot_conserved_graph(label_nodes=True)
#     c.plot_conserved_graph(label_nodes=False)
#     c.plot_difference(label_nodes=True)
#     c.plot_difference(label_nodes=False)
#     # c.logger.info('Calculation completed\n'+'-'*20)
#     c.logger.info('Calculation completed\n'+'-'*20)


# file = ['/Users/evabertalan/Documents/protein_graph_analyzer/test_trajs/workfolder/.graph_objects/jsr1_opt_water_wire_wba_object.pickle', '/Users/evabertalan/Documents/protein_graph_analyzer/test_trajs/workfolder/.graph_objects/sim1_l20_water_wire_wba_object.pickle', '/Users/evabertalan/Documents/protein_graph_analyzer/test_trajs/workfolder/.graph_objects/jsr1_y126a_water_wire_wba_object.pickle']
# _target_folder='/Users/evabertalan/Documents/protein_graph_analyzer/test_trajs/'

# c_dcd = ConservedGraph(type_option='dcd', sequance_identity_threshold=0.75, target_folder=_target_folder)
# c_dcd._load_exisitng_graphs(graph_files=file ,graph_type='water_wire')
# _plot_conserved_graphs(c_dcd,occupancy=0.8)



