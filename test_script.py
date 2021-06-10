from cgraph import WaterClusters
from cgraph import ConservedGraph

test_files = {
  0: {# simple test
    'pdb_root_folder': '/Users/evabertalan/Documents/cgrap_test/GlplG',
    'reference_pdb': '/Users/evabertalan/Documents/cgrap_test/GlplG/2irv_b.pdb'
  },
  1: {# test adenosine GPCRs
    'pdb_root_folder': '/Users/evabertalan/Documents/cgrap_test/adenosin_TEST',
    'reference_pdb': '/Users/evabertalan/Documents/cgrap_test/4eiy_opm_clean.pdb'
  },


# # test multiple chains monomer --> SARS-COV
# '/Users/evabertalan/Documents/cgrap_test/cov_test'
# '/Users/evabertalan/Documents/cgrap_test/cov_test/6m17A_sup.pdb'

# # test multiple chains dimer--> SARS-COV
# '/Users/evabertalan/Documents/cgrap_test/corona_structs'
# '/Users/evabertalan/Documents/cgrap_test/corona_structs/6lzg.pdb'

# # test bovine structures
# '/Users/evabertalan/Documents/cgrap_test/bovine'
# '/Users/evabertalan/Documents/cgrap_test/bovine/1u19_sup.pdb'

# # test squid structures
# '/Users/evabertalan/Documents/cgrap_test/squid'
# '/Users/evabertalan/Documents/cgrap_test/squid/2z73_sup.pdb'

# # test just 1 structure --> spider
# '/Users/evabertalan/Documents/cgrap_test/spider'
# '/Users/evabertalan/Documents/cgrap_test/spider/6i9k_sup.pdb'

# #Kappa opioid
# '/Users/evabertalan/Documents/cgrap_test/kappa'
# '/Users/evabertalan/Documents/cgrap_test/6b73_opm.pdb'

# # test JSR1 sim
# worfolder = '/Users/evabertalan/Documents/cgrap_test/jsr1_tests/'
# psf1 '/Users/evabertalan/Documents/cgrap_test/jsr1_tests/9cis_m103a/read_protein_membrane_7_9cis_m103a_3_2x.psf'
# dcd1 ['/Users/evabertalan/Documents/cgrap_test/jsr1_tests/9cis_m103a/9cis_m103a_last_20frames_pbc.dcd']

# psf2 '/Users/evabertalan/Documents/cgrap_test/jsr1_tests/9cis_optimized/read_protein_membrane_7_opt_3_2x.psf'
# dcd2 ['/Users/evabertalan/Documents/cgrap_test/jsr1_tests/9cis_optimized/9cis_optimized_last_20frames_pbc.dcd']

# psf3 '/Users/evabertalan/Documents/cgrap_test/jsr1_tests/9cis_y126a/read_protein_membrane_7_9cis_y126a_3_2x.psf'
# dcd3 ['/Users/evabertalan/Documents/cgrap_test/jsr1_tests/9cis_y126a/9cis_y126a_last_20frames_pbc.dcd']

# # test opioid MD
# psf1 '/Users/evabertalan/Documents/cgrap_test/opioid_kappa_md/4djh/step5_assembly.xplor_ext.psf'
# dcd1 ['/Users/evabertalan/Documents/cgrap_test/opioid_kappa_md/4djh/step7.20_production.dcd-pbc.dcd']

# psf2 '/Users/evabertalan/Documents/cgrap_test/opioid_kappa_md/6b73/step5_assembly.xplor_ext.psf'
# dcd2 ['/Users/evabertalan/Documents/cgrap_test/opioid_kappa_md/6b73/step7.26_production.dcd-pbc.dcd']


# # test big molecule
# '/Users/evabertalan/Documents/cgrap_test/big_comp'
# '/Users/evabertalan/Documents/cgrap_test/big_comp/6vqr_sup.pdb'
}

# I. water cluster and conserved H-bon graph test
def test_H_bond_with_water_clusters(pdb_root_folder, reference_pdb):
  graph_type = 'hbond'
  w = WaterClusters(pdb_root_folder, reference_pdb=reference_pdb)
  if len(w.superimposed_files) > 1:
    w.evaluate_parameters()
    w.calculate_cluster_centers()
    w.write_cluster_center_coordinates()
    w.draw_clusters_centers_chimera()
    ref_coordinates = w.reference_coordinates
    c = ConservedGraph(pdb_root_folder, reference_pdb=reference_pdb, reference_coordinates=ref_coordinates)
  else:   c = ConservedGraph(pdb_root_folder, reference_pdb=reference_pdb)

  c.calculate_graphs(graph_type=graph_type)
  c.plot_graphs(label_nodes=True)
  c.plot_linear_lenghts()
  c.get_conserved_graph()
  c.plot_conserved_graph(label_nodes=True, xlabel='PCA projected membrane plane')
  c.plot_conserved_graph(label_nodes=False, ylabel='Membrane normal ($\AA$)')
  c.plot_difference(label_nodes=True, xlabel='PCA projected membrane plane')
  c.plot_difference(label_nodes=False, ylabel='Membrane normal ($\AA$)')


# II. conserved water wire static
def test_PDB_water_wire(pdb_root_folder, reference_pdb):
  graph_type = 'water_wire'
  c = ConservedGraph(pdb_root_folder, reference_pdb=reference_pdb)
  c.calculate_graphs(graph_type=graph_type, selection='protein', max_water=3)
  c.plot_graphs(label_nodes=True, xlabel='PCA projected membrane plane')
  c.plot_graphs(label_nodes=False, ylabel='Membrane normal ($\AA$)')
  c.plot_linear_lenghts()
  c.get_conserved_graph()
  c.plot_conserved_graph(label_nodes=True, xlabel='PCA projected membrane plane')
  c.plot_conserved_graph(label_nodes=False, ylabel='Membrane normal ($\AA$)')
  c.plot_difference(label_nodes=True, xlabel='PCA projected membrane plane')
  c.plot_difference(label_nodes=False, ylabel='Membrane normal ($\AA$)')


# III. conserved water wire MD
def test_MD_conserved_graph(occupancy):
  pass


# compare 2 static
def test_compare2_pdb():
  pass


# compare 2 MD
def test_compare2_MD():
  pass


for files in test_files.values():
  test_H_bond_with_water_clusters(files['pdb_root_folder'], files['reference_pdb'])
  test_PDB_water_wire(files['pdb_root_folder'], files['reference_pdb'])

print('.'*40)
print('TEST COMPLETED WITHOUT ISSUES')
print('='*40)
