import os
from cgraph import WaterClusters
from cgraph import ConservedGraph
from cgraph import CompareTwo
import timeit

start = timeit.default_timer()

test_files = {
  # 0: {# simple test
  #   'pdb_root_folder': '/Users/evabertalan/Documents/cgrap_test/GlplG',
  #   'reference_pdb': '/Users/evabertalan/Documents/cgrap_test/GlplG/2irv_b.pdb',
  #   'test': ['water_cluster', 'hbond', 'water_wire', 'compare2']
  # },
  # 1: {# test adenosine GPCRs
  #   'pdb_root_folder': '/Users/evabertalan/Documents/cgrap_test/adenosin_TEST',
  #   'reference_pdb': '/Users/evabertalan/Documents/cgrap_test/4eiy_opm_clean.pdb',
  #   'test': ['water_cluster', 'hbond', 'water_wire', 'compare2']
  # },

  # 2: {# test multiple chains monomer --> SARS-COV
  #   'pdb_root_folder': '/Users/evabertalan/Documents/cgrap_test/cov_test',
  #   'reference_pdb': '/Users/evabertalan/Documents/cgrap_test/cov_test/6m17A_sup.pdb',
  #   'test': ['hbond', 'compare2']
  # },

  # 3: {# test multiple chains dimer--> SARS-COV
  #   'pdb_root_folder': '/Users/evabertalan/Documents/cgrap_test/corona_structs',
  #   'reference_pdb': '/Users/evabertalan/Documents/cgrap_test/corona_structs/6lzg.pdb',
  #   'test': ['hbond', 'compare2']
  # },

  # 4: {# test bovine structures
  #   'pdb_root_folder': '/Users/evabertalan/Documents/cgrap_test/bovine',
  #   'reference_pdb': '/Users/evabertalan/Documents/cgrap_test/bovine/1u19_sup.pdb',
  #   'test': ['water_cluster', 'hbond', 'water_wire', 'compare2']
  # },

  5: {# test squid structures
    'pdb_root_folder': '/Users/evabertalan/Documents/cgrap_test/squid',
    'reference_pdb': '/Users/evabertalan/Documents/cgrap_test/squid/2z73_sup.pdb',
    # 'test': ['water_cluster', 'water_wire', 'compare2']
    'test': ['water_cluster']
  },

  # 6: {#Kappa opioid
  #   'pdb_root_folder': '/Users/evabertalan/Documents/cgrap_test/kappa',
  #   'reference_pdb': '/Users/evabertalan/Documents/cgrap_test/6b73_opm.pdb',
  #   # 'test': ['water_cluster', 'hbond', 'water_wire', 'compare2']
  #   'test': ['water_wire', 'compare2']
  # },

  # 7: {# test JSR1 sim
  #   'worfolder': '/Users/evabertalan/Documents/cgrap_test/jsr1_tests/',
  #   'psf1': '/Users/evabertalan/Documents/cgrap_test/jsr1_tests/9cis_m103a/read_protein_membrane_7_9cis_m103a_3_2x.psf',
  #   'dcd1': ['/Users/evabertalan/Documents/cgrap_test/jsr1_tests/9cis_m103a/9cis_m103a_last_20frames_pbc.dcd'],

  #   'psf2': '/Users/evabertalan/Documents/cgrap_test/jsr1_tests/9cis_optimized/read_protein_membrane_7_opt_3_2x.psf',
  #   'dcd2': ['/Users/evabertalan/Documents/cgrap_test/jsr1_tests/9cis_optimized/9cis_optimized_last_20frames_pbc.dcd'],

  #   'psf3': '/Users/evabertalan/Documents/cgrap_test/jsr1_tests/9cis_y126a/read_protein_membrane_7_9cis_y126a_3_2x.psf',
  #   'dcd3': ['/Users/evabertalan/Documents/cgrap_test/jsr1_tests/9cis_y126a/9cis_y126a_last_20frames_pbc.dcd'],
  #   'test': ['sim_water_wire', 'sim_compare2']
  # },

  # 8: {# test opioid MD
  #   'psf1': '/Users/evabertalan/Documents/cgrap_test/opioid_kappa_md/4djh/step5_assembly.xplor_ext.psf',
  #   'dcd1': ['/Users/evabertalan/Documents/cgrap_test/opioid_kappa_md/4djh/step7.20_production.dcd-pbc.dcd'],

  #   'psf2': '/Users/evabertalan/Documents/cgrap_test/opioid_kappa_md/6b73/step5_assembly.xplor_ext.psf',
  #   'dcd2': ['/Users/evabertalan/Documents/cgrap_test/opioid_kappa_md/6b73/step7.26_production.dcd-pbc.dcd'],
  #   'test': ['sim_water_wire', 'sim_compare2']
  # },

  # 9: {# PDB created from simulation
  #   'pdb_root_folder': '/Users/evabertalan/Documents/cgrap_test/sim_pdb',
  #   'reference_pdb': '/Users/evabertalan/Documents/cgrap_test/sim_pdb/9cis_m103a.pdb',
  #   # 'test': ['water_cluster', 'water_wire', 'compare2']
  #   'test': ['compare2']
  # },

  # 10: {# test big molecule
  #   'pdb_root_folder': '/Users/evabertalan/Documents/cgrap_test/big_comp',
  #   'reference_pdb': '/Users/evabertalan/Documents/cgrap_test/big_comp/6vqr_sup.pdb',
  #   'test': ['hbond', 'compare2']
  # },
  # 11: {# test spider molecule -- what happens with only 1 structures
  #   'pdb_root_folder': '/Users/evabertalan/Documents/cgrap_test/spider',
  #   'reference_pdb': '/Users/evabertalan/Documents/cgrap_test/spider/6i9k_sup.pdb',
  #   'test': ['water_cluster', 'hbond', 'water_wire']
  # }
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
def test_H_bond(pdb_root_folder, reference_pdb):
  graph_type = 'hbond'
  c = ConservedGraph(pdb_root_folder, reference_pdb=reference_pdb)
  c.calculate_graphs(graph_type=graph_type)
  c.plot_graphs(label_nodes=True)
  c.plot_linear_lenghts()
  c.get_conserved_graph()
  c.plot_conserved_graph(label_nodes=True, xlabel='PCA projected membrane plane')
  c.plot_conserved_graph(label_nodes=False, ylabel='Membrane normal ($\AA$)')
  c.plot_difference(label_nodes=True, xlabel='PCA projected membrane plane')
  c.plot_difference(label_nodes=False, ylabel='Membrane normal ($\AA$)')


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
def test_compare2_pdb(comp_type, pdb1, pdb2, target_folder, color1='#1b3ede',color2='#21c25f'):
  comp = CompareTwo('pdb', pdb1=pdb1, pdb2=pdb2, target_folder=target_folder)
  comp.calculate_graphs(graph_type=comp_type, max_water=4, include_backbone_sidechain=False)
  comp.plot_graph_comparison(color1=color1, color2=color2, label_nodes=True)
  comp.plot_graph_comparison(color1=color1, color2=color2, label_nodes=False)


# compare 2 MD
def test_compare2_MD():
  pass


for files in test_files.values():
  print('.'*50)
  print('TEST: ',files['pdb_root_folder'])
  if 'water_cluster' in files['test']:
    test_H_bond_with_water_clusters(files['pdb_root_folder'], files['reference_pdb'])
  if 'hbond' in files['test']:
    test_H_bond(files['pdb_root_folder'], files['reference_pdb'])
  if 'water_wire' in files['test']:
    test_PDB_water_wire(files['pdb_root_folder'], files['reference_pdb'])
  if 'compare2' in files['test']:
    target_folder = files['pdb_root_folder']
    pdb1 = files['pdb_root_folder']+'/'+[file for file in os.listdir(files['pdb_root_folder']) if file.endswith('.pdb')][0]
    pdb2 = files['pdb_root_folder']+'/'+[file for file in os.listdir(files['pdb_root_folder']) if file.endswith('.pdb')][1]
    test_compare2_pdb('hbond', pdb1, pdb2, target_folder)
    test_compare2_pdb('water_wire', pdb1, pdb2, target_folder)



print('.'*40)
stop = timeit.default_timer()
print('Run time: ', stop - start)
print('TEST COMPLETED WITHOUT ISSUES')
print('='*40)
