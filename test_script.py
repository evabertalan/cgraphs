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

  # 5: {# test squid structures
  #   'pdb_root_folder': '/Users/evabertalan/Documents/cgrap_test/squid',
  #   'reference_pdb': '/Users/evabertalan/Documents/cgrap_test/squid/2z73_sup.pdb',
  #   # 'test': ['water_cluster', 'water_wire', 'compare2']
  #   'test': ['water_cluster']
  # },

  # 6: {#Kappa opioid
  #   'pdb_root_folder': '/Users/evabertalan/Documents/cgrap_test/kappa',
  #   'reference_pdb': '/Users/evabertalan/Documents/cgrap_test/6b73_opm.pdb',
  #   # 'test': ['water_cluster', 'hbond', 'water_wire', 'compare2']
  #   'test': ['water_wire', 'compare2']
  # },

  7: {# test JSR1 sim
    'worfolder': '/Users/evabertalan/Documents/cgrap_test/jsr1_tests/',
    'psfs': ['/Users/evabertalan/Documents/cgrap_test/jsr1_tests/9cis_m103a/read_protein_membrane_7_9cis_m103a_3_2x.psf', '/Users/evabertalan/Documents/cgrap_test/jsr1_tests/9cis_optimized/read_protein_membrane_7_opt_3_2x.psf', '/Users/evabertalan/Documents/cgrap_test/jsr1_tests/9cis_y126a/read_protein_membrane_7_9cis_y126a_3_2x.psf'],
    'dcds': [['/Users/evabertalan/Documents/cgrap_test/jsr1_tests/9cis_m103a/9cis_m103a_last_20frames_pbc.dcd'], ['/Users/evabertalan/Documents/cgrap_test/jsr1_tests/9cis_optimized/9cis_optimized_last_20frames_pbc.dcd'], ['/Users/evabertalan/Documents/cgrap_test/jsr1_tests/9cis_y126a/9cis_y126a_last_20frames_pbc.dcd']],
    'names': ['m103a', '9cis_opt', 'y126a'],
    # 'test': ['sim_water_wire', 'sim_compare2']
    'test': ['sim_compare2']
  },

  # 8: {# test opioid MD
  #   'worfolder': '/Users/evabertalan/Documents/cgrap_test/opioid_kappa_md/',
  #   'psfs': ['/Users/evabertalan/Documents/cgrap_test/opioid_kappa_md/4djh/step5_assembly.xplor_ext.psf','/Users/evabertalan/Documents/cgrap_test/opioid_kappa_md/6b73/step5_assembly.xplor_ext.psf'],
  #   'dcds': [['/Users/evabertalan/Documents/cgrap_test/opioid_kappa_md/4djh/step7.20_production.dcd-pbc.dcd'], ['/Users/evabertalan/Documents/cgrap_test/opioid_kappa_md/6b73/step7.26_production.dcd-pbc.dcd']],
  #   'names': ['4djh', '6b73'],
  #   # 'test': ['sim_water_wire', 'sim_compare2']
  #   'test': ['sim_water_wire']
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
def test_MD_conserved_graph(psfs, dcds, names, target_folder):
  c_dcd = ConservedGraph(type_option='dcd', dcd_files=dcds, psf_files=psfs, sim_names=names, target_folder=target_folder)
  for max_water in [1, 3]:
    c_dcd.calculate_graphs(graph_type='water_wire', max_water=max_water)
    for occupancy in [0.1, 0.4, 0.8]:
      c_dcd.get_conserved_graph(conservation_threshold=0.8, occupancy=occupancy)
      c_dcd.plot_graphs(label_nodes=True, xlabel='PCA projected membrane plane', ylabel='Membrane normal ($\AA$)')
      c_dcd.plot_linear_lenghts()
      c_dcd.plot_conserved_graph(label_nodes=True, xlabel='PCA projected membrane plane', ylabel='Membrane normal ($\AA$)')
      c_dcd.plot_difference(label_nodes=True, xlabel='PCA projected membrane plane', ylabel='Membrane normal ($\AA$)')

# compare 2 static
def test_compare2_pdb(comp_type, pdb1, pdb2, target_folder, color1='#1b3ede',color2='#21c25f'):
  comp = CompareTwo('pdb', pdb1=pdb1, pdb2=pdb2, target_folder=target_folder)
  comp.calculate_graphs(graph_type=comp_type, max_water=4, include_backbone_sidechain=False)
  comp.plot_graph_comparison(color1=color1, color2=color2, label_nodes=True)
  comp.plot_graph_comparison(color1=color1, color2=color2, label_nodes=False)


# compare 2 MD
def test_compare2_MD(psf1, psf2, dcd1, dcd2, target_folder):
  comp = CompareTwo('dcd', psf1=psf1, psf2=psf2, dcd1=dcd1, dcd2=dcd2, target_folder=target_folder, name1='name1', name2='name2')
  for max_water in [1, 3]:
    comp.calculate_graphs(graph_type='water_wire', max_water=max_water, distance=3, cut_angle=59)
    for occupancy in [0.1, 0.4, 0.8]:
      comp.construct_comparison_objects(occupancy=occupancy)
      comp.plot_graph_comparison(color1='orange', color2='green', label_nodes=True)
      comp.plot_graph_comparison(color1='orange', color2='green', label_nodes=False)


for files in test_files.values():
  print('.'*50)
  print('TEST: ',files[list(files.keys())[0]])
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

  if 'sim_water_wire' in files['test']:
    test_MD_conserved_graph(files['psfs'], files['dcds'], files['names'], files['worfolder'])
  if 'sim_compare2' in files['test']:
    test_compare2_MD(files['psfs'][0], files['psfs'][1], files['dcds'][0], files['dcds'][1], files['worfolder'])


print('.'*40)
stop = timeit.default_timer()
print('Run time: ', stop - start)
print('TEST COMPLETED WITHOUT ISSUES')
print('='*40)
