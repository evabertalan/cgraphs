import os
import ast
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.transformations import wrap
import pdb

PSF = '/Users/evabertalan/Documents/projects/JSR1_pka_proj/scripts/example_files/test_traj/read_protein_membrane_7_9cis_y126a_3_2x.psf'
DCD = '/Users/evabertalan/Documents/projects/JSR1_pka_proj/scripts/example_files/test_traj/9cis_y126a_last_20frames_pbc.dcd'

CG_INPUT = '/Users/evabertalan/Documents/projects/JSR1_pka_proj/scripts/example_files/test_traj/workfolder/3_water_wires/sim_atomwise/sim_atomwise_max_3_water_bridges_min_occupancy_0.1_water_wire_graph_info.txt'

class DistOverTraj:
    def __init__(self, psf, dcd):

        _, ext = os.path.splitext(psf)
        if ext != '.psf':
            raise ValueError(f'The first argument has to be a PSF file. The provided file is a {ext}.')

        self.psf = psf
        self.dcd = dcd
        self.u = mda.Universe(self.psf, self.dcd)
        self.u.trajectory.add_transformations(wrap(self.u.atoms))


    def _read_cgraphs_input(self, cgraphs_input):
        if os.path.exists(cgraphs_input):
            with open(cgraphs_input) as f:
                for line in f.readlines():
                    if line.startswith('List of edges:'):
                        edges = ast.literal_eval(line.split(': ')[-1])
                return edges

        else:
            raise FileNotFoundError(f'The {cgraphs_input} file does not exist.')

    def calculate_distances(self, cgraphs_input):
        self.edges = self._read_cgraphs_input(cgraphs_input)

        # e1 = [self.u.select_atoms(f"segid {edge[0].split('-')[0]} and resname {edge[0].split('-')[1]} and resid {edge[0].split('-')[2]} and name CA") for edge in self.edges]
        # e2 = [self.u.select_atoms(f"segid {edge[1].split('-')[0]} and resname {edge[1].split('-')[1]} and resid {edge[1].split('-')[2]} and name CA") for edge in self.edges]

        e1 = [self.u.select_atoms(f"segid {edge[0].split('-')[0]} and resname {edge[0].split('-')[1]} and resid {edge[0].split('-')[2]} and name {edge[0].split('-')[3]}") for edge in self.edges]
        e2 = [self.u.select_atoms(f"segid {edge[1].split('-')[0]} and resname {edge[1].split('-')[1]} and resid {edge[1].split('-')[2]} and name {edge[1].split('-')[3]}") for edge in self.edges]

        self.distances_over_time = []
        #add option to only read every 10th frame
        for ts in self.u.trajectory:
            frame_distances = []
            for res1, res2 in zip(e1, e2):
                dist = distance_array(res1.positions, res2.positions)[0][0]
                frame_distances.append(dist)
            self.distances_over_time.append(frame_distances)


    def write_results_to_df(self):
        distances_df = pd.DataFrame(np.array(self.distances_over_time), columns=[f'{e[0]} - {e[1]}' for e in self.edges])
        distances_df.to_csv("residue_pair_distances.csv", index=False)


    def plot_results(self):
        pass
# def process_input():
# cgraphs_input

#         if cgraphs_input:
#             if not cgraphs_input.endswith('_info.txt'):
#                 raise ValueError(f'The --cgraphs_input path has to point to the _info.txt output file of cgrpahs.')
#             else:
#                 residue_selection_string = self._read_cgraphs_input(cgraphs_input)
#                 self._u = self._u.select_atoms(residue_selection_string)


dist_traj = DistOverTraj(PSF, DCD)
dist_traj.calculate_distances(CG_INPUT)
dist_traj.write_results_to_df()
dist_traj.plot_results()

