import os
import ast
import argparse
import glob
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.transformations import wrap
import matplotlib.pyplot as plt

import pdb

class DistOverTraj:
    def __init__(self, psf, dcd, output_folder):

        _, ext = os.path.splitext(psf)
        if ext != '.psf':
            raise ValueError(f'The first argument has to be a PSF file. The provided file is a {ext}.')

        self.psf = psf
        self.dcd = dcd
        self.u = mda.Universe(self.psf, self.dcd)
        self.u.trajectory.add_transformations(wrap(self.u.atoms))

        base = os.path.basename(psf)
        self.base_name, ext = os.path.splitext(base)
        self.output_folder = output_folder


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

        bonding_groups = ('OH2', 'OW', 'NE', 'NH1', 'NH2', 'ND2', 'SG', 'NE2', 'ND1', 'NZ', 'OG', 'OG1', 'NE1', 'OH', 'OE1', 'OE2', 'N16', 'OD1', 'OD2', 'SD')
        selection_bonding_group = ' or name '.join(bonding_groups)

        if len(self.edges[0][0].split('-')) == 3:
            # if Bridge analysis was done with residuewise=True, for the distance calculation we take the CA atom
            # not recommended for distance calumniation, using residuewise=False gives the distance we are interested in
            e1 = [self.u.select_atoms(f"segid {edge[0].split('-')[0]} and resname {edge[0].split('-')[1]} and resid {edge[0].split('-')[2]} and name CA") for edge in self.edges]
            e2 = [self.u.select_atoms(f"segid {edge[1].split('-')[0]} and resname {edge[1].split('-')[1]} and resid {edge[1].split('-')[2]} and name CA") for edge in self.edges]

        else: #for residuewise=False, distance can be accurately calculated between the actual H-bonding atom pairs
            e1 = [self.u.select_atoms(f"segid {edge[0].split('-')[0]} and resname {edge[0].split('-')[1]} and resid {edge[0].split('-')[2]} and name {edge[0].split('-')[3]}") for edge in self.edges]
            e2 = [self.u.select_atoms(f"segid {edge[1].split('-')[0]} and resname {edge[1].split('-')[1]} and resid {edge[1].split('-')[2]} and name {edge[1].split('-')[3]}") for edge in self.edges]

        waters = self.u.select_atoms('resname TIP3 and name OH2')
        self.distances_over_time = []
        self.waters_over_time = []
        self.total_waters_around_res_over_time = []
        threshold_distance = 3.5
        box = self.u.dimensions

        #add option to only read every 10th frame
        for ts in self.u.trajectory:
            frame_distances = []
            water_numbers = []
            total_waters_around_res = []
            for res1, res2 in zip(e1, e2):
                box = self.u.dimensions
                dist = distance_array(res1.positions, res2.positions, box=box)[0][0]
                frame_distances.append(dist)

                res1_w = (distance_array(res1.positions, waters.positions, box=box) < threshold_distance).sum()
                res2_w = (distance_array(res2.positions, waters.positions, box=box) < threshold_distance).sum()
                water_numbers.append(f'{res1_w} - {res2_w}')

                selected_water_numb = []
                for r in [res1, res2]:
                    sstring = f'(segid {r.segids[0]} and resname {r.resnames[0]} and resid {r.resids[0]}) and (name {selection_bonding_group})'
                    selected_water_numb.append(len(self.u.select_atoms(f'(resname TIP3 and name OH2) and (around {threshold_distance} {sstring})')))
                total_waters_around_res.append(f'{selected_water_numb[0]} - {selected_water_numb[1]}')


            self.distances_over_time.append(frame_distances)
            self.waters_over_time.append(water_numbers)
            self.total_waters_around_res_over_time.append(total_waters_around_res)


    def write_results_to_df(self):
        self.distances_df = pd.DataFrame(np.array(self.distances_over_time), columns=[f'{e[0]} - {e[1]}' for e in self.edges])
        # add frame as index column
        self.distances_df.to_csv(f'{self.output_folder}/{self.base_name}_pair_distances.csv', index=False)

        self.waters_df = pd.DataFrame(np.array(self.waters_over_time), columns=[f'{e[0]} - {e[1]}' for e in self.edges])
        self.waters_df.to_csv(f'{self.output_folder}/{self.base_name}_waters_within_3_5.csv', index=False)

        self.total_water_df = pd.DataFrame(np.array(self.total_waters_around_res_over_time), columns=[f'{e[0]} - {e[1]}' for e in self.edges])
        self.total_water_df.to_csv(f'{self.output_folder}/{self.base_name}_all_waters_within_3_5.csv', index=False)

    def plot_results(self):
        self.distances_df = pd.read_csv(f'{self.output_folder}/{self.base_name}_pair_distances.csv')
        fig, axes = plt.subplots(nrows=len(self.distances_df.columns), ncols=1, figsize=(20, 80), sharex=True)

        for i, column in enumerate(self.distances_df.columns):
            axes[i].plot(self.distances_df.index, self.distances_df[column])
            axes[i].spines['right'].set_visible(False)
            axes[i].spines['top'].set_visible(False)
            axes[i].set_title(f"{column}")
            axes[i].set_ylabel("Distance")

        # Set the x-axis label only on the last subplot
        axes[-1].set_xlabel("Frame")

        # Adjust layout to prevent overlap
        plt.tight_layout()
        fig.savefig('./dist_time_series_per_res.png')


# def process_input():
# cgraphs_input

#         if cgraphs_input:
#             if not cgraphs_input.endswith('_info.txt'):
#                 raise ValueError(f'The --cgraphs_input path has to point to the _info.txt output file of cgrpahs.')
#             else:
#                 residue_selection_string = self._read_cgraphs_input(cgraphs_input)
#                 self._u = self._u.select_atoms(residue_selection_string)

def main():
    parser = argparse.ArgumentParser(description='Calculate distances from molecular dynamics trajectories.')
    parser.add_argument('psf', help='Path to the PSF file')
    parser.add_argument('dcd', nargs='+', help='Path to the DCD files. The path can contain regex to select multiple files by a name pattern.')
    parser.add_argument('--cgraphs_input', help='Path to an _info.txt cgraphs file, which will be read as input for the distance calculation.')
    parser.add_argument('--output_folder', help='Path to the output file for pKa data')
    #add selection input
    parser.add_argument('--selection', default='protein', help='Atom selection for pKa calculation in MD Analysis syntax.')


    #clarify the relationship in the code: which is selected first C-Graphs or selection? --> if someone is using cgraphs input, all thos residues will be selected and these nodes can be further restricted with the --selection argument.
    #if only --selection is given, those are the selected residues

    args = parser.parse_args()


    dcd_files = []
    for dcd_file in args.dcd:
            dcd_files += glob.glob(dcd_file)
    dcd_files.sort()

    dist_traj = DistOverTraj(args.psf, dcd_files, args.output_folder)
    dist_traj.calculate_distances(args.cgraphs_input)
    dist_traj.write_results_to_df()
    # dist_traj.plot_results()


if __name__ == '__main__':
    main()
