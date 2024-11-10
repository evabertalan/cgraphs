import os
import argparse
import ast
import glob
from propkatraj import PropkaTraj
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import warnings

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    warnings.filterwarnings("ignore", category=UserWarning)
    import MDAnalysis as mda

class PkaFromTraj:
    def __init__(self, psf, dcd, cgraphs_input):
        """
        Initialize PkaFromTraj with a PSF and DCD file.

        Parameters:
        psf (str): Path to the PSF file.
        dcd (list): A list with the paths to the DCD files.
        """

        _, ext = os.path.splitext(psf)
        if ext != '.psf':
            raise ValueError(f'The first argument has to be a PSF file. The provided file is a {ext}.')

        self.psf = psf
        self.dcd = dcd
        self._u = mda.Universe(self.psf, self.dcd)
        self.pkas = pd.DataFrame()

        if cgraphs_input:
            if not cgraphs_input.endswith('_info.txt'):
                raise ValueError(f'The --cgraphs_input path has to point to the _info.txt output file of cgrpahs.')
            else:
                residue_selection_string = self._read_cgraphs_input(cgraphs_input)
                self._u = self._u.select_atoms(residue_selection_string)
        

    def print_selection_info(self):
        print('Number of atoms in selection:', len(self._u.atoms))
        print('Number of residues in selection:',len(self._u.residues))
        print('Number of frames:', len(self._u.trajectory))

    def compute_pka_for_traj(self, selection='protein', start=None, stop=None, step=None):
        """
        Compute pKa values for the trajectory.
        
        Parameters:
        selection (str): Atom selection for pKa calculation.
        start (int): Starting frame index.
        stop (int): Stopping frame index.
        step (int): Step between frames.
        """
        self.selection = selection

        self._u.add_TopologyAttr('record_types')
        self._u = self._u.select_atoms(self.selection)

        hetatm_residues = self._u.select_atoms('not protein')
        for atom in hetatm_residues:
            atom.record_type = 'HETATM'

        hse_atoms = self._u.select_atoms('resname HSE')
        for hse in hse_atoms:
            hse.residue.resname = 'HIS'

        try:
            self.pkatraj = PropkaTraj(self._u, select=selection, skip_failure=True)
            self.pkatraj.run(start, stop, step)
            self.pkas = self.pkatraj.results.pkas.sort_values(by='time')
            
        except Exception as e:
            print(f'Error computing pKa for trajectory: {e}')

    def get_pka_for_frame(self, write_to_file=None):
        """
        Retrieve pKa values for the trajectory frame.
        
        Parameters:
        write_to_file (str): If provided, writes the pKa values to the specified file.
        
        Returns:
        DataFrame: pKa values for the current frame.
        """

        if write_to_file:
            residues = [self._u.select_atoms(f'resid {col}').residues[0] for col in self.pkas.columns]
            col_names = [f'{res.segid}-{res.resname}-{res.resid}'for res in residues]
            df = self.pkas.copy()
            df.columns = col_names
            df.to_csv(write_to_file)
        return self.pkas

    def get_pka_statistic(self, selection=None, write_to_file=None):
        """
        Get statistical description of pKa values for the selected residues.
        
        Parameters:
        selection (str): Atom selection string for computing statistics.
        write_to_file (str): If provided, writes the statistics to the specified file.
        
        Returns:
        DataFrame: Statistical description of pKa values.
        """
        if self.pkas.empty:
            print('No pKa data available. Please run compute_pka_for_traj first.')
            return None
        resids = self._get_selected_res(selection)
        stats = self.pkas[resids].describe()
        if write_to_file:
            stats.to_csv(write_to_file)
        return stats

    def plot_pka_time_series_for_selection(self, selection=None, write_to_file=None, figsize=None):
        if self.pkas.empty:
            print('No pKa data available. Please run compute_pka_for_traj first.')
            return None, None
        resids = self._get_selected_res(selection)

        title = selection if selection else self.selection
        fig, ax = self._create_plot(title=title, xlabel='Time (ns)', ylabel=r'p$K_a$', figsize=figsize)
        
        for col in resids:
            ax.plot(self.pkas[col], label=col)
        ax.legend()

        if write_to_file:
            fig.savefig(write_to_file)
        return fig, ax

    def plot_pka_statistic_for_selection(self, selection=None, write_to_file=None, figsize=None):
        if self.pkas.empty:
            print('No pKa data available. Please run compute_pka_for_traj first.')
            return None, None
        resids = self._get_selected_res(selection)
        
        title = selection if selection else self.selection
        fig, ax = self._create_plot(title=title, xlabel='Res ID', ylabel=r'p$K_a$', figsize=figsize)
        sns.boxplot(data=self.pkas[resids], ax=ax)

        if write_to_file:
            fig.savefig(write_to_file)
        return fig, ax

    def write_pka_to_external_data_file(self, write_to_file):
        if self.pkas.empty:
            print('No pKa data available. Please run compute_pka_for_traj first.')
            return None
        assert write_to_file.endswith('_data.txt'), 'Name of external data file should end with _data.txt'
        
        stats = self.get_pka_statistic()

        res_name_map = {}
        prot = self._u.select_atoms(self.selection)
        for res in prot.residues:
            res_id = res.resid
            res_name = res.resname
            seg_id = res.segid
            res_name_map.update({f'{res_id}': {'res_name': res_name, 'seg_id': seg_id}})

        with open(write_to_file, 'w') as f:
            for res_id, avg_pka in zip(stats.columns, stats.iloc[1]):
                line = ' '.join([res_name_map[str(res_id)]['res_name'], str(res_id), res_name_map[str(res_id)]['seg_id'], str(round(avg_pka, 3))])
                f.write(line + '\n')

    def write_last_frame_as_pdb(self, write_to_file):
        assert write_to_file.endswith('.pdb'), 'File name has to end with .pdb'
        
        self._u.trajectory[-1]
        with mda.Writer(write_to_file, multiframe=True, bonds=None, n_atoms=self._u.atoms.n_atoms) as PDB:
            PDB.write(self._u.atoms)

    def _get_selected_res(self, selection):
        try:
            if selection and selection != 'protein':
                selected_atoms = self._u.select_atoms(selection)
                if not selected_atoms:
                    raise ValueError(f'No atoms found for selection: {selection}')
                resids = selected_atoms.residues.resids
            else:
                resids = self.pkas.columns
            return resids
        except Exception as e:
            raise ValueError(f'Error in getting selected residues: {e}')

    def _read_cgraphs_input(self, cgraphs_input):
        if os.path.exists(cgraphs_input):
            with open(cgraphs_input) as f:
                for line in f.readlines():
                    if line.startswith('List of nodes:'):
                        nodes = ast.literal_eval(line.split(': ')[-1])

            selection = []
            for n in nodes:
                seg_id, res_name, res_id = n.split('-')
                substr = f'(segid {seg_id} and resname {res_name} and resid {res_id})'
                selection.append(substr)
            selection_str = ' or '.join(selection)
            return selection_str
        else:
            raise FileNotFoundError(f"The {cgraphs_input} file does not exist.")

    def _create_plot(self, title='', xlabel='', ylabel='', figsize=None):
        figsize = figsize if figsize else (12, 8)
        fig, ax = plt.subplots(figsize=figsize)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_title(title, fontsize=12)
        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        return fig, ax


def main():
    parser = argparse.ArgumentParser(description='Process pKa from molecular dynamics trajectories.')
    parser.add_argument('psf', help='Path to the PSF file')
    parser.add_argument('dcd', nargs='+', help='Path to the DCD files. The path can contain regex to select multiple files by a name pattern.')
    parser.add_argument('--selection', default='protein', help='Atom selection for pKa calculation in MD Analysis syntax.')
    parser.add_argument('--start', type=int, help='Starting frame index')
    parser.add_argument('--stop', type=int, help='Stopping frame index')
    parser.add_argument('--step', type=int, help='Step between frames')
    parser.add_argument('--output_folder', help='Path to the output file for pKa data')
    parser.add_argument('--plot', action='store_true', help='Plot time_series and statistic')
    parser.add_argument('--cgraphs_input', help='Path to an _info.txt cgraphs file, which will be read as input for the pKa calculation. pKa values are calculated for titrable residues of the calculated graph nodes..')

    #clarify the relationship in the code: which is selected first C-Graphs or selection? --> if someone is using cgraphs input, all thos residues will be selected and these nodes can be further restricted with the --selection argumetn.
    #if only --selection is given, those are the selected residues
    args = parser.parse_args()
    base = os.path.basename(args.psf)
    base_name, ext = os.path.splitext(base)

    output_folder = args.output_folder if args.output_folder else os.path.dirname(args.psf)
    #check if output folder exists, then creati t

    dcd_files = []
    for dcd_file in args.dcd:
            dcd_files += glob.glob(dcd_file)
    dcd_files.sort()

    pka_traj = PkaFromTraj(args.psf, dcd_files, args.cgraphs_input)
    pka_traj.compute_pka_for_traj(args.selection, args.start, args.stop, args.step)

    pka_frame_filename = os.path.join(output_folder, f'pkas_for_frames_{base_name}.csv')
    pka_traj.get_pka_for_frame(write_to_file=pka_frame_filename)

    stats_filename = os.path.join(output_folder, f'pkas_stats_{base_name}.csv')
    pka_traj.get_pka_statistic(write_to_file=stats_filename)

    external_data_file_name =  os.path.join(output_folder, f'{base_name}_data.txt')
    pka_traj.write_pka_to_external_data_file(external_data_file_name)

    if args.plot:
            time_series_plot_name =  os.path.join(output_folder, f"ts_{args.selection.join('_')}_{base_name}.png")
            pka_traj.plot_pka_time_series_for_selection(selection=args.selection, write_to_file=time_series_plot_name)
            stats_plot_name = os.path.join(output_folder, f"stats_{args.selection.join('_')}_{base_name}.png")
            pka_traj.plot_pka_statistic_for_selection(selection=args.selection, write_to_file=stats_plot_name)



if __name__ == '__main__':
    main()
