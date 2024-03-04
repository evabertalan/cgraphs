import MDAnalysis as mda
from propkatraj import PropkaTraj
import seaborn as sns
import matplotlib.pyplot as plt

class pKa_from_traj:
    def __init__(self, psf, dcd):
        self.psf = psf
        self.dcd = dcd
        self.u = mda.Universe(self.psf, self.dcd)

        print('Number of atoms in selection:', len(self.u.atoms))
        print('Number of residues in selection:', len(self.u.residues))
        print('Number of frames', len(self.u.trajectory))


    def compute_pka_for_traj(self, selection='protein'):
        self.selection = selection
        print('Selection', self.selection)
        self.pkatraj = PropkaTraj(self.u, select=selection)
        self.pkatraj.run()


    def get_pka_for_frame(self, write_to_file=None):
        pkas = self.pkatraj.pkas
        if write_to_file:
            pkas.to_csv(write_to_file)
        return pkas


    def get_pka_statistic(self, selection=None, write_to_file=None):
        resids = self._get_selected_res(selection)
        stats = self.pkatraj.pkas[resids].describe()
        if write_to_file:
             stats.to_csv(write_to_file)
        return stats


    def plot_pka_time_series_for_selection(self, selection=None, write_to_file=None, figsize=None):
        resids = self._get_selected_res(selection)
        t = selection if selection else self.selection
        fig, ax = self._create_plot(title=t, xlabel='Time (ns)', ylabel=r'p$K_a$')
        for col in resids:
            ax.plot(self.pkatraj.pkas[col], label=col)
        ax.legend()

        if write_to_file:
            fig.savefig(write_to_file)


    def plot_pka_statistic_for_selection(self, selection=None, write_to_file=None, figsize=None):
        resids = self._get_selected_res(selection)
        t = selection if selection else self.selection
        fig, ax = self._create_plot(title=t, xlabel='Res ID', ylabel=r'p$K_a$', figsize=figsize)
        sns.boxplot(data=self.pkatraj.pkas[resids], ax=ax)

        if write_to_file:
            fig.savefig(write_to_file)


    def write_pka_to_external_data_file(self, write_to_file, selection=None):
        assert write_to_file.endswith('_data.txt'), 'Name of external data file should end with _data.txt/n The first part of the name should match the PDB file analized e.g: 6b73_data.txt'

        write_to_file.endswith('_data.txt')
        stats = self.get_pka_statistic(selection)

        res_name_map = {}
        prot = self.u.select_atoms('protein')
        for res in prot.residues:
            res_id = res.resid
            res_name = res.resname
            seg_id = res.segid
            res_name_map.update({f'{res_id}': {'res_name': res_name, 'seg_id': seg_id} })

        with open(write_to_file, 'w') as f:
            for res_id, avg_pka in zip(stats.columns, stats.iloc[1]):
                line = ' '.join([res_name_map[str(res_id)]['res_name'], str(res_id), res_name_map[str(res_id)]['seg_id'], str(avg_pka)])
                f.write(line + '\n')


    def write_last_frame_as_pdb(self, write_to_file):
        assert write_to_file.endswith('.pdb'), 'File name has to end to .pdb\n e.g.: 6b73.pdb'

        write_to_file.endswith('.pdb')
        self.u.trajectory[-1]
        with mda.Writer(write_to_file, multiframe=True, bonds=None, n_atoms=self.u.atoms.n_atoms) as PDB:
            PDB.write(self.u.atoms)


    def _get_selected_res(self, selection):
        if selection:
            resids = self.u.select_atoms(selection).residues.resids
        else:
            resids = self.pkatraj.pkas.columns
        return resids


    def _create_plot(self, title='', xlabel='', ylabel='', figsize=None):
        figsize = figsize if figsize is not None else (12, 8)
        fig, ax = plt.subplots(figsize=figsize)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_title(title, fontsize=12)
        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
    #     ax.tick_params(axis='x', labelsize=plot_parameters['plot_tick_fontsize'])
    #     ax.tick_params(axis='y', labelsize=plot_parameters['plot_tick_fontsize'])
        return fig, ax
