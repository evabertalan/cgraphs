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

    def plot_pka_time_series_for_selection(self, selection=None, write_to_file=None):
        resids = self._get_selected_res(selection)
        t = selection if selection else self.selection
        fig, ax = self._create_plot(title=t, xlabel='Time (ns)', ylabel=r'p$K_a$')
        for col in resids:
            ax.plot(self.pkatraj.pkas[col], label=col)
        ax.legend()

        if write_to_file:
            fig.savefig(write_to_file)

    def plot_pka_statistic_for_selection(self, selection=None, write_to_file=None):
        resids = self._get_selected_res(selection)
        t = selection if selection else self.selection
        fig, ax = self._create_plot(title=t, xlabel='Res ID', ylabel=r'p$K_a$')
        sns.boxplot(data=self.pkatraj.pkas[resids], ax=ax)

        if write_to_file:
            fig.savefig(write_to_file)

    def write_pka_to_external_data_file(self, write_to_file, selection=None):
        pass

    def _get_selected_res(self, selection):
        if selection:
            resids = self.u.select_atoms(selection).residues.resids
        else:
            resids = self.pkatraj.pkas.columns
        return resids

    def _create_plot(self, title='', xlabel='', ylabel='',figsize=(12, 8)):
        fig, ax = plt.subplots(figsize=figsize)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_title(title, fontsize=12)
        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
    #     ax.tick_params(axis='x', labelsize=plot_parameters['plot_tick_fontsize'])
    #     ax.tick_params(axis='y', labelsize=plot_parameters['plot_tick_fontsize'])
        return fig, ax
