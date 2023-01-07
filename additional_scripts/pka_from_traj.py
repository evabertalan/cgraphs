import MDAnalysis as mda
from propkatraj import PropkaTraj

class pKa_from_traj:
    def __init__(self, psf, dcd):
        self.psf = psf
        self.dcd = dcd
        self.u = mda.Universe(self.psf, self.dcd)
        print('Number of atoms in selection:', len(self.u.atoms))
        print('Number of residues in selection:', len(self.u.residues))
        print('Number of frames', len(self.u.trajectory))
        
    def compute_pka_for_traj(self, selection='protein'):
        print('Selection', selection)
        self.pkatraj = PropkaTraj(self.u, select=selection)
        self.pkatraj.run()
    
    def get_pka_for_frame(self, write_to_file=None):
        pkas = self.pkatraj.pkas
        if write_to_file:
            pkas.to_csv(write_to_file)
        return pkas
    
    def get_pka_statistic(self, write_to_file=None):
        pass
    
    def plot_pka_time_series_for_selection(self, selection, write_to_file=None):
        pass
    
    def plot_pka_time_series_for_residue(self, selection, write_to_file=None):
        pass
    
    def plot_pka_statistic_for_selection(self, selection, write_to_file=None):
        pass
