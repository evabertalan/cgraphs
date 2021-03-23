import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from . import crystal_strucutre_analyser_view as csa
from . import trajectory_analyser_view as ta
from .waterclusters import WaterClusters
from .conservedgraph import ConservedGraph
from .proteingraphanalyser import ProteinGraphAnalyser


class View:
    def __init__(self, master):
        self.master = master
        self.ipadx = 1
        self.ipady = 1
        self.padx = 1
        self.pady =1
        self.button_width = 1

    def main_modal(self):
        if hasattr(self, 'mainframe'):
            self._destroy_frame()

        self.master.title('Protein graph analyser')
        self.master.geometry('950x700')
        self._create_frame()

        # self.pdb_root_folder = '/Users/evabertalan/Documents/protein_graph_analyzer/workfolder/test_files_GlplG'
        # self.reference_pdb = '/Users/evabertalan/Documents/protein_graph_analyzer/workfolder/2irv_aout.pdb'

        self.pdb_root_folder = '/Users/evabertalan/Documents/JSR/rhodopsin_crystal_structures/squid'
        self.reference_pdb = '/Users/evabertalan/Documents/JSR/rhodopsin_crystal_structures/squid/2z73_sup.pdb'


        csa.csa_view(self)

        # self.psf_files = {}
        self.psf_files = ['/Users/evabertalan/Documents/protein_graph_analyzer/test_trajs/read_protein_membrane_7_opt_3_2x.psf']
        self.dcd_files = [('/Users/evabertalan/Documents/protein_graph_analyzer/test_trajs/1-pbc.dcd','/Users/evabertalan/Documents/protein_graph_analyzer/test_trajs/9cis_optimized_last_20frames_pbc.dcd')]
        self.sim_names = []
        self._target_folder = '/Users/evabertalan/Documents/protein_graph_analyzer/test_trajs/'
        ta.ta_view(self)

# -------------------- crystal_strucutre_analyser_view ------------


    def _select_root_folder(self):
        # self.pdb_root_folder = filedialog.askdirectory(initialdir = "../")
        self._configure_entry_field(self._input_folder, self.pdb_root_folder)

    def _select_reference_file(self):
        # self.reference_pdb = filedialog.askopenfilename(initialdir = "../")
        self._configure_entry_field(self._input_pdb, self.reference_pdb)

        # tk.Label(self.mainframe, anchor='w', text='All the generated files can be found in:\n'+self.pdb_root_folder+'/workfolder/\n\n Plots are located in:\n'+self.pdb_root_folder+'/workfolder/plots/', wraplength=520).grid(row=10, column=0, columnspan=3)

    def _perform_parameter_analysis(self):
        sst = int(self.sequance_identity_threshold.get())/100
        # valudate sst
        self.w = WaterClusters(self.pdb_root_folder, reference_pdb=self.reference_pdb, sequance_identity_threshold=sst)
        self.w.fit_parameters()
        # tk.Label(self.waterClusterFrame, text='See results of the parameter analysis are in:\n'+self.pdb_root_folder+'/workfolder/water_clusters/parameter_analysis').grid(row=1, sticky="EW")

    def _init_water_clusters(self):
        if not hasattr(self, 'w'):
            print('inintalaize warer clister')
            sst = int(self.sequance_identity_threshold.get())/100
            # valudate sst
            self.w = WaterClusters(self.pdb_root_folder, reference_pdb=self.reference_pdb, sequance_identity_threshold=sst)
        self.w.evaluate_parameters(eps=float(self.eps.get()))
        self.w.calculate_cluster_centers()
        self.w.write_cluster_center_coordinates()
        self.w.draw_clusters_centers_chimera()
        self.ref_coordinates = self.w.reference_coordinates
        tk.Label(self.waterClusterFrame, text='There are '+str(len(self.w.water_coordinates))+' water molecules in the '+str(len(self.w.superimposed_files))+' uperimposed files.\n The algorithm found '+str(self.w.n_clusters_)+' water clusters.').grid(row=5, column=0)
        self.w.logger.info('Water cluster calculation is completed\n'+'-'*20)

    def _init_conserved_graph_analysis(self, graph_type):
        sst = int(self.sequance_identity_threshold.get())/100
        self._update_lable_text('') #FIX THIS, not working
        self.completed.grid_forget()
        ebb = False
        # ebb = not self.include_backbone_backbone.get()
        ieb = self.include_backbone_sidechain.get()
        if self.useWaterCoords.get(): _ref_coord = self.ref_coordinates
        else: _ref_coord=None
        c = ConservedGraph(self.pdb_root_folder, reference_pdb=self.reference_pdb, reference_coordinates=_ref_coord, sequance_identity_threshold=sst)
        if graph_type == 'water_wire': c.calculate_graphs(graph_type=graph_type, max_water=int(self.max_water.get()))
        else: c.calculate_graphs(graph_type=graph_type, exclude_backbone_backbone=ebb, include_backbone_sidechain=ieb)
        c.plot_graphs(label_nodes=True)
        c.plot_graphs(label_nodes=False)
        c.plot_linear_lenghts()
        cth = int(self.conservation_threshold.get())/100
        c.get_conserved_graph(conservation_threshold=cth)
        c.plot_conserved_graph(label_nodes=True)
        c.plot_conserved_graph(label_nodes=False)
        c.plot_difference(label_nodes=True)
        c.plot_difference(label_nodes=False)
        self._update_lable_text('Calculation completed')
        self.completed.configure(fg='green')
        self.completed.grid()
        # c.logger.info('Calculation completed\n'+'-'*20)
        c.logger.info('Calculation completed\n'+'-'*20)

#--------------------- trajectory_analyser_view ------------

    def _select_psf_file(self):
        psf_file = filedialog.askopenfilename(initialdir = "../", title='Select protein structure file file', filetypes=[('psf', '.psf')])
        # self.psf_files.update( { name: psf_file  } )
        self.psf_files.append(psf_file)
        self._configure_entry_field(self._input_psf, psf_file)

    def _select_dcd_files(self):
        dcd_files = filedialog.askopenfilenames(initialdir = "../", title='Select trajectory files', filetypes=[('dcd', '.dcd')])
        # self.dcd_files.update( { name: dcd_files  } )
        self.dcd_files.append(dcd_files)
        self._configure_entry_field(self._input_dcd, dcd_files)

    def _select_target_folder(self):
        # self._target_folder = filedialog.askdirectory(initialdir = "../")
        self._configure_entry_field(self._input_target, self._target_folder)

    def _constract_sim_graphs(self):
        _sim_names = []
        for i in range(len(self.sim_names)):
            _sim_names.append(self.sim_names[i].get())
        # p = ProteinGraphAnalyser(type_option='dcd', dcd_files=self.dcd_files,psf_files=self.psf_files, sim_names=self.sim_names[i].get(), target_folder=self._target_folder)
        p = ProteinGraphAnalyser(type_option='dcd', dcd_files=self.dcd_files,psf_files=self.psf_files, sim_names=_sim_names, target_folder=self._target_folder)
        p.calculate_graphs(graph_type='water_wire', max_water=int(self.sim_max_water.get()))
#--------------------- COMMON ---------------------

    def _configure_entry_field(self, field, value):
        field.configure(state='normal')
        field.insert(0, str(value))
        field.configure(state='disabled')

    def _add_horisontal_scroll(self, target, row=1, column=0):
        scroll = tk.Scrollbar(target, orient='horizontal')
        scroll.grid(row=row, column=column, sticky='EW')
        return scroll

    def _update_lable_text(self, text):
        self.completedText.set(text)

    def _destroy_frame(self):
        self.mainframe.destroy()

    def _create_frame(self):
        tab_parnt = ttk.Notebook(self.master)
        self.mainframe = ttk.Frame(tab_parnt)
        self.dcdframe = ttk.Frame(tab_parnt)

        tab_parnt.add(self.mainframe, text='Crystal structure analysis')
        tab_parnt.add(self.dcdframe, text='MD trajectory analysis')
        tab_parnt.place(relx=0.5, rely=0.5, anchor=tk.CENTER)


def start():
    root = tk.Tk()
    view = View(root)
    view.main_modal()
    root.mainloop()
