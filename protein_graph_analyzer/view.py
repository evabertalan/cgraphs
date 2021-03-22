import tkinter as tk
from tkinter import ttk
from . import crystal_strucutre_analyser_view as csa
from .waterclusters import WaterClusters
from .conservedgraph import ConservedGraph


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

        self.pdb_root_folder = '/Users/evabertalan/Documents/protein_graph_analyzer/workfolder/test_files_GlplG'
        self.reference_pdb = '/Users/evabertalan/Documents/protein_graph_analyzer/workfolder/2irv_aout.pdb'

    # self.pdb_root_folder = '/Users/evabertalan/Documents/protein_graph_analyzer/workfolder/test_files_GPCR'
    # self.reference_pdb = '/Users/evabertalan/Documents/protein_graph_analyzer/workfolder/test_files_GPCR/4eiy_opm.pdb'


        csa.csa_view(self)

    def _destroy_frame(self):
        self.mainframe.destroy()

    def _create_frame(self):
        tab_parnt = ttk.Notebook(self.master)
        self.mainframe = ttk.Frame(tab_parnt)
        self.dcdframe = ttk.Frame(tab_parnt)

        tab_parnt.add(self.mainframe, text='Crystal structure analysis')
        tab_parnt.add(self.dcdframe, text='MD trajectory analysis')
        tab_parnt.place(relx=0.5, rely=0.5, anchor=tk.CENTER)

    def _select_root_folder(self):
        # self.pdb_root_folder = filedialog.askdirectory(initialdir = "../")
        self._input_folder.configure(state='normal')
        self._input_folder.insert(0, str(self.pdb_root_folder))
        self._input_folder.configure(state='disabled')

    def _select_reference_file(self):
        # self.reference_pdb = filedialog.askopenfilename(initialdir = "../")
        self._input_pdb.configure(state='normal')
        self._input_pdb.insert(0, str(self.reference_pdb))
        self._input_pdb.configure(state='disabled')

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
        # c.logger.info('Calculation completed\n'+'-'*20)
        c.logger.info('Calculation completed\n'+'-'*20)

    def _add_horisontal_scroll(self, target, row=1, column=0):
        scroll = tk.Scrollbar(target, orient='horizontal')
        scroll.grid(row=row, column=column, sticky='EW')
        return scroll

    def _update_lable_text(self, text):
        self.completedText.set(text)


def start():
    root = tk.Tk()
    view = View(root)
    view.main_modal()
    root.mainloop()
