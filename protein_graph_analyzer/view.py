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

        # self.pdb_root_folder = '/Users/evabertalan/Documents/JSR/rhodopsin_crystal_structures/squid'
        # self.reference_pdb = '/Users/evabertalan/Documents/JSR/rhodopsin_crystal_structures/squid/2z73_sup.pdb'

        # self.psf_file = '/Users/evabertalan/Documents/protein_graph_analyzer/test_trajs/read_protein_membrane_7_opt_3_2x.psf'
        # self.dcd_files = ('/Users/evabertalan/Documents/protein_graph_analyzer/test_trajs/1-pbc.dcd','/Users/evabertalan/Documents/protein_graph_analyzer/test_trajs/9cis_optimized_last_20frames_pbc.dcd')
        # self._target_folder = '/Users/evabertalan/Documents/protein_graph_analyzer/test_trajs/'

        csa.csa_view(self)
        self.dcd_load_button= None
        self.graph_files = None
        self.DcdInfoFrame = None
        ta.ta_view(self)

# -------------------- crystal_strucutre_analyser_view ------------


    def _select_root_folder(self):
        self.pdb_root_folder = filedialog.askdirectory(initialdir = '../', parent=self.mainframe)
        self._configure_entry_field(self._input_folder, self.pdb_root_folder)

    def _select_reference_file(self):
        self.reference_pdb = filedialog.askopenfilename(initialdir = '../', filetypes=[('pdb', '.pdb')], parent=self.mainframe)
        self._configure_entry_field(self._input_pdb, self.reference_pdb)

    def _perform_parameter_analysis(self):
        sst = int(self.sequance_identity_threshold.get())/100
        # valudate sst
        self.w = WaterClusters(self.pdb_root_folder, reference_pdb=self.reference_pdb, sequance_identity_threshold=sst)
        self.w.fit_parameters()

    def _init_water_clusters(self):
        if not hasattr(self, 'w'):
            print('inintalaize warer clister')
            sst = int(self.sequance_identity_threshold.get())/100
            # valudate sst
            self.w = WaterClusters(self.pdb_root_folder, reference_pdb=self.reference_pdb, sequance_identity_threshold=sst)
        # self.w.evaluate_parameters(eps=float(self.eps.get())) #TEMPORARY FOR TESTING
        self.w.evaluate_parameters(eps=1.4)
        self.w.calculate_cluster_centers()
        self.w.write_cluster_center_coordinates()
        self.w.draw_clusters_centers_chimera()
        self.ref_coordinates = self.w.reference_coordinates
        tk.Label(self.waterClusterFrame, text='There are '+str(len(self.w.water_coordinates))+' water molecules in the '+str(len(self.w.superimposed_files))+' uperimposed files.\n The algorithm found '+str(self.w.n_clusters_)+' water clusters.').grid(row=5, column=0)
        self.w.logger.info('Water cluster calculation is completed\n'+'-'*20)

    def _init_pdb_conserved_graph_analysis(self, graph_type):
        sst = int(self.sequance_identity_threshold.get())/100
        ebb = False
        # ebb = not self.include_backbone_backbone.get()
        ieb = self.include_backbone_sidechain.get()
        if self.useWaterCoords.get(): _ref_coord = self.ref_coordinates
        else: _ref_coord=None
        c = ConservedGraph(self.pdb_root_folder, reference_pdb=self.reference_pdb, reference_coordinates=_ref_coord, sequance_identity_threshold=sst)
        if graph_type == 'water_wire': c.calculate_graphs(graph_type=graph_type, max_water=int(self.max_water.get()))
        else: c.calculate_graphs(graph_type=graph_type, exclude_backbone_backbone=ebb, include_backbone_sidechain=ieb)
        self._plot_conserved_graphs(c)

#--------------------- trajectory_analyser_view ------------

    def _select_target_folder(self):
        self._target_folder = filedialog.askdirectory(initialdir = '../',  parent=self.DcdWaterWireFrame)
        self._configure_entry_field(self._input_target, self._target_folder)

    def _select_psf_file(self):
        self.psf_file = filedialog.askopenfilename(initialdir = '../', title='Select protein structure file file', filetypes=[('psf', '.psf')], parent=self.DcdWaterWireFrame)
        self._configure_entry_field(self._input_psf, self.psf_file)

    def _select_dcd_files(self):
        self.dcd_files = filedialog.askopenfilenames(initialdir = '../', title='Select trajectory files', filetypes=[('dcd', '.dcd')],  parent=self.DcdWaterWireFrame)
        self._configure_entry_field(self._input_dcd, self.dcd_files)

    def _construct_sim_graphs(self):
        if self.DcdInfoFrame: self.DcdInfoFrame.destroy()
        p = ProteinGraphAnalyser(type_option='dcd', dcd_files=self.dcd_files, psf_file=self.psf_file, sim_name=self.sim_name.get(), target_folder=self._target_folder)
        p.calculate_graphs(graph_type='water_wire', max_water=int(self.sim_max_water.get()))
        self.DcdInfoFrame = tk.Frame(self.selectSimFrame)
        self.DcdInfoFrame .grid(row=11, column=1, columnspan=2, sticky="EW")
        tk.Label(self.DcdInfoFrame, text='Calculation completed', fg='green').grid(row=12, column=0, sticky='W')
        tk.Label(self.DcdInfoFrame, text='Now you can calculate the water wire network or costruct \ngraphs from other simuliation and then calculate the conserved network.').grid(row=13, column=0, sticky='W')


    def _init_dcd_conserved_graph_analysis(self):
        c_dcd = ConservedGraph(type_option='dcd', target_folder=self._target_folder)
        c_dcd._load_exisitng_graphs(graph_files=self.graph_files, graph_type='water_wire')

        if self.dcd_load_button: self.dcd_load_button.destroy()
        self.DcdOptionsFrame = tk.Frame(self.LoadGraphFrame)
        self.DcdOptionsFrame.grid(row=self.row+1, column=0, columnspan=2)
        self.conservation_threshold_dcd = tk.StringVar(value='90')
        tk.Label(self.DcdOptionsFrame, text='Conservation of H-bonding groups across structures (%)', anchor='w').grid(row=self.row+2, column=0, sticky='W')
        ttk.Spinbox(self.DcdOptionsFrame, textvariable=self.conservation_threshold_dcd, from_=1, to=100).grid(row=self.row+2, column=1, sticky="EW")
        self.min_occupancy = tk.StringVar(value='10')
        tk.Label(self.DcdOptionsFrame, text='Minimum H-bond occupancy (%)', anchor='w').grid(row=self.row+3, column=0, sticky='W')
        ttk.Spinbox(self.DcdOptionsFrame, textvariable=self.min_occupancy, from_=1, to=100).grid(row=self.row+3, column=1, sticky="EW")

        self.dcd_calc_button = tk.Button(self.LoadGraphFrame, text='Calculate conserved network', command=lambda:self._plot_conserved_graphs(c_dcd, int(self.min_occupancy.get())/100), width=self.button_width)
        self.dcd_calc_button.grid(self._create_big_button_grid(self.row+4))


    def _load_graph_files(self, row):
        if self.dcd_load_button:
            self.LoadGraphFrame.destroy()
            self.graph_files = None
        self.graph_files = filedialog.askopenfilenames(initialdir =self._target_folder+'/workfolder/.graph_objects/' , title='Select simulation graphs', filetypes=[('pickle', '.pickle')], parent=self.DcdWaterWireFrame)
        if self.graph_files:
            self.LoadGraphFrame = tk.Frame(self.DcdWaterWireFrame)
            self.LoadGraphFrame.grid(self._crate_frame_grid(row, columnspan=2))
            sh = tk.Scrollbar(self.LoadGraphFrame, orient='vertical')
            sh.grid(row=row+2, column=1)
            self.graph_text_field = tk.Text(self.LoadGraphFrame, height=4, width=80, yscrollcommand=sh.set)
            self.graph_text_field.grid(row=row+2,  column=0, sticky="EW")
            sh.config(command=self.graph_text_field.yview)
            tk.Label(self.LoadGraphFrame, text='Selected simulations for conserved network calculation: ', anchor='w').grid(row=row+1, column=0, sticky='W')

            for graph_file in self.graph_files:
                sim_name = graph_file.split('/')[-1].split('_water_wire')[0]
                self.graph_text_field.insert('end', sim_name+'\n')
                row += 1
            self.graph_text_field.configure(state='disabled')
            self.dcd_load_button = tk.Button(self.LoadGraphFrame, text='Load graphs', command=self._init_dcd_conserved_graph_analysis, width=self.button_width)
            self.dcd_load_button.grid(self._create_big_button_grid(row+3))
            self.row = row+4


#--------------------- COMMON ---------------------

    def _plot_conserved_graphs(self, c, occupancy=None):
        c.plot_graphs(label_nodes=True, occupancy=occupancy)
        c.plot_graphs(label_nodes=False, occupancy=occupancy)
        c.plot_linear_lenghts(occupancy=occupancy)
        cth = int(self.conservation_threshold.get())/100
        c.get_conserved_graph(conservation_threshold=cth, occupancy=occupancy)
        c.plot_conserved_graph(label_nodes=True)
        c.plot_conserved_graph(label_nodes=False)
        c.plot_difference(label_nodes=True)
        c.plot_difference(label_nodes=False)
        c.logger.info('Calculation completed\n'+'-'*20)

    def _configure_entry_field(self, field, value=None):
        field.configure(state='normal')
        if value: field.insert(0, str(value))
        else: field.delete(0, 'end')
        field.configure(state='disabled')

    def _add_horisontal_scroll(self, target, row=1, column=0):
        scroll = tk.Scrollbar(target, orient='horizontal')
        scroll.grid(row=row, column=column, sticky='EW')
        return scroll

    def _destroy_frame(self):
        self.mainframe.destroy()

    def _crate_frame_grid(self, row, columnspan=3):
        return {
            'row': row,
            'columnspan': columnspan,
            'sticky': 'EW',
            'padx': (5, 5),
            'pady': (5, 5)
        }

    def _create_big_button_grid(self, row, column=0):
        return {
            'row': row,
            'column': column,
            'padx': (5, 5),
            'pady': (5, 5),
            'sticky': "EW"
        }

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
