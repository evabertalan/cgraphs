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
        # self.ifnum_cmd = (self.master.register(lambda x, y :self.VaidateNum('%S', '%P', x, y)), '%S', '%P')
        self.ifnum_cmd = self.master.register(self.VaidateNum)


        self.pdb_root_folder= '/Users/evabertalan/Documents/c_test_files/test_new_features_GPCR/'
        self.reference_pdb='/Users/evabertalan/Documents/c_test_files/4eiy_opm.pdb'

    def main_modal(self):
        if hasattr(self, 'mainframe'):
            self._destroy_frame()

        self.master.title('C-Graphs - Protein Conserved Graph Analyser')
        self.master.geometry('900x750')
        self._create_frame()

        csa.csa_view(self)
        self.dcd_load_button= None
        self.graph_files = None
        self.DcdInfoFrame = None
        self.is_linear_lenght_plot_dcd = tk.BooleanVar()
        self.is_induvidual_graph_dcd = tk.BooleanVar()
        self.is_difference_graph_dcd = tk.BooleanVar()
        ta.ta_view(self)

# -------------------- crystal_strucutre_analyser_view ------------


    def _select_root_folder(self):
        self.pdb_root_folder = filedialog.askdirectory(parent=self.mainframe)
        self._configure_entry_field(self._input_folder, self.pdb_root_folder)

    def _select_reference_file(self):
        self.reference_pdb = filedialog.askopenfilename(filetypes=[('pdb', '.pdb')], parent=self.mainframe)
        self._configure_entry_field(self._input_pdb, self.reference_pdb)

    def _perform_parameter_analysis(self):
        sst = int(self.sequance_identity_threshold.get())/100
        self.w = WaterClusters(self.pdb_root_folder, reference_pdb=self.reference_pdb, sequance_identity_threshold=sst)
        self.w.fit_parameters()

    def _init_water_clusters(self):
        sst = int(self.sequance_identity_threshold.get())/100
        self.w = WaterClusters(self.pdb_root_folder, reference_pdb=self.reference_pdb, sequance_identity_threshold=sst)
        # self.w.evaluate_parameters(eps=float(self.eps.get())) #TEMPORARY FOR TESTING
        self.w.evaluate_parameters(eps=1.4)
        self.w.calculate_cluster_centers()
        self.w.write_cluster_center_coordinates()
        self.w.draw_clusters_centers_chimera()
        self.ref_coordinates = self.w.reference_coordinates
        tk.Label(self.waterClusterFrame, text=' There are '+str(len(self.w.water_coordinates))+' water molecules in the '+str(len(self.w.superimposed_files))+' superimposed files.\n The algorithm found '+str(self.w.n_clusters_)+' water clusters.', anchor='w', justify=tk.LEFT).grid(row=5, column=0, sticky='w')
        self.w.logger.info('Water cluster calculation is completed\n'+'-'*20)

    def _init_pdb_conserved_graph_analysis(self, graph_type):
        sst = int(self.sequance_identity_threshold.get())/100
        ebb = True
        # ebb = not self.include_backbone_backbone.get()
        ieb = self.include_backbone_sidechain.get()
        if self.useWaterCoords.get(): _ref_coord = self.ref_coordinates
        else: _ref_coord=None
        c = ConservedGraph(self.pdb_root_folder, reference_pdb=self.reference_pdb, reference_coordinates=_ref_coord, sequance_identity_threshold=sst)
        if graph_type == 'water_wire': c.calculate_graphs(graph_type=graph_type, max_water=int(self.max_water.get()), distance=float(self.c_distance.get()), cut_angle=float(self.c_cut_angle.get()), selection=self.selection_string.get())
        else: c.calculate_graphs(graph_type=graph_type, exclude_backbone_backbone=ebb, include_backbone_sidechain=ieb, include_waters=self.include_waters_hbond.get(), distance=float(self.c_distance.get()), cut_angle=float(self.c_cut_angle.get()), selection=self.selection_string.get())
        self._plot_conserved_graphs(c, self.is_linear_lenght_plot.get(), self.is_induvidual_graph.get(), self.is_difference_graph.get(), cth=int(self.conservation_threshold.get())/100)

#--------------------- trajectory_analyser_view ------------

    def _select_target_folder(self):
        self._target_folder = filedialog.askdirectory(parent=self.DcdWaterWireFrame)
        self._configure_entry_field(self._input_target, self._target_folder)

    def _select_psf_file(self):
        self.psf_file = filedialog.askopenfilename(title='Select protein structure file file', filetypes=[('psf', '.psf')], parent=self.DcdWaterWireFrame)
        self._configure_entry_field(self._input_psf, self.psf_file)

    def _select_dcd_files(self):
        self.dcd_files = filedialog.askopenfilenames(title='Select trajectory files', filetypes=[('dcd', '.dcd')],  parent=self.DcdWaterWireFrame)
        self._configure_entry_field(self._input_dcd, self.dcd_files)

    def _construct_sim_graphs(self):
        if not hasattr(self, '_target_folder'): print('WARNING: Please select a folder to Save results to!')
        if self.DcdInfoFrame: self.DcdInfoFrame.destroy()
        p = ProteinGraphAnalyser(type_option='dcd', dcd_files=self.dcd_files, psf_file=self.psf_file, sim_name=self.sim_name.get(), target_folder=self._target_folder)
        p.calculate_graphs(graph_type='water_wire', max_water=int(self.sim_max_water.get()), distance=float(self.sim_distance.get()), cut_angle=float(self.sim_cut_angle.get()), selection=self.sim_selection_string.get())
        self.DcdInfoFrame = tk.Frame(self.selectSimFrame)
        self.DcdInfoFrame.grid(row=12, column=1, columnspan=2, sticky="EW")
        tk.Label(self.DcdInfoFrame, text='Calculation completed for '+self.sim_name.get(), fg='green', anchor='w').grid(row=13, column=0, sticky='W')


    def _init_dcd_conserved_graph_analysis(self):
        c_dcd = ConservedGraph(type_option='dcd', target_folder=self._target_folder)
        c_dcd._load_exisitng_graphs(graph_files=self.graph_files, graph_type='water_wire')

        if self.dcd_load_button: self.dcd_load_button.destroy()
        self.DcdOptionsFrame = tk.Frame(self.LoadGraphFrame)
        self.DcdOptionsFrame.grid(row=self.row+1, column=0, columnspan=2)
        self.conservation_threshold_dcd = tk.StringVar(value='90')
        tk.Label(self.DcdOptionsFrame, text='Conservation of H-bonding groups across structures (%)', anchor='w').grid(row=self.row+2, column=0, sticky='W')
        ttk.Spinbox(self.DcdOptionsFrame, textvariable=self.conservation_threshold_dcd, from_=1, to=100, validate="key", validatecommand=(self.ifnum_cmd, '%S', '%P', 0, 100)).grid(row=self.row+2, column=1, sticky="EW")
        self.min_occupancy = tk.StringVar(value='10')
        tk.Label(self.DcdOptionsFrame, text='Minimum H-bond occupancy (%)', anchor='w').grid(row=self.row+3, column=0, sticky='W')
        ttk.Spinbox(self.DcdOptionsFrame, textvariable=self.min_occupancy, from_=1, to=100, validate="key", validatecommand=(self.ifnum_cmd, '%S', '%P', 0, 100)).grid(row=self.row+3, column=1, sticky="EW")

        tk.Checkbutton(self.DcdOptionsFrame, text='Plot network for each structure', variable=self.is_induvidual_graph_dcd, anchor="w").grid(self._create_big_button_grid(self.row+4))
        tk.Checkbutton(self.DcdOptionsFrame, text='Plot difference graph for each structure', variable=self.is_difference_graph_dcd, anchor="w").grid(self._create_big_button_grid(self.row+5))
        tk.Checkbutton(self.DcdOptionsFrame, text='Plot linear lengths of continuous networks for each structure', variable=self.is_linear_lenght_plot_dcd, anchor="w").grid(self._create_big_button_grid(self.row+6))

        self.dcd_calc_button = tk.Button(self.LoadGraphFrame, text='Calculate conserved network', command=lambda:self._plot_conserved_graphs(c_dcd, self.is_linear_lenght_plot_dcd.get(), self.is_induvidual_graph_dcd.get(), self.is_difference_graph_dcd.get(), cth=int(self.conservation_threshold_dcd.get())/100, occupancy=int(self.min_occupancy.get())/100), width=self.button_width)
        self.dcd_calc_button.grid(self._create_big_button_grid(self.row+7))


    def _load_graph_files(self, row):
        if not hasattr(self, '_target_folder'): print('WARNING: Please select a folder to Save results to!')
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

            prev_water_number = None
            for graph_file in self.graph_files:
                _split = graph_file.split('/')[-1].split('_water_wire')[0]
                sim_name = _split[0:-2]
                water_number = _split[-1]
                if prev_water_number != None and prev_water_number != water_number:
                    print('WARNING: Calculating conserved graphs is only possible from networks, in which the maximum number of water molecules allowed in the bridge is the same. Please select graphs to compare with the same allowed maximum waters.')
                    self.LoadGraphFrame.destroy()
                    self.graph_files = None
                    return
                prev_water_number = water_number
                self.graph_text_field.insert('end', sim_name+'\n')
                row += 1
            self.graph_text_field.configure(state='disabled')
            self.dcd_load_button = tk.Button(self.LoadGraphFrame, text='Load graphs', command=self._init_dcd_conserved_graph_analysis, width=self.button_width)
            self.dcd_load_button.grid(self._create_big_button_grid(row+3))
            self.row = row+4


#--------------------- COMMON ---------------------

    def _plot_conserved_graphs(self, c, plot_linear_length, plot_induvidual_graph, plot_difference_graph, cth=0.9, occupancy=None):
        c.get_conserved_graph(conservation_threshold=cth, occupancy=occupancy)
        c.plot_conserved_graph(label_nodes=True)
        c.plot_conserved_graph(label_nodes=False)
        if plot_linear_length:
            c.plot_linear_lenghts(occupancy=occupancy)
        if plot_induvidual_graph:
            c.plot_graphs(label_nodes=True, occupancy=occupancy)
            c.plot_graphs(label_nodes=False, occupancy=occupancy)
        if plot_difference_graph:
            c.plot_difference(label_nodes=True)
            c.plot_difference(label_nodes=False)
        c.logger.info('Calculation completed\n'+'-'*20)

    def _configure_entry_field(self, field, value=None):
        field.configure(state='normal')
        field.delete(0, 'end')
        if value: field.insert(0, str(value))
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
        style = ttk.Style()
        gray = '#f5f5f5'
        style.theme_create( 'style', settings={
            '.': {'configure': {'background': 'white', 'relief': 'flat', 'takefocus':'false'}},
            'TNotebook': {'configure': {'tabmargins': [2, 5, 0, 0] } },
            'TFrame': {'configure': {'relief': 'flat', 'padding': [30,8,30,10]}},
            'TLabelframe': {'configure': {'relief': 'flat', 'padding': [30,8,30,10]}},
            'TLabelframe.Label': {'configure': {'font': ('Helvetica', 13, 'bold')}},
            'TNotebook.Tab': {
                    'configure': {'padding': [8, 4], 'background': gray },
                    'map': {'background': [('selected', 'white')]}}
        } )
        style.theme_use('style')

        tab_parnt = ttk.Notebook(self.master)
        self.mainframe = ttk.Frame(tab_parnt)
        self.dcdframe = ttk.Frame(tab_parnt)

        self.dcdframe.grid_columnconfigure(0, weight=1)

        tab_parnt.add(self.mainframe, text='Crystal structure analysis')
        tab_parnt.add(self.dcdframe, text='MD trajectory analysis')
        tab_parnt.place(relx=0.5, rely=0.5, anchor=tk.CENTER)

    def VaidateNum(self, S, P, _min, _max):
        '''S If the call was due to an insertion or deletion, this argument will be the text being inserted or deleted. '''
        ''' P The value that the text will have if the change is allowed. '''
        valid = (S.isdigit() or S == '.') and P != '' and float(P)>=float(_min) and float(P)<=float(_max)
        if not valid:
            print('Please select a number between '+ str(_min)+' and '+str(_max))
            self.master.bell()
        return valid

def start():
    root = tk.Tk()
    view = View(root)
    view.main_modal()
    root.mainloop()
