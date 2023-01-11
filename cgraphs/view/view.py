import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from tkinter import colorchooser
from . import crystal_strucutre_analyser_view as csa
from . import trajectory_analyser_view as ta
from . import compare_2_view as comp
from . import plot_settings_view as ps
from ..waterclusters import WaterClusters
from ..conservedgraph import ConservedGraph
from ..proteingraphanalyser import ProteinGraphAnalyser
from ..comparetwo import CompareTwo
import re

class popupWindow(object):
    def __init__(self, master):
        self.top = tk.Toplevel(master, bg='white')
        self.top.geometry('650x200')
        self.top.columnconfigure(1, weight=1)

    def custom_selection_sting(self, selection_entry, selected_donors, selected_acceptors):
        tk.Label(self.top, text='Customize selection string and additional donors and acceptors', bg='white', fg='black', pady=4).grid(row=0, column=0, sticky="EW", columnspan=3)
        tk.Label(self.top, text='Selection string:', bg='white', fg='black').grid(row=1, column=0, sticky="W")
        scroll = tk.Scrollbar(self.top, orient='horizontal')
        self.sel_string = tk.Entry(self.top, xscrollcommand=scroll.set, bg='white', fg='black', highlightbackground='white', insertbackground='black')
        self.sel_string.grid(row=1, column=1, sticky="EW", columnspan=2)
        self.sel_string.insert(0, str(selection_entry.get()))
        scroll.grid(row=2, column=1, sticky='EW', columnspan=2)
        scroll.configure(command=self.sel_string.xview, bg='white')

        tk.Label(self.top, text='List of additional donors:', bg='white', fg='black').grid(row=4, column=0, sticky="W")
        self.sel_donors = tk.Entry(self.top, bg='white', fg='black', highlightbackground='white', insertbackground='black')
        self.sel_donors.grid(row=4, column=1, sticky="EW",  columnspan=2)
        self.sel_donors.insert(0, str(selected_donors.get()))

        tk.Label(self.top, text='List of additional acceptors:', bg='white', fg='black').grid(row=5, column=0, sticky="W")
        self.sel_acceptors = tk.Entry(self.top,  bg='white', fg='black', highlightbackground='white', insertbackground='black')
        self.sel_acceptors.grid(row=5, column=1, sticky="EW",  columnspan=2)
        self.sel_acceptors.insert(0, str(selected_acceptors.get()))

        ok_button = tk.Button(self.top, text='Ok', command=self.cleanup, highlightbackground='white', bg='white', fg='black')
        ok_button.grid(row=6, column=0, sticky="EW", columnspan=3)
        self.top.protocol("WM_DELETE_WINDOW", self.cleanup)

    def node_color_selection(self, selected_nodes_for_color, selected_color):
        tk.Label(self.top, text='Specify residues for node coloring', bg='white', fg='black', pady=4).grid(row=0, column=0, sticky="EW", columnspan=3)
        tk.Label(self.top, text='Selection string:', bg='white', fg='black').grid(row=1, column=0, sticky="W")
        scroll = tk.Scrollbar(self.top, orient='horizontal')

        self.nc_sel_string = tk.Entry(self.top, xscrollcommand=scroll.set, bg='white', fg='black', highlightbackground='white', insertbackground='black')
        self.nc_sel_string.grid(row=1, column=1, sticky="EW", columnspan=2)
        self.nc_sel_string.insert(0, str(selected_nodes_for_color.get()))
        scroll.grid(row=2, column=1, sticky='EW', columnspan=2)
        scroll.configure(command=self.nc_sel_string.xview, bg='white')

        color_options = ('plasma', 'viridis', 'RdBu', 'Blues', 'Greens', 'Reds', 'YlGn', 'PiYG', 'coolwarm', 'bwr', 'RdYlBu', 'PRGn')
        self.nc_color_map = tk.StringVar()
        self.nc_color_map.set(selected_color.get())
        tk.Label(self.top, text='Color map:', bg='white', fg='black').grid(row=2, column=0, sticky="W")
        tk.OptionMenu(self.top, self.nc_color_map, *color_options).grid(row=2, column=1, sticky="EW", columnspan=2)

        ok_button = tk.Button(self.top, text='Ok', command=self.close, highlightbackground='white', bg='white', fg='black')
        ok_button.grid(row=6, column=0, sticky="EW", columnspan=3)
        self.top.protocol("WM_DELETE_WINDOW", self.close)

    def cleanup(self):
        self._sel_string=self.sel_string.get()
        self._sel_donors=self.sel_donors.get()
        self._sel_acceptors=self.sel_acceptors.get()
        self.top.destroy()

    def close(self):
        self._nc_sel_string = self.nc_sel_string.get()
        self._nc_sel_color = self.nc_color_map.get()
        self.top.destroy()



class View:
    def __init__(self, master):
        self.master = master
        self.padx = 1
        self.pady = 1
        self.button_width = 1
        self.ifnum_cmd = self.master.register(self.VaidateNum)
        self.gray = '#f5f5f5'
        self.plot_parameters = {}

    def main_modal(self):
        if hasattr(self, 'mainframe'):
            self._destroy_frame()

        self.master.title('C-Graphs - Protein Conserved Graph Analyser')
        self.master.geometry('900x780')
        self._create_frame()

        csa.csa_view(self)
        self.dcd_load_button= None
        self.graph_files = None
        self.DcdInfoFrame = None
        self.water_wire_comp_frame_dcd = None
        self.is_linear_lenght_plot_dcd = tk.BooleanVar()
        self.is_induvidual_graph_dcd = tk.BooleanVar()
        self.is_difference_graph_dcd = tk.BooleanVar()
        ta.ta_view(self)
        comp.compare_view(self)
        ps.plot_settings(self)

# -------------------- crystal_strucutre_analyser_view ------------
    def _select_pdb_root_folder(self, field):
        self.pdb_root_folder = filedialog.askdirectory(parent=self.mainframe)
        self._configure_entry_field(field, self.pdb_root_folder)

    def _select_reference_pdb(self, field):
        self.reference_pdb = filedialog.askopenfilename(filetypes=[('pdb', '.pdb')], parent=self.mainframe)
        self._configure_entry_field(field, self.reference_pdb)

    def _perform_parameter_analysis(self):
        sst = int(self.sequance_identity_threshold.get())/100
        self.w = WaterClusters(self.pdb_root_folder, reference_pdb=self.reference_pdb, sequance_identity_threshold=sst,superimposition_threshold=float(self.superimposition_threshold.get()))
        if self.w.valid_structures_for_clustering:
            self.w.fit_parameters()

    def _init_water_clusters(self):
        if not hasattr(self, 'pdb_root_folder') or self.pdb_root_folder is None:
            print('WARNING: Please select a PDB folder and reference file!')
        else:
            self.w = WaterClusters(
                    self.pdb_root_folder,
                    reference_pdb=self.reference_pdb,
                    sequance_identity_threshold=int(self.sequance_identity_threshold.get())/100,
                    superimposition_threshold=float(self.superimposition_threshold.get()))
            if self.w.valid_structures_for_clustering:
                self.w.evaluate_parameters(eps=float(self.eps.get()))
                # self.w.evaluate_parameters(eps=1.4)
                self.w.calculate_cluster_centers()
                self.w.write_cluster_center_coordinates()
                self.w.draw_clusters_centers_chimera()
                self.ref_coordinates = self.w.reference_coordinates
                tk.Label(self.waterClusterFrame, text=' There are '+str(len(self.w.water_coordinates))+' water molecules in the '+str(len(self.w.superimposed_files))+' superimposed files.\n The algorithm found '+str(self.w.n_clusters_)+' water clusters.', anchor='w', justify=tk.LEFT, fg='green', bg='white').grid(row=5, column=0, sticky='w')
                self.w.logger.info('Water cluster calculation is completed\n'+'-'*20)

    def _init_pdb_conserved_graph_analysis(self, graph_type):
        if not hasattr(self, 'pdb_root_folder') or self.pdb_root_folder is None:
            print('WARNING: Please select a PDB folder and reference file!')
        else:
            sst = int(self.sequance_identity_threshold.get())/100
            ebb = True
            # ebb = not self.include_backbone_backbone.get()
            ieb = self.include_backbone_sidechain.get()
            if self.useWaterCoords.get() and (not hasattr(self, 'ref_coordinates') or self.ref_coordinates is None):
                print('WARNING: There are no water cluster coordinates!')
            else:
                additional_donors = self._create_list_from_sting(self.selected_donors_pdb.get())
                additional_acceptors = self._create_list_from_sting(self.selected_acceptors_pdb.get())
                if self.useWaterCoords.get():
                    _ref_coord = self.ref_coordinates
                    eps =float(self.eps.get())
                else:
                    _ref_coord=None
                    eps = 1.5
                c = ConservedGraph(self.pdb_root_folder, reference_pdb=self.reference_pdb, reference_coordinates=_ref_coord, sequance_identity_threshold=sst, plot_parameters=self.plot_parameters)
                if graph_type == 'water_wire': c.calculate_graphs(graph_type=graph_type, max_water=int(self.max_water.get()), distance=float(self.c_distance.get()), cut_angle=float(self.c_cut_angle.get()), check_angle=self.c_use_angle.get(), selection=self.selection_string.get(), additional_donors=additional_donors, additional_acceptors=additional_acceptors)
                else: c.calculate_graphs(graph_type=graph_type, exclude_backbone_backbone=ebb, include_backbone_sidechain=ieb, include_waters=self.include_waters_hbond.get(), distance=float(self.c_distance.get()), cut_angle=float(self.c_cut_angle.get()), check_angle=self.c_use_angle.get(), selection=self.selection_string.get(),  additional_donors=additional_donors, additional_acceptors=additional_acceptors)
                self._plot_conserved_graphs(c, self.is_linear_lenght_plot.get(), self.is_induvidual_graph.get(), self.is_difference_graph.get(), self.is_bond_distance_plot.get(), cth=int(self.conservation_threshold.get())/100, eps=eps)

#--------------------- trajectory_analyser_view ------------

    def _select_dcd_workfolder(self, field):
        self._target_folder = filedialog.askdirectory(parent=self.dcdframe)
        self._configure_entry_field(field, self._target_folder)

    def _select_psf_file(self):
        self.psf_file = filedialog.askopenfilename(title='Select protein structure file file', filetypes=[('psf', '.psf'), ('pdb', '.pdb'), ('gro', '.gro')], parent=self.DcdWaterWireFrame)
        self._configure_entry_field(self._input_psf, self.psf_file)

    def _select_dcd_files(self):
        self.dcd_files = filedialog.askopenfilenames(title='Select trajectory files', filetypes=[('dcd', '.dcd'), ('xtc', '.xtc'), ('trr', '.trr')],  parent=self.DcdWaterWireFrame)
        self._configure_entry_field(self._input_dcd, self.dcd_files)

    def _construct_sim_graphs(self):
        if not hasattr(self, '_target_folder') or self._target_folder is None:
            print('WARNING: Please select the location of the workfolder!')
        else:
            if self.DcdInfoFrame: self.DcdInfoFrame.destroy()
            additional_donors = self._create_list_from_sting(self.sim_selected_donors_pdb.get())
            additional_acceptors = self._create_list_from_sting(self.sim_selected_acceptors_pdb.get())
            p = ProteinGraphAnalyser(type_option='dcd', dcd_files=[self.dcd_files], psf_files=[self.psf_file], sim_names=[self.sim_name.get()], target_folder=self._target_folder, plot_parameters=self.plot_parameters)
            p.calculate_graphs(graph_type='water_wire', max_water=int(self.sim_max_water.get()), distance=float(self.sim_distance.get()), cut_angle=float(self.sim_cut_angle.get()), check_angle=self.sim_use_angle.get(), selection=self.sim_selection_string.get(), additional_donors=additional_donors, additional_acceptors=additional_acceptors)
            self.DcdInfoFrame = tk.Frame(self.selectSimFrame, bg='white')
            self.DcdInfoFrame.grid(row=12, column=1, columnspan=2, sticky="EW")
            tk.Label(self.DcdInfoFrame, text='Calculation completed for '+self.sim_name.get(), fg='green', anchor='w', bg='white').grid(row=13, column=0, sticky='W')


    def _init_dcd_conserved_graph_analysis(self):
        if not hasattr(self, '_target_folder') or self._target_folder is None:
            print('WARNING: Please select the location of the workfolder!')
        else:
            c_dcd = ConservedGraph(type_option='dcd', target_folder=self._target_folder, plot_parameters=self.plot_parameters)
            c_dcd._load_exisitng_graphs(graph_files=self.graph_files, graph_type='water_wire', selection=self.sim_selection_string.get())

            if self.dcd_load_button: self.dcd_load_button.destroy()
            self.DcdOptionsFrame = tk.Frame(self.LoadGraphFrame, bg='white')
            self.DcdOptionsFrame.grid(row=self.row+1, column=0, columnspan=2)
            self.conservation_threshold_dcd = tk.DoubleVar(value=90)
            tk.Label(self.DcdOptionsFrame, text='Conservation of H-bonding groups across structures (%)', anchor='w', bg='white', fg='black').grid(row=self.row+2, column=0, sticky='W')
            ttk.Spinbox(self.DcdOptionsFrame, textvariable=self.conservation_threshold_dcd, from_=1, to=100, width=5,validate="key", validatecommand=(self.ifnum_cmd, '%S', '%P', 0, 100)).grid(row=self.row+2, column=1, sticky="EW")
            self.min_occupancy = tk.DoubleVar(value=10)
            tk.Label(self.DcdOptionsFrame, text='Minimum H-bond occupancy (%)', anchor='w', bg='white', fg='black').grid(row=self.row+3, column=0, sticky='W')
            ttk.Spinbox(self.DcdOptionsFrame, textvariable=self.min_occupancy, from_=1, to=100, width=5, validate="key", validatecommand=(self.ifnum_cmd, '%S', '%P', 0, 100)).grid(row=self.row+3, column=1, sticky="EW")

            tk.Label(self.DcdOptionsFrame, text='Plot for each structure:', anchor="w", bg='white', fg='black').grid(row=self.row+4, column=0, sticky='W')
            each_plots_dcd = tk.Frame(self.DcdOptionsFrame, bg='white')
            each_plots_dcd.grid(row=self.row+4, column=1, columnspan=3, sticky="EW")

            tk.Checkbutton(each_plots_dcd, text='Individual network    ', variable=self.is_induvidual_graph_dcd, anchor="w", bg='white', fg='black').grid(row=self.row+4, column=1)
            tk.Checkbutton(each_plots_dcd, text='Difference graph    ', variable=self.is_difference_graph_dcd, anchor="w", bg='white', fg='black').grid(row=self.row+4, column=2)
            tk.Checkbutton(each_plots_dcd, text='Linear lengths', variable=self.is_linear_lenght_plot_dcd, anchor="w", bg='white', fg='black').grid(row=self.row+4, column=3)

            self.dcd_calc_button = tk.Button(self.LoadGraphFrame, text='Calculate conserved network', command=lambda:self._plot_conserved_graphs(c_dcd, self.is_linear_lenght_plot_dcd.get(), self.is_induvidual_graph_dcd.get(), self.is_difference_graph_dcd.get(), False, cth=int(self.conservation_threshold_dcd.get())/100, occupancy=int(self.min_occupancy.get())/100), width=self.button_width, bg='white', fg='black', highlightbackground='white')
            self.dcd_calc_button.grid(self._create_big_button_grid(self.row+7))


    def _load_graph_files(self, row):
        if not hasattr(self, '_target_folder') or self._target_folder is None:
            print('WARNING: Please select the location of the workfolder!')
        else:
            if self.dcd_load_button:
                self.LoadGraphFrame.destroy()
                self.graph_files = None
            self.graph_files = filedialog.askopenfilenames(initialdir =self._target_folder+'/workfolder/graph_objects/', title='Select simulation graphs', filetypes=[('pickle', '.pickle')], parent=self.DcdWaterWireFrame)
            if self.graph_files:
                self.LoadGraphFrame = tk.Frame(self.DcdWaterWireFrame, bg='white')
                self.LoadGraphFrame.grid(self._crate_frame_grid(row, columnspan=2), sticky="EW")
                self.LoadGraphFrame.columnconfigure(0, weight=1)
                tk.Label(self.LoadGraphFrame, text='Selected simulations for conserved network calculation: ', anchor='w', bg='white', fg='black').grid(row=row+1, column=0, sticky='EW')
                sh = tk.Scrollbar(self.LoadGraphFrame, orient='vertical', bg='white')
                sh.grid(row=row+2, column=1)
                self.graph_text_field = tk.Text(self.LoadGraphFrame, height=4, yscrollcommand=sh.set, bg='white', fg='black', highlightbackground='white')
                self.graph_text_field.grid(row=row+2,  column=0, sticky="EW")
                sh.config(command=self.graph_text_field.yview)

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
                self.dcd_load_button = tk.Button(self.LoadGraphFrame, text='Load graphs', command=self._init_dcd_conserved_graph_analysis, width=self.button_width, bg='white', fg='black', highlightbackground='white')
                self.dcd_load_button.grid(self._create_big_button_grid(row+3))
                self.row = row+4

#--------------------- COMPARE 2 STRUCTURES ---------------------
    def _choose_color(self, color, label_field, var):
        color = colorchooser.askcolor(title="Choose color")[1]
        label_field.configure(bg=color)
        setattr(self, var, color)

    def _select_pdb1(self, field):
        self.pdb_1 = filedialog.askopenfilename(filetypes=[('pdb', '.pdb')], parent=self.compframe)
        self._configure_entry_field(field, self.pdb_1)

    def _choose_color1(self, color, label_field):
        color = colorchooser.askcolor(title="Choose color")[1]
        label_field.configure(bg=color)
        self.color1 = color

    def _select_pdb2(self, field):
        self.pdb_2 = filedialog.askopenfilename(filetypes=[('pdb', '.pdb')], parent=self.compframe)
        self._configure_entry_field(field, self.pdb_2)

    def _choose_color2(self, color, label_field):
        color = colorchooser.askcolor(title="Choose color")[1]
        label_field.configure(bg=color)
        self.color2 = color

    def _select_psf1_compare(self, field):
        self.psf_1 = filedialog.askopenfilename(title='Select protein structure file file', filetypes=[('psf', '.psf'), ('pdb', '.pdb'), ('gro', '.gro')], parent=self.DcdWaterWireFrame)
        self._configure_entry_field(field, self.psf_1)

    def _select_dcd1_compare(self, field):
        self.dcd_1 = filedialog.askopenfilenames(title='Select trajectory files', filetypes=[('dcd', '.dcd'), ('xtc', '.xtc'), ('trr', '.trr')],  parent=self.compframe)
        self._configure_entry_field(field, self.dcd_1)

    def _choose_dcd_color1(self, color, label_field):
        color = colorchooser.askcolor(title="Choose color")[1]
        label_field.configure(bg=color)
        self.color_dcd1 = color

    def _select_psf2_compare(self, field):
        self.psf_2 = filedialog.askopenfilename(title='Select protein structure file file', filetypes=[('psf', '.psf'), ('pdb', '.pdb'), ('gro', '.gro')], parent=self.DcdWaterWireFrame)
        self._configure_entry_field(field, self.psf_2)

    def _select_dcd2_compare(self, field):
        self.dcd_2 = filedialog.askopenfilenames(title='Select trajectory files', filetypes=[('dcd', '.dcd'), ('xtc', '.xtc'), ('trr', '.trr')],  parent=self.compframe)
        self._configure_entry_field(field, self.dcd_2)

    def _choose_dcd_color2(self, color, label_field):
        color = colorchooser.askcolor(title="Choose color")[1]
        label_field.configure(bg=color)
        self.color_dcd2 = color

    def _set_compare_result_folder(self, field):
        self.compare_results_folder = filedialog.askdirectory(parent=self.compframe)
        self._configure_entry_field(field, self.compare_results_folder)

    def _init_pdb_comparison(self, comp_type, pdb1, pdb2, color1='#1b3ede',color2='#21c25f'):
        if not hasattr(self, 'compare_results_folder') or self.compare_results_folder is None:
            print('WARNING: Please select the location of the workfolder!')
        else:
            comp = CompareTwo('pdb', pdb1=pdb1, pdb2=pdb2, target_folder=self.compare_results_folder, plot_parameters=self.plot_parameters)
            additional_donors = self._create_list_from_sting(self.dcd_comp_selected_donors_pdb.get())
            additional_acceptors = self._create_list_from_sting(self.dcd_comp_selected_acceptors_pdb.get())
            comp.calculate_graphs(graph_type=comp_type, max_water=self.max_water_comp.get(), include_backbone_sidechain=self.include_backbone_sidechain_comp.get(), include_waters=self.include_waters_comp.get(), distance=self.comp_distance.get(), cut_angle=self.comp_cut_angle.get(), check_angle=self.comp_use_angle.get(), selection=self.pdb_comp_selection_string.get(), additional_donors=additional_donors, additional_acceptors=additional_acceptors)
            comp.plot_graph_comparison(color1=color1, color2=color2, label_nodes=True, label_edges=True, calcualte_distance=self.calculate_distance_differences_comp.get())
            comp.plot_graph_comparison(color1=color1, color2=color2, label_nodes=False, label_edges=False)

            if self.color_propka_on_compare.get():
                comp.plot_graph_comparison(color1=color1, color2=color2, label_nodes=True, label_edges=True, color_propka=True, node_color_selection=self.selected_nodes_for_color_on_compare.get(), node_color_map=self.selected_color_map_on_compare.get(), calcualte_distance=self.calculate_distance_differences_comp.get())
                comp.plot_graph_comparison(color1=color1, color2=color2, label_nodes=False, label_edges=False, color_propka=True, node_color_selection=self.selected_nodes_for_color_on_compare.get(), node_color_map=self.selected_color_map_on_compare.get())
            if self.color_data_on_compare.get():
                comp.plot_graph_comparison(color1=color1, color2=color2, label_nodes=True, label_edges=True, color_data=True, node_color_selection=self.selected_nodes_for_color_on_compare.get(), node_color_map=self.selected_color_map_on_compare.get(), calcualte_distance=self.calculate_distance_differences_comp.get())
                comp.plot_graph_comparison(color1=color1, color2=color2, label_nodes=False, label_edges=False, color_data=True, node_color_selection=self.selected_nodes_for_color_on_compare.get(), node_color_map=self.selected_color_map_on_compare.get())

            comp.logger.info('Calculation completed')


    def _construct_compare_graphs(self, psf1, psf2, dcd1, dcd2, ):
        if not hasattr(self, 'compare_results_folder') or self.compare_results_folder is None:
            print('WARNING: Please select the location of the workfolder!')
        else:
            self.comp = CompareTwo('dcd', psf1=psf1, psf2=psf2, dcd1=dcd1, dcd2=dcd2, target_folder=self.compare_results_folder, name1=self.compare_dcd1_name.get(), name2=self.compare_dcd2_name.get(), plot_parameters=self.plot_parameters)
            additional_donors = self._create_list_from_sting(self.pdb_comp_selected_donors_pdb.get())
            additional_acceptors = self._create_list_from_sting(self.pdb_comp_selected_acceptors_pdb.get())
            self.comp.calculate_graphs(graph_type='water_wire', max_water=self.max_water_comp_dcd.get(), distance=self.comp_distance.get(), cut_angle=self.comp_cut_angle.get(), check_angle=self.comp_use_angle.get(), selection=self.dcd_comp_selection_string.get(), additional_donors=additional_donors, additional_acceptors=additional_acceptors)

                # -------------------water_wire_frame -----------------------
            if self.water_wire_comp_frame_dcd: self.water_wire_comp_frame_dcd.destroy()
            self.water_wire_comp_frame_dcd = ttk.LabelFrame(self.dcd_compare_tab, text='Water wire network')
            self.water_wire_comp_frame_dcd.grid(self._crate_frame_grid(self.compare_row))
            self.water_wire_comp_frame_dcd.columnconfigure(0, weight=1)
            self.water_wire_comp_frame_dcd.columnconfigure(1, weight=1)

            self.min_occupancy_comp = tk.DoubleVar(value=10)
            tk.Label(self.water_wire_comp_frame_dcd, text='Minimum H-bond occupancy (%)', anchor='w', bg='white', fg='black').grid(row=1, column=0, sticky='W')
            ttk.Spinbox(self.water_wire_comp_frame_dcd, textvariable=self.min_occupancy_comp, from_=1, to=100, validate="key", validatecommand=(self.ifnum_cmd, '%S', '%P', 0, 100)).grid(row=1, column=1, sticky="EW")
            tk.Button(self.water_wire_comp_frame_dcd, text='Compare water wire network', command=lambda:self._plot_dcd_comparison(color1=self.color_dcd1, color2=self.color_dcd2), width=self.button_width, bg='white', fg='black', highlightbackground='white').grid(self._create_big_button_grid(2), columnspan=2)


    def _plot_dcd_comparison(self, color1='#1b3ede',color2='#21c25f'):
        self.comp.construct_comparison_objects(occupancy=float(self.min_occupancy_comp.get())/100)
        self.comp.plot_graph_comparison(color1=color1, color2=color2, label_nodes=True, label_edges=True)
        self.comp.plot_graph_comparison(color1=color1, color2=color2, label_nodes=False, label_edges=False)
        self.comp.logger.info('Calculation completed')


#--------------------- COMMON ---------------------

    def _plot_conserved_graphs(self, c, plot_linear_length, plot_induvidual_graph, plot_difference_graph, plot_distance_graph, cth=0.9, occupancy=None, eps=1.5):
        c.get_conserved_graph(conservation_threshold=cth, occupancy=occupancy, eps=eps)
        c.plot_conserved_graph(label_nodes=True, label_edges=True)
        c.plot_conserved_graph(label_nodes=False, label_edges=False)
        if plot_linear_length:
            c.plot_linear_lenghts(occupancy=occupancy, label_nodes=True)
            c.plot_linear_lenghts(occupancy=occupancy, label_nodes=False)
        if plot_induvidual_graph:
            c.plot_graphs(label_nodes=True, label_edges=True, occupancy=occupancy)
            c.plot_graphs(label_nodes=False, label_edges=False, occupancy=occupancy)
        if plot_distance_graph:
            c.plot_graphs(label_nodes=True, label_edges=True, occupancy=occupancy, calcualte_distances=True)
        if self.color_propka.get():
            c.plot_graphs(label_nodes=True, label_edges=True, occupancy=occupancy, color_propka=True, node_color_selection=self.selected_nodes_for_color.get(), node_color_map=self.selected_color_map.get(), calcualte_distances=plot_distance_graph)
            c.plot_graphs(label_nodes=False, label_edges=False, occupancy=occupancy, color_propka=True, node_color_selection=self.selected_nodes_for_color.get(), node_color_map=self.selected_color_map.get())
        if self.color_data.get():
            c.plot_graphs(label_nodes=True, label_edges=True, occupancy=occupancy, color_data=True, node_color_selection=self.selected_nodes_for_color.get(), node_color_map=self.selected_color_map.get(), calcualte_distances=plot_distance_graph)
            c.plot_graphs(label_nodes=False, label_edges=False, occupancy=occupancy, color_data=True, node_color_selection=self.selected_nodes_for_color.get(), node_color_map=self.selected_color_map.get())

        if plot_difference_graph:
            c.plot_difference(label_nodes=True, label_edges=True)
            c.plot_difference(label_nodes=False, label_edges=False)
        c.logger.info('Calculation completed\n'+'-'*20)

    def custom_selection_string(self, parent_frame, row):
        # self._selection_string = tk.StringVar(value='protein')
        selection_entry = tk.Entry(parent_frame, state='disabled', bg='white', fg='black', highlightbackground='white', disabledbackground=self.gray, disabledforeground='black')
        self._configure_entry_field(selection_entry, value='protein')
        selection_entry.grid(row=row, column=1, sticky="EW")
        selected_donors = tk.StringVar()
        selected_acceptors = tk.StringVar()
        self.custom_selection_button = tk.Button(parent_frame, text='Selection', command=lambda:self.selection_string_popup(selection_entry, selected_donors, selected_acceptors), takefocus=False, bg='white', fg='black', highlightbackground='white')
        self.custom_selection_button.grid(row=row, column=0, sticky="W")
        return selection_entry, selected_donors, selected_acceptors

    def _save_plot_settings(self):
        formats = ['png']
        if self.eps_format.get(): formats.append('eps')
        if self.svg.get(): formats.append('svg')

        self.plot_parameters = {
                'edge_width': float(self.edge_width.get()),
                'node_label_size': float(self.node_label_size.get()),
                'edge_label_size': float(self.edge_label_size.get()),
                'node_size': float(self.node_size.get()),
                'graph_color': self.graph_color,
                'difference_graph_color':self.difference_graph_color,
                'water_node_color': self.water_node_color,
                'non_prot_color': self.non_prot_color,
                'plot_title_fontsize':float(self.plot_title_fontsize.get()),
                'plot_label_fontsize':float(self.plot_label_fontsize.get() ),
                'plot_tick_fontsize':float(self.plot_tick_fontsize.get()),
                'figsize': (int(self.plot_width.get()), int(self.plot_height.get())),
                'plot_resolution':float(self.plot_resolution.get()),
                'formats': formats,
                'show_chain_label': self.show_chain_label.get()
            }

    #SOURCE: https://stackoverflow.com/questions/10020885/creating-a-popup-message-box-with-an-entry-field
    def selection_string_popup(self, selection_entry, selected_donors, selected_acceptors):
        self.popup=popupWindow(self.master)
        self.popup.custom_selection_sting(selection_entry, selected_donors, selected_acceptors)
        self.custom_selection_button['state'] = 'disabled'
        self.master.wait_window(self.popup.top)
        self.custom_selection_button['state'] = 'normal'
        self._configure_entry_field(selection_entry, value=self.popup._sel_string)
        selected_donors.set(self.popup._sel_donors)
        selected_acceptors.set(self.popup._sel_acceptors)

    def node_color_selelection_pop_up(self, selected_nodes_for_color, selected_color):
        self.popup=popupWindow(self.master)
        self.popup.node_color_selection(selected_nodes_for_color, selected_color)
        self.custom_selection_button['state'] = 'disabled'
        self.master.wait_window(self.popup.top)
        self.custom_selection_button['state'] = 'normal'
        selected_nodes_for_color.set(self.popup._nc_sel_string)
        selected_color.set(self.popup._nc_sel_color)


    def _configure_entry_field(self, field, value=None):
        field.configure(state='normal', bg='white', fg='black')
        field.delete(0, 'end')
        if value: field.insert(0, str(value))
        field.configure(state='disabled')

    def _add_horisontal_scroll(self, target, row=1, column=0):
        scroll = tk.Scrollbar(target, orient='horizontal')
        scroll.grid(row=row, column=column, sticky='EW')
        return scroll

    def _create_list_from_sting(self, atom_list):
        return list(filter(None, re.split('\s|;|,|\*|\n|\.', atom_list)))

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
            'pady': (5, 5),
            'sticky': "EW"
        }

    def _create_frame(self):
        style = ttk.Style()
        style.theme_create( 'style', settings={
            '.': {'configure': {'background': 'white', 'relief': 'flat', 'takefocus':'false'}},
            'TNotebook': {'configure': {'tabmargins': [2, 5, 0, 0], 'background':'white'} },
            'TFrame': {'configure': {'relief': 'flat', 'padding': [6,5,6,15], 'background':'white', 'highlightbackground':'white'}},
            'TLabelframe': {'configure': {'relief': 'flat', 'padding': [6,5,6,15], 'background':'white', 'highlightbackground':'white'}},
            'TLabelframe.Label': {'configure': {'font': ('Helvetica', 13, 'bold')}, 'background':'white'},
            'TSpinbox': {'configure': {'background':'white', 'foreground':'black', 'insertcolor': 'black'}},
            'TCombobox': {'configure': {'background':'white', 'foreground':'black'}},
            'TNotebook.Tab': {
                    'configure': {'padding': [8, 4], 'background': self.gray },
                    'map': {'background': [('selected', 'white')]}}
        } )
        style.theme_use('style')

        # https://blog.teclado.com/tkinter-scrollable-frames/
        canvas = tk.Canvas(self.master, bg='white', highlightthickness=0)

        tab_parnt = ttk.Notebook(canvas)
        self.mainframe = ttk.Frame(tab_parnt)
        self.dcdframe = ttk.Frame(tab_parnt)
        self.compframe = ttk.Frame(tab_parnt)
        self.plotsettings = ttk.Frame(tab_parnt)

        tab_parnt.add(self.mainframe, text='Crystal structure analysis')
        tab_parnt.add(self.dcdframe, text='MD trajectory analysis')
        tab_parnt.add(self.compframe, text='Compare 2 structures')
        tab_parnt.add(self.plotsettings, text='Plot settings')

        scrollbar = tk.Scrollbar(self.master, orient="vertical", command=canvas.yview, bg='white')

        tab_parnt.bind(
            "<Configure>",
            lambda e: canvas.configure(
                scrollregion=canvas.bbox("all")
            )
        )

        canvas.create_window((50, 0), window=tab_parnt, anchor='nw')

        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")



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
    root.configure(bg='white')
    view = View(root)
    view.main_modal()
    root.mainloop()
