import tkinter as tk
from tkinter import ttk

def csa_view(self):

    # ----------------------- inputFrame -----------------------

    self.inputFrame = ttk.LabelFrame(self.mainframe, text='Input locations')
    self.inputFrame.grid(self._crate_frame_grid(0))
    self.inputFrame.columnconfigure(0, weight=1)
    self.inputFrame.columnconfigure(1, weight=1)

    s1 = self._add_horisontal_scroll(self.inputFrame, row=2, column=1)
    self._input_folder = tk.Entry(self.inputFrame, state='disabled', xscrollcommand=s1.set, bg='white', fg='black', highlightbackground='white', disabledbackground=self.gray, disabledforeground='black')
    tk.Button(self.inputFrame, text='Select PDB folder', command=lambda:self._select_pdb_root_folder(self._input_folder), takefocus=False, bg='white', fg='black', highlightbackground='white').grid(row=1, column=0, sticky="EW")
    self._input_folder.grid(row=1, column=1, sticky="EW", columnspan=2)
    s1.configure(command=self._input_folder.xview, bg='white')

    tk.Button(self.inputFrame, text='Select reference file', command=lambda:self._select_reference_pdb(self._input_pdb),  highlightbackground='white', bg='white', fg='black').grid(row=4, column=0, sticky="EW")
    s2 = self._add_horisontal_scroll(self.inputFrame, row=5, column=1)
    self._input_pdb = tk.Entry(self.inputFrame, state='disabled', xscrollcommand=s2.set, bg='white', fg='black',  highlightbackground='white',disabledbackground=self.gray, disabledforeground='black')
    self._input_pdb.grid(row=4, column=1, sticky="EW")
    s2.configure(command=self._input_pdb.xview, bg='white')


    self.sequance_identity_threshold = tk.DoubleVar(value=75)
    tk.Label(self.inputFrame, text='Minimum sequence identity (%)', bg='white', fg='black').grid(row=6, column=0, sticky='W')
    ttk.Spinbox(self.inputFrame, textvariable=self.sequance_identity_threshold, from_=1, to=100, validate="key", validatecommand=(self.ifnum_cmd, '%S', '%P', 0, 100)).grid(row=6, column=1, sticky="EW")

    # ----------------------- waterClusterFrame -----------------------

    self.waterClusterFrame = ttk.LabelFrame(self.mainframe, text='Water cluster analysis with DBSCAN')
    self.waterClusterFrame.grid(self._crate_frame_grid(7))
    self.waterClusterFrame.columnconfigure(0, weight=1)

    # tk.Button(self.waterClusterFrame, text='Parameter analysis', command=self._perform_parameter_analysis, width=self.button_width).grid(row=0, column=0, padx=(self.padx,self.padx), pady=(self.pady,self.pady), sticky="EW")
    waterClusterParameterFrame = tk.Frame(self.waterClusterFrame, bg='white')
    waterClusterParameterFrame.grid(row=1, column=0, columnspan=4, sticky="EW")
    lwcp = tk.Label(waterClusterParameterFrame, text='For the DBSCAN analysis the default parameters are recommended to use.', anchor='w', bg='white', fg='black')
    lwcp.grid(row=0, columnspan=4, sticky='EW')
    lwcp.config(font=("Helvetica", 11), bg='white', fg='black')

    tk.Label(waterClusterParameterFrame, text='DBSCAN eps:', bg='white', fg='black').grid(row=1, column=0, sticky="W")
    self.eps = tk.DoubleVar(value=1.4)
    ttk.Spinbox(waterClusterParameterFrame, textvariable=self.eps, width=5, validate="key", validatecommand=(self.ifnum_cmd, '%S', '%P', 0, 100)).grid(row=1, column=1, sticky="W")

    tk.Label(waterClusterParameterFrame, text='Superimposition RMS threshold:', bg='white', fg='black').grid(row=1, column=2, sticky="E")
    self.superimposition_threshold = tk.DoubleVar(value=5)
    ttk.Spinbox(waterClusterParameterFrame, textvariable=self.superimposition_threshold, width=5, validate="key", validatecommand=(self.ifnum_cmd, '%S', '%P', 0, 100)).grid(row=1, column=3, sticky="E")

    tk.Button(self.waterClusterFrame, text='Calculate water clusters', command=self._init_water_clusters, width=self.button_width, bg='white', fg='black',  highlightbackground='white').grid(self._create_big_button_grid(2))


    # ----------------------- conservedNetworkFrame -----------------------

    self.conservedNetworkFrame = ttk.LabelFrame(self.mainframe, text='Conserved network analysis with Bridge')
    self.conservedNetworkFrame.grid(self._crate_frame_grid(8))
    self.conservedNetworkFrame.columnconfigure(0, weight=1)
    self.conservedNetworkFrame.columnconfigure(1, weight=1)

    selsting_frame = tk.Frame(self.conservedNetworkFrame, bg='white')
    selsting_frame.grid(row=8, column=0, columnspan=2, sticky="EW")
    selsting_frame.columnconfigure(1, weight=1)
    self.selection_string, self.selected_donors_pdb, self.selected_acceptors_pdb = self.custom_selection_strin(selsting_frame, 8)


    tk.Label(self.conservedNetworkFrame, text='H-bond criteria ', anchor="w", bg='white', fg='black').grid(row=9, column=0, sticky='W')
    hcritera_frame = tk.Frame(self.conservedNetworkFrame, bg='white')
    hcritera_frame.grid(row=9, column=1, columnspan=5, sticky="EW")
    self.c_distance = tk.DoubleVar(value=3.5)
    self.c_cut_angle = tk.DoubleVar(value=60)
    ttk.Spinbox(hcritera_frame, textvariable=self.c_distance, from_=0, to=5, width=5, validate="key", validatecommand=(self.ifnum_cmd, '%S', '%P', 0, 5)).grid(row=9, column=1, sticky='W')
    tk.Label(hcritera_frame, text='Å distance     ', anchor="w", bg='white', fg='black').grid(row=9, column=2, sticky='W')
    self.c_use_angle = tk.BooleanVar()
    self.c_use_angle.set(False)
    tk.Checkbutton(hcritera_frame, variable=self.c_use_angle, anchor="w", bg='white', fg='black').grid(row=9, column=3, sticky='E')
    ttk.Spinbox(hcritera_frame, textvariable=self.c_cut_angle, from_=0, to=180, width=5, validate="key", validatecommand=(self.ifnum_cmd, '%S', '%P', 0, 180)).grid(row=9, column=4, sticky='W')
    tk.Label(hcritera_frame, text='degrees angle', anchor="w", bg='white', fg='black').grid(row=9, column=5, sticky='W')

    self.conservation_threshold = tk.DoubleVar(value=90)
    tk.Label(self.conservedNetworkFrame, text='Conservation of H-bonding groups across structures', anchor="w", bg='white', fg='black').grid(row=10, column=0, sticky='W')
    ttk.Spinbox(self.conservedNetworkFrame, textvariable=self.conservation_threshold, from_=1, to=100, width=5, validate="key", validatecommand=(self.ifnum_cmd, '%S', '%P', 0, 100)).grid(row=10, column=1, sticky="W")

    tk.Label(self.conservedNetworkFrame, text='Color nodes by:', anchor="w", bg='white', fg='black').grid(row=11, column=0, sticky='W')
    color_plots_crystal = tk.Frame(self.conservedNetworkFrame, bg='white')
    color_plots_crystal.grid(row=11, column=1, columnspan=3, sticky="EW", pady=(0,10))

    self.color_propka = tk.BooleanVar()
    self.color_propka.set(False)
    tk.Checkbutton(color_plots_crystal, text='Propka file   ', variable=self.color_propka, anchor="w", bg='white', fg='black').grid(row=11, column=1, sticky='E')

    self.color_data = tk.BooleanVar()
    self.color_data.set(False)
    tk.Checkbutton(color_plots_crystal, text='Custom _data.txt    ', variable=self.color_data, anchor="w", bg='white', fg='black').grid(row=11, column=2, sticky='E')
    # tk.Label(color_plots_crystal, text='Restrict amino acid type', anchor="w", bg='white', fg='black').grid(row=11, column=3, sticky='E')

    tk.Label(self.conservedNetworkFrame, text='Plot for each structure:', anchor="w", bg='white', fg='black').grid(row=12, column=0, sticky='W')
    each_plots_crystal = tk.Frame(self.conservedNetworkFrame)
    each_plots_crystal.grid(row=12, column=1, columnspan=3, sticky="EW", pady=(0,10))
    self.is_induvidual_graph = tk.BooleanVar()
    tk.Checkbutton(each_plots_crystal, text='Individual network    ', variable=self.is_induvidual_graph, anchor="w", bg='white', fg='black').grid(row=12, column=1, sticky='E')

    self.is_difference_graph = tk.BooleanVar()
    tk.Checkbutton(each_plots_crystal, text='Difference graph    ', variable=self.is_difference_graph, anchor="w", bg='white', fg='black').grid(row=12, column=2, sticky='E')

    self.is_linear_lenght_plot = tk.BooleanVar()
    tk.Checkbutton(each_plots_crystal, text='Linear lengths', variable=self.is_linear_lenght_plot, anchor="w", bg='white', fg='black').grid(row=12, column=3, sticky='E')



    # ----------------------- HbondNetworkFrame -----------------------

    self.HbondNetworkFrame = ttk.LabelFrame(self.conservedNetworkFrame, text='H-bond network')
    self.HbondNetworkFrame.grid(self._crate_frame_grid(14))
    self.HbondNetworkFrame.columnconfigure(0, weight=1)
    # self.HbondNetworkFrame.columnconfigure(1, weight=1)

    self.include_waters_hbond = tk.BooleanVar()
    self.include_waters_hbond.set(True)
    tk.Checkbutton(self.HbondNetworkFrame, text='Include crystallographic waters', variable=self.include_waters_hbond, anchor="w", bg='white', fg='black').grid(self._create_big_button_grid(14))

    self.useWaterCoords = tk.BooleanVar()
    tk.Checkbutton(self.HbondNetworkFrame, text='Use water cluster coordinates', variable=self.useWaterCoords, anchor="w", bg='white', fg='black').grid(self._create_big_button_grid(15))

    self.include_backbone_sidechain = tk.BooleanVar()
    tk.Checkbutton(self.HbondNetworkFrame, text='Include sidechain-backbone interactions', variable=self.include_backbone_sidechain, anchor="w", bg='white', fg='black').grid(self._create_big_button_grid(16))

    tk.Button(self.HbondNetworkFrame, text='Calculate conserved H-bond network', command=lambda:self._init_pdb_conserved_graph_analysis('hbond'), width=self.button_width, bg='white', fg='black', highlightbackground='white').grid(self._create_big_button_grid(17))


    # ----------------------- WaterWireFrame -----------------------

    self.WaterWireFrame = ttk.LabelFrame(self.conservedNetworkFrame, text='Water wire network')
    self.WaterWireFrame.grid(self._crate_frame_grid(16))
    self.WaterWireFrame.columnconfigure(0, weight=1)
    self.WaterWireFrame.columnconfigure(1, weight=1)

    self.max_water = tk.IntVar(value=3)
    tk.Label(self.WaterWireFrame, text=' Maximum number of waters in the bridge', anchor="w", bg='white', fg='black').grid(row=17, column=0, sticky='W')
    ttk.Combobox(self.WaterWireFrame, textvariable=self.max_water, values=[1,2,3,4,5], state='readonly').grid(row=17, column=1, sticky="EW")
    tk.Button(self.WaterWireFrame, text='Calculate conserved water wire network', command=lambda:self._init_pdb_conserved_graph_analysis('water_wire'), width=self.button_width, bg='white', fg='black',  highlightbackground='white').grid(self._create_big_button_grid(18), columnspan=2)

    # self.completedText = tk.StringVar(value='')
    # self.completed = tk.Label(self.conservedNetworkFrame, textvariable=self.completedText)
    # self.completed.grid(row=17, column=0)

    return
