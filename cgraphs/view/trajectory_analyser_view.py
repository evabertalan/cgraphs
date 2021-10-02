import tkinter as tk
from tkinter import ttk

def ta_view(self):

   #--------target folder select---------------------
    inital_sim_settings = ttk.LabelFrame(self.dcdframe, text='Parameter set up for the analysis')
    inital_sim_settings.grid(self._crate_frame_grid(0))
    inital_sim_settings.columnconfigure(0, weight=1)
    inital_sim_settings.columnconfigure(1, weight=1)

    self._target_folder = None
    s1 = self._add_horisontal_scroll(inital_sim_settings, row=1, column=1)
    self._input_target = tk.Entry(inital_sim_settings, state='disabled', xscrollcommand=s1.set, bg='white', fg='black',  highlightbackground='white', disabledbackground=self.gray, disabledforeground='black')
    tk.Button(inital_sim_settings, text='Location of workfolder', bg='white', fg='black', command=lambda:self._select_dcd_workfolder(self._input_target), takefocus=False,  highlightbackground='white').grid(row=0, column=0, sticky="EW",)
    self._input_target.grid(row=0, column=1, sticky="EW")
    s1.configure(command=self._input_target.xview, bg='white')

    self.sim_max_water = tk.IntVar(value=3)
    tk.Label(inital_sim_settings, text='Maximum number of waters in the bridge', anchor='w',  bg='white', fg='black').grid(row=2, column=0, sticky='W')
    ttk.Combobox(inital_sim_settings, textvariable=self.sim_max_water, values=[1,2,3,4,5], state='readonly').grid(row=2, column=1, sticky="EW")

    tk.Label(inital_sim_settings, text='H-bond criteria ', anchor="w",  bg='white', fg='black').grid(row=3, column=0, sticky='W')
    sim_critera_frame = tk.Frame(inital_sim_settings,  bg='white')
    sim_critera_frame.grid(row=3, column=1, columnspan=4, sticky="EW")
    self.sim_distance = tk.DoubleVar(value=3.5)
    self.sim_cut_angle = tk.DoubleVar(value=60)
    ttk.Spinbox(sim_critera_frame, textvariable=self.sim_distance, from_=0, to=5, width=5, validate="key", validatecommand=(self.ifnum_cmd, '%S', '%P', 0, 5)).grid(row=3, column=1, sticky='W')
    tk.Label(sim_critera_frame, text='Å distance     ', anchor="w",  bg='white', fg='black').grid(row=3, column=2, sticky='W')
    self.sim_use_angle = tk.BooleanVar()
    self.sim_use_angle.set(True)
    tk.Checkbutton(sim_critera_frame, variable=self.sim_use_angle, anchor="w",  bg='white', fg='black').grid(row=3, column=3, sticky='E')
    ttk.Spinbox(sim_critera_frame, textvariable=self.sim_cut_angle, from_=0, to=180, width=5, validate="key", validatecommand=(self.ifnum_cmd, '%S', '%P', 0, 180)).grid(row=3, column=4, sticky='W')
    tk.Label(sim_critera_frame, text='degrees angle', anchor="w",  bg='white', fg='black').grid(row=3, column=5, sticky='W')

    selsting_frame = tk.Frame(inital_sim_settings, bg='white')
    selsting_frame.grid(row=4, column=0, columnspan=2, sticky="EW")
    selsting_frame.columnconfigure(1, weight=1)
    self.sim_selection_string, self.sim_selected_donors_pdb, self.sim_selected_acceptors_pdb = self.custom_selection_strin(selsting_frame, 2)
    #--------------------------- dcd select------------

    self.selectSimFrame = ttk.LabelFrame(self.dcdframe, text='Select simulation')
    self.selectSimFrame.grid(self._crate_frame_grid(4))
    self.selectSimFrame.columnconfigure(0, weight=0)
    self.selectSimFrame.columnconfigure(1, weight=1)
    sl = tk.Label(self.selectSimFrame, text='Compute H-bond graph for one simulation at a time. Once the graph is computed, the water wire for the current simulation may be calculated.\nAlternatively, the graph for another simulation may be computed.', anchor='w', justify='left',  bg='white', fg='black')
    sl.grid(row=4, columnspan=2, sticky='W')
    sl.config(font=("Helvetica", 11))

    tk.Button(self.selectSimFrame, text='Select PSF', command=self._select_psf_file, bg='white', fg='black', highlightbackground='white').grid(row=5, column=0, sticky="EW")
    s2 = self._add_horisontal_scroll(self.selectSimFrame, row=6, column=1)
    self._input_psf = tk.Entry(self.selectSimFrame, state='disabled', xscrollcommand=s2.set, bg='white', fg='black', highlightbackground='white', disabledbackground=self.gray, disabledforeground='black')
    self._input_psf.grid(row=5, column=1, sticky="EW")
    s2.configure(command=self._input_psf.xview,  bg='white')

    tk.Button(self.selectSimFrame, text='Select DCDs', command=self._select_dcd_files, bg='white', fg='black', highlightbackground='white').grid(row=7, column=0, sticky="EW")
    s3 = self._add_horisontal_scroll(self.selectSimFrame, row=8, column=1)
    self._input_dcd = tk.Entry(self.selectSimFrame, state='disabled', xscrollcommand=s3.set, bg='white', fg='black', highlightbackground='white', disabledbackground=self.gray, disabledforeground='black')
    self._input_dcd.grid(row=7, column=1, sticky="EW")
    s3.configure(command=self._input_dcd.xview, bg='white')

    tk.Label(self.selectSimFrame, text='Name as: ', anchor='w', bg='white', fg='black').grid(row=9, column=0, sticky='W')
    # self.sim_name1 = tk.StringVar(value='sim1')
    self.sim_name = tk.Entry(self.selectSimFrame, bg='white', fg='black', highlightbackground='white', insertbackground='black')
    self.sim_name.insert(0, 'sim1') # remove when test resolved
    self.sim_name.grid(row=9, column=1, sticky="EW")
    l1 = tk.Label(self.selectSimFrame, text='Each simulation must be given a unique name, else C-Graphs results will be overwritten.', anchor='w',bg='white', fg='black')
    l1.grid(row=10, column=1, sticky='EW')
    l1.config(font=("Helvetica", 11))

    tk.Button(self.selectSimFrame, text='Construct graph', command=self._construct_sim_graphs, bg='white', fg='black', highlightbackground='white').grid(self._create_big_button_grid(11, column=1))

    # ----------------------- DcdWaterWireFrame -----------------------
    self.DcdWaterWireFrame = ttk.LabelFrame(self.dcdframe, text='Water wire network')
    self.DcdWaterWireFrame.grid(self._crate_frame_grid(15))
    self.DcdWaterWireFrame.columnconfigure(0, weight=1)
    # self.DcdWaterWireFrame.columnconfigure(1, weight=1)
    # self.DcdWaterWireFrame.columnconfigure(2, weight=1)

    self.row=16
    tk.Button(self.DcdWaterWireFrame, text='Select graphs to compare', command=lambda:self._load_graph_files(self.row),  bg='white', fg='black', highlightbackground='white').grid(self._create_big_button_grid(15))

    return
