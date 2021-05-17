import tkinter as tk
from tkinter import ttk

def ta_view(self):

   #--------target folder select---------------------
    inital_sim_settings = ttk.LabelFrame(self.dcdframe, text='Parameter set up for the analysis')
    inital_sim_settings.grid(self._crate_frame_grid(0))
    inital_sim_settings.columnconfigure(1, weight=1)

    self._target_folder = None
    s1 = self._add_horisontal_scroll(inital_sim_settings, row=2, column=1)
    self._input_target = tk.Entry(inital_sim_settings, state='disabled', xscrollcommand=s1.set)
    tk.Button(inital_sim_settings, text='Location of workfolder', command=lambda:self._select_dcd_workfolder(self._input_target), takefocus=False).grid(row=1, column=0, sticky="EW",)
    self._input_target.grid(row=1, column=1, sticky="EW")
    s1.configure(command=self._input_target.xview)

    self.sim_max_water = tk.IntVar(value=3)
    tk.Label(inital_sim_settings, text='Maximum number of water molecules allowed in the bridge', anchor='w').grid(row=3, column=0, sticky='W')
    ttk.Combobox(inital_sim_settings, textvariable=self.sim_max_water, values=[1,2,3,4,5], state='readonly').grid(row=3, column=1, sticky="EW")


    self.sim_selection_string = tk.StringVar(value='protein')
    # tk.Label(selsting_frame, text='  Selection string', anchor="w").grid(row=8, column=0, sticky='W')
    # tk.Entry(selsting_frame, textvariable=self.selection_string).grid(row=8, column=1, sticky="EW")

    sim_critera_frame = tk.Frame(inital_sim_settings)
    sim_critera_frame.grid(row=4, column=0, columnspan=5, sticky="EW")
    self.sim_distance = tk.DoubleVar(value=3.5)
    self.sim_cut_angle = tk.DoubleVar(value=60)
    tk.Label(sim_critera_frame, text='H-bond cut-off criteria ', anchor="w").grid(row=4, column=0, sticky='E')
    ttk.Spinbox(sim_critera_frame, textvariable=self.sim_distance, from_=0, to=5, width=13, validate="key", validatecommand=(self.ifnum_cmd, '%S', '%P', 0, 5)).grid(row=4, column=1, sticky='W')
    tk.Label(sim_critera_frame, text='Å distance and ', anchor="w").grid(row=4, column=2, sticky='W')
    ttk.Spinbox(sim_critera_frame, textvariable=self.sim_cut_angle, from_=0, to=180, width=13, validate="key", validatecommand=(self.ifnum_cmd, '%S', '%P', 0, 180)).grid(row=4, column=3, sticky='W')
    tk.Label(sim_critera_frame, text='degrees angle', anchor="w").grid(row=4, column=4, sticky='W')

    #--------------------------- dcd select------------

    self.selectSimFrame = ttk.LabelFrame(self.dcdframe, text='Select simulation')
    self.selectSimFrame.grid(self._crate_frame_grid(4))
    self.selectSimFrame.columnconfigure(0, weight=0)
    self.selectSimFrame.columnconfigure(1, weight=1)
    tk.Label(self.selectSimFrame, text='Select and compute H-bond graph for one simulation at a time. After the calculation is completed \nyou can construct the water wire network or compute graphs from other simulations.', anchor='w', justify='left').grid(row=4, columnspan=2, sticky='W')

    tk.Button(self.selectSimFrame, text='Select PSF', command=self._select_psf_file).grid(row=5, column=0, sticky="EW")
    s2 = self._add_horisontal_scroll(self.selectSimFrame, row=6, column=1)
    self._input_psf = tk.Entry(self.selectSimFrame, state='disabled', xscrollcommand=s2.set)
    self._input_psf.grid(row=5, column=1, sticky="EW")
    s2.configure(command=self._input_psf.xview)

    tk.Button(self.selectSimFrame, text='Select DCDs', command=self._select_dcd_files).grid(row=7, column=0, sticky="EW")
    s3 = self._add_horisontal_scroll(self.selectSimFrame, row=8, column=1)
    self._input_dcd = tk.Entry(self.selectSimFrame, state='disabled', xscrollcommand=s3.set)
    self._input_dcd.grid(row=7, column=1, sticky="EW")
    s3.configure(command=self._input_dcd.xview)

    tk.Label(self.selectSimFrame, text='Name as: ', anchor='w').grid(row=9, column=0, sticky='W')
    # self.sim_name1 = tk.StringVar(value='sim1')
    self.sim_name = tk.Entry(self.selectSimFrame)
    self.sim_name.insert(0, 'sim1') # remove when test resolved
    self.sim_name.grid(row=9, column=1, sticky="EW")
    l1 = tk.Label(self.selectSimFrame, text='Give a unique name to your simulation. Calculations with the same name are overwritten.', anchor='w')
    l1.grid(row=10, column=1, sticky='EW')
    l1.config(font=("Helvetica", 11))

    tk.Button(self.selectSimFrame, text='Construct graph', command=self._construct_sim_graphs).grid(self._create_big_button_grid(11, column=1))

    # ----------------------- DcdWaterWireFrame -----------------------
    self.DcdWaterWireFrame = ttk.LabelFrame(self.dcdframe, text='Water wire network')
    self.DcdWaterWireFrame.grid(self._crate_frame_grid(15))
    self.DcdWaterWireFrame.columnconfigure(0, weight=1)
    # self.DcdWaterWireFrame.columnconfigure(1, weight=1)
    # self.DcdWaterWireFrame.columnconfigure(2, weight=1)

    self.row=16
    tk.Button(self.DcdWaterWireFrame, text='Select graphs to compare', command=lambda:self._load_graph_files(self.row)).grid(self._create_big_button_grid(15))

    return
