import tkinter as tk
from tkinter import ttk

def ta_view(self):


   #--------target folder select---------------------
    self.inital_sim_settings = tk.Frame(self.dcdframe)
    self.inital_sim_settings.grid(row=0, columnspan=3, sticky='EW', padx=(self.padx,self.padx), pady=(self.pady,self.pady), ipadx=self.ipadx, ipady=self.ipady)
    # self.inital_sim_settings.columnconfigure(0, weight=1)
    # self.inital_sim_settings.columnconfigure(1, weight=1)
    # self.inital_sim_settings.columnconfigure(2, weight=1)

    tk.Button(self.inital_sim_settings, text='Save results to:', command=self._select_target_folder).grid(row=1, column=0, sticky="EW",)
    s5 = self._add_horisontal_scroll(self.inital_sim_settings, row=2, column=1)
    self._input_target = tk.Entry(self.inital_sim_settings, state='disabled', xscrollcommand=s5.set)
    self._input_target.grid(row=1, column=1, sticky="EW")
    s5.configure(command=self._input_target.xview)

    self.sim_max_water = tk.StringVar(value='3')
    tk.Label(self.inital_sim_settings, text='Maximum number of water molecules allowed in the bridge', anchor="w").grid(row=3, column=0)
    ttk.Combobox(self.inital_sim_settings, textvariable=self.sim_max_water, values=['1','2','3','4','5']).grid(row=3, column=1, sticky="EW")

    #--------------------------- dcd select------------

    self.selectSimFrame = tk.LabelFrame(self.dcdframe, text='Select simulation')
    self.selectSimFrame.grid(row=4, columnspan=3, sticky='EW', padx=(self.padx,self.padx), pady=(self.pady,self.pady), ipadx=self.ipadx, ipady=self.ipady)
    self.selectSimFrame.columnconfigure(0, weight=0)
    self.selectSimFrame.columnconfigure(1, weight=1)
    # self.selectSimFrame.columnconfigure(2, weight=1)

    # tk.Button(self.selectSimFrame, text='Select protein PDB', command=self._select_sim_pdb_file).grid(row=5, column=0, sticky="EW")
    # s6 = self._add_horisontal_scroll(self.selectSimFrame, row=6, column=1)
    # self._input_sim_pdb = tk.Entry(self.selectSimFrame, state='disabled', xscrollcommand=s6.set)
    # self._input_sim_pdb.grid(row=5, column=1, sticky="EW")
    # s6.configure(command=self._input_sim_pdb.xview)

    tk.Button(self.selectSimFrame, text='Select PSF', command=self._select_psf_file).grid(row=5, column=0, sticky="EW")
    s3 = self._add_horisontal_scroll(self.selectSimFrame, row=6, column=1)
    self._input_psf = tk.Entry(self.selectSimFrame, state='disabled', xscrollcommand=s3.set)
    self._input_psf.grid(row=5, column=1, sticky="EW")
    s3.configure(command=self._input_psf.xview)

    tk.Button(self.selectSimFrame, text='Select DCDs', command=self._select_dcd_files).grid(row=7, column=0, sticky="EW")
    s4 = self._add_horisontal_scroll(self.selectSimFrame, row=8, column=1)
    self._input_dcd = tk.Entry(self.selectSimFrame, state='disabled', xscrollcommand=s4.set)
    self._input_dcd.grid(row=7, column=1, sticky="EW")
    s4.configure(command=self._input_dcd.xview)

    tk.Label(self.selectSimFrame, text='name as: ').grid(row=9, column=0, sticky="EW")
    # self.sim_name1 = tk.StringVar(value='sim1')
    self.sim_name = tk.Entry(self.selectSimFrame)
    self.sim_name.insert(0, 'sim1') # remove when test resolved
    self.sim_name.grid(row=9, column=1, sticky="EW")

    tk.Button(self.selectSimFrame, text='Construct graph', command=self._construct_sim_graphs).grid(row=10, column=1, padx=(self.padx,self.padx), pady=(self.pady,self.pady), sticky="EW")

    # self.dcdComputeInfo = tk.StringVar()
    self.dcd_compute = tk.Label(self.selectSimFrame)
    self.dcd_compute.grid(row=10, column=0)
    self.dcd_compute2 = tk.Label(self.selectSimFrame, text='Now you can calculate the water wire network or costruct \ngraphs from other simuliation and then calculate the conserved network.')
    self.dcd_compute2.grid(row=11, column=0)
    self.dcd_compute2.grid_forget()

    # ----------------------- DcdWaterWireFrame -----------------------
    self.DcdWaterWireFrame = tk.LabelFrame(self.dcdframe, text='Water wire network')
    self.DcdWaterWireFrame.grid(row=12, columnspan=3, sticky='EW', padx=(self.padx,self.padx), pady=(self.pady,self.pady), ipadx=self.ipadx, ipady=self.ipady)
    # self.DcdWaterWireFrame.columnconfigure(0, weight=0)
    # self.DcdWaterWireFrame.columnconfigure(1, weight=1)
    # self.DcdWaterWireFrame.columnconfigure(2, weight=1)

    # tk.Button(self.DcdWaterWireFrame, text='Select reference pdb', command=self._select_dcd_reference_file).grid(row=12, column=0, sticky="EW")
    # s6 = self._add_horisontal_scroll(self.DcdWaterWireFrame, row=13, column=1)
    # self._input_pdb_dcd = tk.Entry(self.DcdWaterWireFrame, state='disabled', xscrollcommand=s6.set)
    # self._input_pdb_dcd.grid(row=12, column=1, sticky="EW")
    # s6.configure(command=self._input_pdb_dcd.xview)

    self.sequance_identity_threshold_dcd = tk.StringVar(value='75')
    tk.Label(self.DcdWaterWireFrame, text='Minimum sequence identity (%)').grid(row=14, column=0)
    ttk.Spinbox(self.DcdWaterWireFrame, textvariable=self.sequance_identity_threshold_dcd, from_=1, to=100).grid(row=14, column=1, sticky="EW")

    self.conservation_threshold_dcd = tk.StringVar(value='90')
    tk.Label(self.DcdWaterWireFrame, text='Conservation of H-bonding groups across structures', anchor="w").grid(row=15, column=0)
    ttk.Spinbox(self.DcdWaterWireFrame, textvariable=self.conservation_threshold_dcd, from_=1, to=100).grid(row=15, column=1, sticky="EW")


    self.min_occupancy = tk.StringVar(value='10')
    tk.Label(self.DcdWaterWireFrame, text='Minimum H-bond occupancy (%').grid(row=16, column=0)
    ttk.Spinbox(self.DcdWaterWireFrame, textvariable=self.min_occupancy, from_=1, to=100).grid(row=16, column=1, sticky="EW")
    row=19
    tk.Button(self.DcdWaterWireFrame, text='Load graphs to compare:', command=lambda:self._load_graph_files(row)).grid(row=17, column=0, sticky="EW")
    self.LoadGraphFrame = tk.Frame(self.DcdWaterWireFrame)
    self.LoadGraphFrame.grid(row=18, columnspan=3, sticky='EW', padx=(self.padx,self.padx), pady=(self.pady,self.pady), ipadx=self.ipadx, ipady=self.ipady)





