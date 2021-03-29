import tkinter as tk
from tkinter import ttk

def ta_view(self):


   #--------target folder select---------------------
    self.inital_sim_settings = ttk.LabelFrame(self.dcdframe, text='Parameter set up for the analysis')
    self.inital_sim_settings.grid(self._crate_frame_grid(0))
    # self.inital_sim_settings.columnconfigure(0, weight=1)
    self.inital_sim_settings.columnconfigure(1, weight=1)
    # self.inital_sim_settings.columnconfigure(2, weight=1)

    tk.Button(self.inital_sim_settings, text='Save results to:', command=self._select_target_folder, takefocus=False).grid(row=1, column=0, sticky="EW",)
    s5 = self._add_horisontal_scroll(self.inital_sim_settings, row=2, column=1)
    self._input_target = tk.Entry(self.inital_sim_settings, state='disabled', xscrollcommand=s5.set)
    self._input_target.grid(row=1, column=1, sticky="EW")
    s5.configure(command=self._input_target.xview)

    self.sim_max_water = tk.StringVar(value='3')
    tk.Label(self.inital_sim_settings, text='Maximum number of water molecules allowed in the bridge', anchor='w').grid(row=3, column=0, sticky='W')
    ttk.Combobox(self.inital_sim_settings, textvariable=self.sim_max_water, values=['1','2','3','4','5']).grid(row=3, column=1, sticky="EW")

    #--------------------------- dcd select------------

    self.selectSimFrame = ttk.LabelFrame(self.dcdframe, text='Select simulation')
    self.selectSimFrame.grid(self._crate_frame_grid(4))
    self.selectSimFrame.columnconfigure(0, weight=0)
    self.selectSimFrame.columnconfigure(1, weight=1)

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

    tk.Label(self.selectSimFrame, text='Name as: ', anchor='w').grid(row=9, column=0, sticky='W')
    # self.sim_name1 = tk.StringVar(value='sim1')
    self.sim_name = tk.Entry(self.selectSimFrame)
    self.sim_name.insert(0, 'sim1') # remove when test resolved
    self.sim_name.grid(row=9, column=1, sticky="EW")

    tk.Button(self.selectSimFrame, text='Construct graph', command=self._construct_sim_graphs).grid(self._create_big_button_grid(10, column=1))

    # ----------------------- DcdWaterWireFrame -----------------------
    self.DcdWaterWireFrame = ttk.LabelFrame(self.dcdframe, text='Water wire network')
    self.DcdWaterWireFrame.grid(self._crate_frame_grid(14))
    self.DcdWaterWireFrame.columnconfigure(0, weight=1)
    # self.DcdWaterWireFrame.columnconfigure(1, weight=1)
    # self.DcdWaterWireFrame.columnconfigure(2, weight=1)

    self.row=15
    tk.Button(self.DcdWaterWireFrame, text='Select graphs to compare', command=lambda:self._load_graph_files(self.row)).grid(self._create_big_button_grid(14))

    return
