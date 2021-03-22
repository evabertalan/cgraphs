import tkinter as tk
from tkinter import ttk

# from tkinter import filedialog


def csa_view(self):



    # ----------------------- inputFrame -----------------------

    self.inputFrame = tk.LabelFrame(self.mainframe, text='Input locations')
    self.inputFrame.grid(row=0, columnspan=3, sticky='EW', padx=(self.padx,self.padx), pady=(self.pady,self.pady), ipadx=self.ipadx, ipady=self.ipady)
    self.inputFrame.columnconfigure(0, weight=1)
    self.inputFrame.columnconfigure(1, weight=1)

    tk.Button(self.inputFrame, text='Select PDB folder', command=self._select_root_folder, width=self.button_width).grid(row=1, column=0, sticky="EW")
    s1 = self._add_horisontal_scroll(self.inputFrame, row=2, column=1)
    self._input_folder = tk.Entry(self.inputFrame, state='disabled', xscrollcommand=s1.set)
    self._input_folder.grid(row=1, column=1, sticky="EW", columnspan=2)
    s1.configure(command=self._input_folder.xview)

    tk.Button(self.inputFrame, text='Select reference file', command=self._select_reference_file, width=self.button_width).grid(row=4, column=0, sticky="EW")
    s2 = self._add_horisontal_scroll(self.inputFrame, row=5, column=1)
    self._input_pdb = tk.Entry(self.inputFrame, state='disabled', xscrollcommand=s2.set)
    self._input_pdb.grid(row=4, column=1, sticky="EW")
    s2.configure(command=self._input_pdb.xview)

    self.sequance_identity_threshold = tk.StringVar(value='75')
    tk.Label(self.inputFrame, text='Minimum sequence identity (%)').grid(row=6, column=0)
    ttk.Spinbox(self.inputFrame, textvariable=self.sequance_identity_threshold, from_=1, to=100).grid(row=6, column=1, sticky="EW")


    # ----------------------- waterClusterFrame -----------------------

    self.waterClusterFrame = tk.LabelFrame(self.mainframe, text='Water cluster analysis')
    self.waterClusterFrame.grid(row=7, columnspan=3, sticky='EW', padx=(self.padx,self.padx), pady=(self.pady,self.pady), ipadx=self.ipadx, ipady=self.ipady)
    self.waterClusterFrame.columnconfigure(0, weight=1)

    tk.Button(self.waterClusterFrame, text='Parameter analysis', command=self._perform_parameter_analysis, width=self.button_width).grid(row=0, column=0, padx=(self.padx,self.padx), pady=(self.pady,self.pady), sticky="EW")

    tk.Label(self.waterClusterFrame, text='DBSCAN eps:').grid(row=0, column=1, sticky="EW")
    self.eps = tk.StringVar(value='1.4')
    tk.Entry(self.waterClusterFrame, textvariable=self.eps).grid(row=0, column=2, sticky="EW")

    tk.Button(self.waterClusterFrame, text='Calculate water clusters', command=self._init_water_clusters, width=self.button_width).grid(row=2, column=0, padx=(self.padx,self.padx), pady=(self.pady,self.pady), sticky="EW")


    # ----------------------- conservedNetworkFrame -----------------------

    self.conservedNetworkFrame = tk.LabelFrame(self.mainframe, text='Conserved network analysis')
    self.conservedNetworkFrame.grid(row=8, columnspan=3, sticky='EW', padx=(self.padx,self.padx), pady=(self.pady,self.pady), ipadx=self.ipadx, ipady=self.ipady)
    self.conservedNetworkFrame.columnconfigure(0, weight=1)
    self.conservedNetworkFrame.columnconfigure(1, weight=1)

    self.conservation_threshold = tk.StringVar(value='90')
    tk.Label(self.conservedNetworkFrame, text='Conservation of H-bonding groups across structures', anchor="w").grid(row=9, column=0)
    ttk.Spinbox(self.conservedNetworkFrame, textvariable=self.conservation_threshold, from_=1, to=100).grid(row=9, column=1, sticky="EW")

    # ----------------------- HbondNetworkFrame -----------------------

    self.HbondNetworkFrame = tk.LabelFrame(self.conservedNetworkFrame, text='H-bond network')
    self.HbondNetworkFrame.grid(row=10, columnspan=3, sticky='EW', padx=(self.padx,self.padx), pady=(self.pady,self.pady), ipadx=self.ipadx, ipady=self.ipady)
    self.HbondNetworkFrame.columnconfigure(0, weight=1)
    self.HbondNetworkFrame.columnconfigure(1, weight=1)

    self.useWaterCoords = tk.BooleanVar()
    tk.Checkbutton(self.HbondNetworkFrame, text='Use water cluster coordinates', variable=self.useWaterCoords, anchor="w").grid(row=10, column=0, padx=(self.padx,self.padx), pady=(self.pady,self.pady), sticky="EW")

    self.include_backbone_sidechain = tk.BooleanVar()
    tk.Checkbutton(self.HbondNetworkFrame, text='Include sidechain-backbone interactions', variable=self.include_backbone_sidechain, anchor="w").grid(row=11, column=0, padx=(self.padx,self.padx), pady=(self.pady,self.pady), sticky="EW")

    # self.include_backbone_backbone = tk.BooleanVar()
    # tk.Checkbutton(self.HbondNetworkFrame, text='Include backbone-backbone interactions', variable=self.include_backbone_backbone, anchor="w").grid(row=12, column=0, padx=(self.padx,self.padx), pady=(self.pady,self.pady), sticky="EW")

    tk.Button(self.HbondNetworkFrame, text='Calculate conserved H-bond network', command=lambda:self._init_conserved_graph_analysis('hbond'), width=self.button_width).grid(row=13, column=0, padx=(self.padx,self.padx), pady=(self.pady,self.pady), sticky="EW")


    # ----------------------- WaterWireFrame -----------------------

    self.WaterWireFrame = tk.LabelFrame(self.conservedNetworkFrame, text='Water wire network')
    self.WaterWireFrame.grid(row=14, columnspan=3, sticky='EW', padx=(self.padx,self.padx), pady=(self.pady,self.pady), ipadx=self.ipadx, ipady=self.ipady)
    self.WaterWireFrame.columnconfigure(0, weight=1)
    self.WaterWireFrame.columnconfigure(1, weight=1)

    self.max_water = tk.StringVar(value='3')
    tk.Label(self.WaterWireFrame, text='Maximum number of water molecules allowed in the bridge', anchor="w").grid(row=15, column=0)
    ttk.Combobox(self.WaterWireFrame, textvariable=self.max_water, values=['1','2','3','4','5']).grid(row=15, column=1, sticky="EW")
    tk.Button(self.WaterWireFrame, text='Calculate conserved water wire network', command=lambda:self._init_conserved_graph_analysis('water_wire'), width=self.button_width).grid(row=16, column=0, padx=(self.padx,self.padx), pady=(self.pady,self.pady), sticky="EW")

    self.completedText = tk.StringVar()
    self.completed = tk.Label(self.conservedNetworkFrame, textvariable=self.completedText)
    self.completed.grid(row=17, column=0)

    return
