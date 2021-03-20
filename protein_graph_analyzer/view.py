import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from .waterclusters import WaterClusters
from .conservedgraph import ConservedGraph


class View:
  def __init__(self, master):
    self.master = master
    self.ipad = 12
    self.pad = 6
    self.button_width = 15

  def main_modal(self):
    if hasattr(self, 'mainframe'):
      self._destroy_frame()

    self.master.title('Protein graph analyser')
    self.master.geometry('800x750')
    self._create_frame()

    self.inputFrame = tk.LabelFrame(self.mainframe, text='Input locations')
    self.inputFrame.grid(row=0, columnspan=3, sticky='EW', padx=(self.pad,self.pad), pady=(self.pad,self.pad), ipadx=self.ipad, ipady=self.ipad)
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
    self.waterClusterFrame.grid(row=7, columnspan=3, sticky='EW', padx=(self.pad,self.pad), pady=(self.pad,self.pad), ipadx=self.ipad, ipady=self.ipad)
    self.waterClusterFrame.columnconfigure(0, weight=1)

    tk.Button(self.waterClusterFrame, text='Calculate water clusters', command=self._init_water_clusters, width=self.button_width).grid(row=4, column=0, padx=(self.pad,self.pad), pady=(self.pad,self.pad), sticky="EW")


    # ----------------------- conservedNetworkFrame -----------------------

    self.conservedNetworkFrame = tk.LabelFrame(self.mainframe, text='Conserved network analysis')
    self.conservedNetworkFrame.grid(row=8, columnspan=3, sticky='EW', padx=(self.pad,self.pad), pady=(self.pad,self.pad), ipadx=self.ipad, ipady=self.ipad)
    self.conservedNetworkFrame.columnconfigure(0, weight=1)
    self.conservedNetworkFrame.columnconfigure(1, weight=1)

    self.conservation_threshold = tk.StringVar(value='90')
    tk.Label(self.conservedNetworkFrame, text='Conservation of H-bonding groups across structures', anchor="w").grid(row=9, column=0)
    ttk.Spinbox(self.conservedNetworkFrame, textvariable=self.conservation_threshold, from_=1, to=100).grid(row=9, column=1, sticky="EW")

    # ----------------------- HbondNetworkFrame -----------------------

    self.HbondNetworkFrame = tk.LabelFrame(self.conservedNetworkFrame, text='H-bond network')
    self.HbondNetworkFrame.grid(row=10, columnspan=3, sticky='EW', padx=(self.pad,self.pad), pady=(self.pad,self.pad), ipadx=self.ipad, ipady=self.ipad)
    self.HbondNetworkFrame.columnconfigure(0, weight=1)
    self.HbondNetworkFrame.columnconfigure(1, weight=1)

    self.useWaterCoords = tk.BooleanVar()
    tk.Checkbutton(self.HbondNetworkFrame, text='Use water cluster coordinates', variable=self.useWaterCoords, anchor="w").grid(row=10, column=0, padx=(self.pad,self.pad), pady=(self.pad,self.pad), sticky="EW")

    self.include_backbone_sidechain = tk.BooleanVar()
    tk.Checkbutton(self.HbondNetworkFrame, text='Include sidechain-backbone interactions', variable=self.include_backbone_sidechain, anchor="w").grid(row=11, column=0, padx=(self.pad,self.pad), pady=(self.pad,self.pad), sticky="EW")

    self.include_backbone_backbone = tk.BooleanVar()
    tk.Checkbutton(self.HbondNetworkFrame, text='Include backbone-backbone interactions', variable=self.include_backbone_backbone, anchor="w").grid(row=12, column=0, padx=(self.pad,self.pad), pady=(self.pad,self.pad), sticky="EW")

    tk.Button(self.HbondNetworkFrame, text='Calculate conserved H-bond network', command=lambda:self._init_conserved_graph_analysis('hbond'), width=self.button_width).grid(row=13, column=0, padx=(self.pad,self.pad), pady=(self.pad,self.pad), sticky="EW")


    # ----------------------- WaterWireFrame -----------------------

    self.WaterWireFrame = tk.LabelFrame(self.conservedNetworkFrame, text='Water wire network')
    self.WaterWireFrame.grid(row=14, columnspan=3, sticky='EW', padx=(self.pad,self.pad), pady=(self.pad,self.pad), ipadx=self.ipad, ipady=self.ipad)
    self.WaterWireFrame.columnconfigure(0, weight=1)
    self.WaterWireFrame.columnconfigure(1, weight=1)

    self.max_water = tk.StringVar(value='3')
    tk.Label(self.WaterWireFrame, text='Maximum number of water molecules allowed in the bridge', anchor="w").grid(row=15, column=0)
    ttk.Combobox(self.WaterWireFrame, textvariable=self.max_water, values=['1','2','3','4','5']).grid(row=15, column=1, sticky="EW")
    tk.Button(self.WaterWireFrame, text='Calculate conserved water wire network', command=lambda:self._init_conserved_graph_analysis('water_wire'), width=self.button_width).grid(row=16, column=0, padx=(self.pad,self.pad), pady=(self.pad,self.pad), sticky="EW")

    self.completedText = tk.StringVar()
    self.completed = tk.Label(self.conservedNetworkFrame, textvariable=self.completedText)
    self.completed.grid(row=17, column=0)

  def _select_root_folder(self):
    self.pdb_root_folder = filedialog.askdirectory(initialdir = "../")
    self._input_folder.configure(state='normal')
    self._input_folder.insert(0, str(self.pdb_root_folder))
    self._input_folder.configure(state='disabled')


  def _select_reference_file(self):
    self.reference_pdb = filedialog.askopenfilename(initialdir = "../")
    self._input_pdb.configure(state='normal')
    self._input_pdb.insert(0, str(self.reference_pdb))
    self._input_pdb.configure(state='disabled')

    # tk.Label(self.mainframe, anchor='w', text='All the generated files can be found in:\n'+self.pdb_root_folder+'/workfolder/\n\n Plots are located in:\n'+self.pdb_root_folder+'/workfolder/plots/', wraplength=520).grid(row=10, column=0, columnspan=3)

  def _init_water_clusters(self):
    sst = int(self.sequance_identity_threshold.get())/100
    # valudate sst
    w = WaterClusters(self.pdb_root_folder, reference_pdb=self.reference_pdb, sequance_identity_threshold=sst)
    w.evaluate_parameters()
    w.calculate_cluster_centers()
    w.write_cluster_center_coordinates()
    w.draw_clusters_centers_chimera()
    self.ref_coordinates = w.reference_coordinates
    tk.Label(self.waterClusterFrame, text='There are '+str(len(w.water_coordinates))+' water molecules in the '+str(len(w.superimposed_files))+' uperimposed files.\n The algorithm found '+str(w.n_clusters_)+' water clusters.').grid(row=5, column=0)
    w.logger.info('Water cluster calculation is completed\n'+'-'*20)

  def _init_conserved_graph_analysis(self, graph_type):
    self._update_lable_text('')
    ebb = not self.include_backbone_backbone.get()
    ieb = self.include_backbone_sidechain.get()
    if self.useWaterCoords.get(): _ref_coord = self.ref_coordinates
    else: _ref_coord=None
    c = ConservedGraph(self.pdb_root_folder, reference_pdb=self.reference_pdb, reference_coordinates=_ref_coord)
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

  def _destroy_frame(self):
    self.mainframe.destroy()

  def _create_frame(self):
    self.mainframe = tk.Frame(self.master)
    self.mainframe.place(relx=0.5, rely=0.5, anchor=tk.CENTER)

def start():
    root = tk.Tk()
    view = View(root)
    view.main_modal()
    root.mainloop()
