import tkinter as tk
from tkinter import filedialog
import waterclusters as wc
import conservedgraph as cg


class View:
  def __init__(self, master):
    self.master = master
    self.ipad = 12
    self.pad = 6

  def start(self):
    if hasattr(self, 'mainframe'):
      self._destroy_frame()

    self.master.title('Protein graph analyser')
    self.master.geometry('550x550')
    self._create_frame()


    self.inputFrame = tk.LabelFrame(self.mainframe, text='Input Locations')
    self.inputFrame.grid(row=0, columnspan=3, sticky='EW', padx=(self.pad,self.pad), pady=(self.pad,self.pad), ipadx=self.ipad, ipady=self.ipad)
    self.inputFrame.columnconfigure(0, weight=1)
    self.inputFrame.columnconfigure(1, weight=1)

    tk.Button(self.inputFrame, text='Select PDB Folder', command=self._select_root_folder).grid(row=1, column=0, sticky="EW")
    s1 = self._add_horisontal_scroll(self.inputFrame, row=2, column=1)
    self._input_folder = tk.Entry(self.inputFrame, state='disabled', xscrollcommand=s1.set)
    self._input_folder.grid(row=1, column=1, sticky="EW", columnspan=2)
    s1.configure(command = self._input_folder.xview)

    tk.Button(self.inputFrame, text='Select refernece file', command=self._select_reference_file).grid(row=4, column=0, sticky="EW")
    s2 = self._add_horisontal_scroll(self.inputFrame, row=5, column=1)
    self._input_pdb = tk.Entry(self.inputFrame, state='disabled', xscrollcommand=s2.set)
    self._input_pdb.grid(row=4, column=1, sticky="EW")
    s2.configure(command = self._input_pdb.xview)


    self.waterClusterFrame = tk.LabelFrame(self.mainframe, text='Water Cluster Analysis')
    self.waterClusterFrame.grid(row=4, columnspan=3, sticky='EW', padx=(self.pad,self.pad), pady=(self.pad,self.pad), ipadx=self.ipad, ipady=self.ipad)
    self.waterClusterFrame.columnconfigure(0, weight=1)

    tk.Button(self.waterClusterFrame, text='Calculate water clustes', command=self._init_water_clusters).grid(row=4, column=0, padx=(self.pad,self.pad), pady=(self.pad,self.pad), sticky="EW")


    self.conservedNetworkFrame = tk.LabelFrame(self.mainframe, text='Conserved network analysis')
    self.conservedNetworkFrame.grid(row=6, columnspan=3, sticky='EW', padx=(self.pad,self.pad), pady=(self.pad,self.pad), ipadx=self.ipad, ipady=self.ipad)
    self.conservedNetworkFrame.columnconfigure(0, weight=1)
    self.conservedNetworkFrame.columnconfigure(1, weight=1)
    tk.Button(self.conservedNetworkFrame, text='Calculate conserved H-bond network', command=lambda:self._init_conserved_graph_analysis('hbond')).grid(row=7, column=0, padx=(self.pad,self.pad), pady=(self.pad,self.pad), sticky="EW")

    self.useWaterCoords = tk.BooleanVar()
    tk.Checkbutton(self.conservedNetworkFrame, text='Use water cluster coordinates', variable=self.useWaterCoords).grid(row=7, column=1, padx=(self.pad,self.pad), pady=(self.pad,self.pad))

    tk.Button(self.conservedNetworkFrame, text='Calculate conserved water wire network', command=lambda:self._init_conserved_graph_analysis('water_wire')).grid(row=8, column=0, padx=(self.pad,self.pad), pady=(self.pad,self.pad), sticky="EW")
    self.completed = tk.Label(self.conservedNetworkFrame, text='', fg='white')
    self.completed.grid(row=9, column=0)


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

    tk.Label(self.mainframe, text='All the generated files can be found in:\n'+self.pdb_root_folder+'/workfolder/\n\n Plots are located in:\n'+self.pdb_root_folder+'/workfolder/plots/', wraplength=520).grid(row=10, column=0, columnspan=3)

  def _init_water_clusters(self):
    w = wc.WaterClusters(self.pdb_root_folder, reference_pdb=self.reference_pdb)
    w.evaluate_parameters()
    w.calculate_cluster_centers()
    self.ref_coordinates = w.reference_coordinates
    tk.Label(self.waterClusterFrame, text='There are '+str(len(w.water_coordinates))+' water molecules in the '+str(len(w.superimposed_files))+' superimpsed files.\n The algorithm found '+str(w.n_clusters_)+' water clusters.').grid(row=5, column=0)

  def _init_conserved_graph_analysis(self, graph_type):
    self.completed.configure(text='', fg='white')
    if self.useWaterCoords.get(): _ref_coord = self.ref_coordinates
    else: _ref_coord=None
    c = cg.ConservedGraph(self.pdb_root_folder, reference_pdb=self.reference_pdb, reference_coordinates=_ref_coord)
    c.calculate_graphs(graph_type=graph_type)
    c.plot_graphs(label_nodes=True)
    c.plot_graphs(label_nodes=False)
    c.get_conserved_graph()
    c.plot_conserved_graph(label_nodes=True)
    c.plot_conserved_graph(label_nodes=False)
    c.plot_difference(label_nodes=True)
    c.plot_difference(label_nodes=False)
    self.completed.configure(text='Calculation completed', fg='green')

  def _add_horisontal_scroll(self, target, row=1, column=0):
    scroll = tk.Scrollbar(target, orient='horizontal')
    scroll.grid(row=row, column=column, sticky='EW')
    return scroll


  def on_enter(self, event, text):
        self.l2.configure(text=text)

  def on_leave(self, enter):
        self.l2.configure(text='')

  def _destroy_frame(self):
    self.mainframe.destroy()

  def _create_frame(self):
    self.mainframe = tk.Frame(self.master)
    self.mainframe.place(relx=0.5, rely=0.5, anchor=tk.CENTER)


root = tk.Tk()
view = View(root)
view.start()
root.mainloop()
