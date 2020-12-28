import tkinter as tk
from tkinter import filedialog
import waterclusters as wc
import matplotlib.pyplot as plt


class View:
  def __init__(self, master):
    self.master = master

  def start(self):
    if hasattr(self, 'mainframe'):
      self._destroy_frame()

    self.master.title('Protein graph analyser')
    self.master.geometry('650x650')
    self._create_frame()
    tk.Label(self.mainframe, text='Please select a folder, which contains all the PDB files you want to include in the analysis').grid(row=0, column=0)
    tk.Button(self.mainframe, text='Select Folder', command=self._select_root_folder).grid(row=0, column=1)

    tk.Label(self.mainframe, text='Select refernece file').grid(row=2, column=0)
    tk.Button(self.mainframe, text='Select refernece file', command=self._select_reference_file).grid(row=2, column=1)

    tk.Label(self.mainframe, text='Water Cluster Analysis').grid(row=4, column=0)
    tk.Button(self.mainframe, text='Calculate water clustes', command=self._init_water_clusters).grid(row=4, column=1)

    tk.Label(self.mainframe, text='Conserved graph analysis').grid(row=6, column=0)
    tk.Button(self.mainframe, text='Calculate conserved H-bond network', command=self._init_conserved_graph_analysis).grid(row=6, column=1)

    useWaterCoords = tk.BooleanVar(self.mainframe, False)
    tk.Checkbutton(self.mainframe, text='Use water cluster coordinates', variable=useWaterCoords).grid(row=7, column=1)

    tk.Button(self.mainframe, text='Calculate conserved water wire network', command=self._init_conserved_graph_analysis).grid(row=6, column=1)

  def _select_root_folder(self):
    self.pdb_root_folder = filedialog.askdirectory(initialdir = "../")
    tk.Label(self.mainframe, text='Selected Folder:'+str(self.pdb_root_folder)).grid(row=1, column=0)

  def _select_reference_file(self):
    self.reference_pdb = filedialog.askopenfilename(initialdir = "../")
    tk.Label(self.mainframe, text='Selected Folder:'+str(self.reference_pdb)).grid(row=3, column=0)

  def _init_water_clusters(self):
    w = wc.WaterClusters(self.pdb_root_folder, reference_pdb=self.reference_pdb)
    tk.Label(self.mainframe, text='There are '+str(len(w.water_coordinates))+' water molecules in the '+str(len(w.superimposed_files))+' superimpsed files').grid(row=5, column=0)
    fig, ax = w.evaluate_parameters()
    plt.show()

  def _init_conserved_graph_analysis(self):
    pass




  def _destroy_frame(self):
    self.mainframe.destroy()

  def _create_frame(self):
    self.mainframe = tk.Frame(self.master)
    self.mainframe.place(relx=0.5, rely=0.5, anchor=tk.CENTER)


root = tk.Tk()
view = View(root)
view.start()
root.mainloop()
