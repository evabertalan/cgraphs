import tkinter as tk
from tkinter import ttk

def compare_view(self):
    file_selector = ttk.LabelFrame(self.mainframe, text='Input locations')
    file_selector.grid(self._crate_frame_grid(0))

    # tk.Button(file_selector, text='PDB 1', command=select_pdb, takefocus=False).grid(row=1, column=0, sticky="EW",)
    # s5 = self._add_horisontal_scroll(self.inital_sim_settings, row=2, column=1)
    # self._input_target = tk.Entry(self.inital_sim_settings, state='disabled', xscrollcommand=s5.set)
    # self._input_target.grid(row=1, column=1, sticky="EW")
    # s5.configure(command=self._input_target.xview)
