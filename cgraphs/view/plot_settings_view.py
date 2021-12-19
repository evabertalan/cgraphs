import tkinter as tk
from tkinter import ttk



def plot_settings(self):
    main_frame = ttk.LabelFrame(self.plotsettings, text='Plot settings')
    main_frame.grid(self._crate_frame_grid(0))
    main_frame.columnconfigure(0, weight=1)
    main_frame.columnconfigure(1, weight=1)
    self.edge_with = tk.Entry(main_frame, bg='white', fg='black', highlightbackground='white', insertbackground='black', validate="key", validatecommand=(self.ifnum_cmd, '%S', '%P', 0, 20))
    self.edge_with.grid(row=1, column=1, sticky='W')
    # self.plot_parameters['edge_width'] = edge_with.get()
    tk.Button(main_frame, text='Save', bg='white', fg='black', command=self._save_plot_settings, takefocus=False,  highlightbackground='white').grid(row=2, column=2, sticky="EW",)


