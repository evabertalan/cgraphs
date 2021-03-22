import tkinter as tk
from tkinter import ttk

def ta_view(self):

    self.selectSimFrame = tk.LabelFrame(self.dcdframe, text='Select simulations')
    self.selectSimFrame.grid(row=0, columnspan=4, sticky='EW', padx=(self.padx,self.padx), pady=(self.pady,self.pady), ipadx=self.ipadx, ipady=self.ipady)
    self.selectSimFrame.columnconfigure(0, weight=1)
    self.selectSimFrame.columnconfigure(1, weight=1)
    self.selectSimFrame.columnconfigure(2, weight=1)
    self.selectSimFrame.columnconfigure(3, weight=1)

    tk.Button(self.selectSimFrame, text='Select PSF', command=self._select_psf_file).grid(row=1, column=0, sticky="EW", columnspan=1)
    s3 = self._add_horisontal_scroll(self.selectSimFrame, row=2, column=1)
    self._input_psf = tk.Entry(self.selectSimFrame, state='disabled', xscrollcommand=s3.set)
    self._input_psf.grid(row=1, column=1, sticky="EW", columnspan=3)
    s3.configure(command=self._input_psf.xview)

    tk.Button(self.selectSimFrame, text='Select DCDs', command=self._select_dcd_files).grid(row=3, column=0, sticky="EW", columnspan=1)
    s4 = self._add_horisontal_scroll(self.selectSimFrame, row=4, column=1,)
    self._input_dcd = tk.Entry(self.selectSimFrame, state='disabled', xscrollcommand=s4.set)
    self._input_dcd.grid(row=3, column=1, sticky="EW", columnspan=3)
    s4.configure(command=self._input_dcd.xview)

    tk.Label(self.selectSimFrame, text='AS ').grid(row=4, column=3, sticky="EW")
    self.sim_name = tk.StringVar(value='sim1')
    tk.Entry(self.selectSimFrame, textvariable=self.sim_name).grid(row=4, column=4, sticky="EW")

    tk.Button(self.selectSimFrame, text='Constract graph', command=self._constract_sim_graphs).grid(row=5, column=0, padx=(self.padx,self.padx), pady=(self.pady,self.pady), sticky="EW")
