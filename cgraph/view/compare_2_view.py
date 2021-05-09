import tkinter as tk
from tkinter import ttk

def compare_view(self):
    file_selector = ttk.LabelFrame(self.compframe, text='Input locations')
    file_selector.grid(self._crate_frame_grid(0), columnspan=3)
    file_selector.columnconfigure(0, weight=1)
    file_selector.columnconfigure(1, weight=1)

    tk.Button(file_selector, text='PDB 1', command=lambda:self._select_pdb1(pdb1_field), takefocus=False).grid(row=0, column=0, sticky="EW")
    s1 = self._add_horisontal_scroll(file_selector, row=1, column=1)
    pdb1_field = tk.Entry(file_selector, state='disabled', xscrollcommand=s1.set)
    pdb1_field.grid(row=0, column=1, sticky="EW")
    s1.configure(command=pdb1_field.xview)
    self.color1 = '#1b3ede'
    color_field1 = tk.Label(file_selector, width=2, bg=self.color1, anchor="w")
    color_field1.grid(row=0, column=2, sticky='W')
    color_field1.bind("<Button-1>", lambda x=self.color1, y=color_field1:self._choose_color1(x, y))

    tk.Button(file_selector, text='PDB 2', command=lambda:self._select_pdb2(pdb2_field)).grid(row=2, column=0, sticky="EW")
    s2 = self._add_horisontal_scroll(file_selector, row=3, column=1)
    pdb2_field = tk.Entry(file_selector, state='disabled', xscrollcommand=s2.set)
    pdb2_field.grid(row=2, column=1, sticky="EW")
    s2.configure(command=pdb2_field.xview)
    self.color2 = '#21c25f'
    color_field2 = tk.Label(file_selector, width=2, bg=self.color2, anchor="w")
    color_field2.grid(row=2, column=2, sticky='W')
    color_field2.bind("<Button-1>", lambda x=self.color2, y=color_field2:self._choose_color2(x, y))

    tk.Button(file_selector, text='Location of workfolder', command=lambda:self._set_compare_result_folder(compare_results_input), takefocus=False).grid(row=4, column=0, sticky="EW")
    s3 = self._add_horisontal_scroll(file_selector, row=5, column=1)
    compare_results_input = tk.Entry(file_selector, state='disabled', xscrollcommand=s3.set)
    compare_results_input.grid(row=4, column=1, sticky="EW")
    s3.configure(command=compare_results_input.xview)

    hcritera_frame = tk.Frame(self.compframe)
    hcritera_frame.grid(row=6, column=0, columnspan=5, sticky="EW")
    self.comp_distance = tk.DoubleVar(value=3.5)
    self.comp_cut_angle = tk.DoubleVar(value=60)
    tk.Label(hcritera_frame, text='  H-bond cut-off criteria ', anchor="w").grid(row=6, column=0, sticky='E')
    ttk.Spinbox(hcritera_frame, textvariable=self.c_distance, from_=0, to=5, width=11, validate="key", validatecommand=(self.ifnum_cmd, '%S', '%P', 0, 5)).grid(row=6, column=1, sticky='W')
    tk.Label(hcritera_frame, text='Å distance and ', anchor="w").grid(row=6, column=2, sticky='W')
    ttk.Spinbox(hcritera_frame, textvariable=self.c_cut_angle, from_=0, to=180, width=11, validate="key", validatecommand=(self.ifnum_cmd, '%S', '%P', 0, 180)).grid(row=6, column=3, sticky='W')
    tk.Label(hcritera_frame, text='degrees angle', anchor="w").grid(row=6, column=4, sticky='W')

    hbond_frame = ttk.LabelFrame(self.compframe, text='H-bond network')
    hbond_frame.grid(self._crate_frame_grid(7))
    hbond_frame.columnconfigure(0, weight=1)
    hbond_frame.columnconfigure(1, weight=2)

    self.include_waters_comp = tk.BooleanVar()
    self.include_waters_comp.set(True)
    tk.Checkbutton(hbond_frame, text='Include crystallographic waters', variable=self.include_waters_comp, anchor="w").grid(self._create_big_button_grid(8))

    self.include_backbone_sidechain_comp = tk.BooleanVar()
    tk.Checkbutton(hbond_frame, text='Include sidechain-backbone interactions', variable=self.include_backbone_sidechain_comp, anchor="w").grid(self._create_big_button_grid(9))

    tk.Button(hbond_frame, text='Compare H-bond networks', command=lambda:self._init_pdb_comparison('hbond', pdb1=self.pdb_1, pdb2=self.pdb_2, color1=self.color1, color2=self.color2), width=self.button_width).grid(self._create_big_button_grid(10))


    # -------------------water_wire_frame -----------------------

    water_wire_frame = ttk.LabelFrame(self.compframe, text='Water wire network')
    water_wire_frame.grid(self._crate_frame_grid(11))
    water_wire_frame.columnconfigure(0, weight=1)
    water_wire_frame.columnconfigure(1, weight=1)

    self.max_water_comp = tk.IntVar(value=3)
    tk.Label(water_wire_frame, text='Maximum number of water molecules allowed in the bridge', anchor="w").grid(row=11, column=0, sticky='W')
    ttk.Combobox(water_wire_frame, textvariable=self.max_water_comp, values=[1,2,3,4,5], state='readonly').grid(row=11, column=1, sticky="EW")
    tk.Button(water_wire_frame, text='Compare water wire network', command=lambda:self._init_pdb_comparison('water_wire',pdb1=self.pdb_1, pdb2=self.pdb_2, color1=self.color1, color2=self.color2), width=self.button_width).grid(self._create_big_button_grid(12))
