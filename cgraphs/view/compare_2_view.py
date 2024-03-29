import tkinter as tk
from tkinter import ttk

def compare_view(self):
    file_selector = ttk.LabelFrame(self.compframe, text='Parameter set up for the analysis')
    file_selector.grid(self._crate_frame_grid(0), columnspan=3)
    file_selector.columnconfigure(0, weight=1)
    file_selector.columnconfigure(1, weight=1)

    tk.Button(file_selector, text='Location of workfolder', command=lambda:self._set_compare_result_folder(compare_results_input), takefocus=False, bg='white', fg='black', highlightbackground='white').grid(row=0, column=0, sticky="EW")
    s3 = self._add_horisontal_scroll(file_selector, row=1, column=1)
    compare_results_input = tk.Entry(file_selector, state='disabled', xscrollcommand=s3.set, highlightbackground='white', bg='white', fg='black', disabledbackground=self.gray, disabledforeground='black')
    compare_results_input.grid(row=0, column=1, sticky="EW")
    s3.configure(command=compare_results_input.xview, bg='white')

    tk.Label(file_selector, text='  H-bond criteria ', anchor="w", bg='white', fg='black').grid(row=2, column=0, sticky='W')
    comp_hcritera_frame = tk.Frame(file_selector,bg='white')
    comp_hcritera_frame.grid(row=2, column=1, columnspan=4, sticky="EW")
    self.comp_distance = tk.DoubleVar(value=3.5)
    self.comp_cut_angle = tk.DoubleVar(value=60)
    ttk.Spinbox(comp_hcritera_frame, textvariable=self.comp_distance, from_=0, to=5, width=5, validate="key", validatecommand=(self.ifnum_cmd, '%S', '%P', 0, 5)).grid(row=0, column=1, sticky='W')
    tk.Label(comp_hcritera_frame, text='Å distance     ', anchor="w", bg='white', fg='black').grid(row=0, column=2, sticky='W')
    self.comp_use_angle = tk.BooleanVar()
    self.comp_use_angle.set(False)
    tk.Checkbutton(comp_hcritera_frame, variable=self.comp_use_angle, anchor="w", bg='white', fg='black').grid(row=0, column=3, sticky='E')
    ttk.Spinbox(comp_hcritera_frame, textvariable=self.comp_cut_angle, from_=0, to=180, width=5, validate="key", validatecommand=(self.ifnum_cmd, '%S', '%P', 0, 180)).grid(row=0, column=4, sticky='W')
    tk.Label(comp_hcritera_frame, text='degrees angle', anchor="w", bg='white', fg='black').grid(row=0, column=5, sticky='W')

    comparison_frame = ttk.LabelFrame(self.compframe, text='Compare networks')
    comparison_frame.grid(self._crate_frame_grid(2))
    compare_tab = ttk.Notebook(comparison_frame)
    pdb_compare_tab = ttk.Frame(compare_tab)
    pdb_compare_tab.grid(columnspan=3)
    self.dcd_compare_tab = ttk.Frame(compare_tab)
    compare_tab.add(pdb_compare_tab, text='PDB crystal structures', padding=(30, 8, 30, 8))
    compare_tab.add(self.dcd_compare_tab, text='MD trajectories', padding=(30, 8, 30, 8))
    compare_tab.grid(self._crate_frame_grid(0))
    pdb_compare_tab.columnconfigure(0, weight=1)
    pdb_compare_tab.columnconfigure(1, weight=1)
    self.dcd_compare_tab.columnconfigure(0, weight=1)
    self.dcd_compare_tab.columnconfigure(1, weight=1)

    # ------------------- PDB tab -----------------------

    tk.Button(pdb_compare_tab, text='PDB 1', command=lambda:self._select_pdb1(pdb1_field), takefocus=False, bg='white', fg='black', highlightbackground='white').grid(row=2, column=0, sticky="EW")
    s1 = self._add_horisontal_scroll(pdb_compare_tab, row=3, column=1)
    pdb1_field = tk.Entry(pdb_compare_tab, state='disabled', xscrollcommand=s1.set, bg='white', fg='black', highlightbackground='white', disabledbackground=self.gray, disabledforeground='black')
    pdb1_field.grid(row=2, column=1, sticky="EW")
    s1.configure(command=pdb1_field.xview, bg='white')
    self.color1 = '#1b3ede'
    color_field1 = tk.Label(pdb_compare_tab, width=2, bg=self.color1, anchor="w")
    color_field1.grid(row=2, column=2, sticky='W')
    color_field1.bind("<Button-1>", lambda x=self.color1, y=color_field1:self._choose_color1(x, y))

    tk.Button(pdb_compare_tab, text='PDB 2', command=lambda:self._select_pdb2(pdb2_field), bg='white', fg='black',highlightbackground='white').grid(row=4, column=0, sticky="EW")
    s2 = self._add_horisontal_scroll(pdb_compare_tab, row=5, column=1)
    pdb2_field = tk.Entry(pdb_compare_tab, state='disabled', xscrollcommand=s2.set, bg='white', fg='black', highlightbackground='white', disabledbackground=self.gray, disabledforeground='black')
    pdb2_field.grid(row=4, column=1, sticky="EW")
    s2.configure(command=pdb2_field.xview, bg='white')
    self.color2 = '#21c25f'
    color_field2 = tk.Label(pdb_compare_tab, width=2, bg=self.color2, anchor="w")
    color_field2.grid(row=4, column=2, sticky='W')
    color_field2.bind("<Button-1>", lambda x=self.color2, y=color_field2:self._choose_color2(x, y))

    selsting_frame = tk.Frame(pdb_compare_tab, bg='white')
    selsting_frame.grid(row=6, column=0, columnspan=2, sticky="EW")
    selsting_frame.columnconfigure(1, weight=1)
    self.pdb_comp_selection_string, self.pdb_comp_selected_donors_pdb, self.pdb_comp_selected_acceptors_pdb = self.custom_selection_strin(selsting_frame, 6)

    # # -------------------hbond -----------------------

    hbond_frame = ttk.LabelFrame(pdb_compare_tab, text='H-bond network')
    hbond_frame.grid(self._crate_frame_grid(7))
    hbond_frame.columnconfigure(0, weight=1)
    hbond_frame.columnconfigure(1, weight=2)

    self.include_waters_comp = tk.BooleanVar()
    self.include_waters_comp.set(True)
    tk.Checkbutton(hbond_frame, text='Include crystallographic waters', variable=self.include_waters_comp, anchor="w", bg='white', fg='black').grid(self._create_big_button_grid(8))

    self.include_backbone_sidechain_comp = tk.BooleanVar()
    tk.Checkbutton(hbond_frame, text='Include sidechain-backbone interactions', variable=self.include_backbone_sidechain_comp, anchor="w", bg='white', fg='black').grid(self._create_big_button_grid(9))

    tk.Button(hbond_frame, text='Compare H-bond networks', command=lambda:self._init_pdb_comparison(comp_type='hbond', pdb1=self.pdb_1, pdb2=self.pdb_2, color1=self.color1, color2=self.color2), width=self.button_width, bg='white', fg='black',highlightbackground='white').grid(self._create_big_button_grid(10), columnspan=2)


    # # -------------------water_wire_frame -----------------------

    water_wire_frame = ttk.LabelFrame(pdb_compare_tab, text='Water wire network')
    water_wire_frame.grid(self._crate_frame_grid(11))
    water_wire_frame.columnconfigure(0, weight=1)
    water_wire_frame.columnconfigure(1, weight=1)

    self.max_water_comp = tk.IntVar(value=3)
    tk.Label(water_wire_frame, text='Maximum number of waters in the bridge', anchor="w", bg='white', fg='black').grid(row=11, column=0, sticky='W')
    ttk.Combobox(water_wire_frame, textvariable=self.max_water_comp, values=[1,2,3,4,5], state='readonly').grid(row=11, column=1, sticky="EW")
    tk.Button(water_wire_frame, text='Compare water wire network', command=lambda:self._init_pdb_comparison(comp_type='water_wire' ,pdb1=self.pdb_1, pdb2=self.pdb_2, color1=self.color1, color2=self.color2), width=self.button_width, bg='white', fg='black', highlightbackground='white').grid(self._create_big_button_grid(12), columnspan=2)

    # # ------------------- TRAJECTORY tab -----------------------

    tk.Button(self.dcd_compare_tab, text='PSF 1', command=lambda:self._select_psf1_compare(psf1_field), takefocus=False, bg='white', fg='black', highlightbackground='white').grid(row=2, column=0, sticky="EW")
    s3 = self._add_horisontal_scroll(self.dcd_compare_tab, row=3, column=1)
    psf1_field = tk.Entry(self.dcd_compare_tab, state='disabled', xscrollcommand=s3.set,  highlightbackground='white', bg='white', fg='black', disabledbackground=self.gray, disabledforeground='black')
    psf1_field.grid(row=2, column=1, sticky="EW")
    s3.configure(command=psf1_field.xview, bg='white')

    tk.Button(self.dcd_compare_tab, text='DCDs 1', command=lambda:self._select_dcd1_compare(dcd1_field), bg='white', fg='black', highlightbackground='white').grid(row=4, column=0, sticky="EW")
    s4 = self._add_horisontal_scroll(self.dcd_compare_tab, row=5, column=1)
    dcd1_field = tk.Entry(self.dcd_compare_tab, state='disabled', xscrollcommand=s4.set, highlightbackground='white', bg='white', fg='black', disabledbackground=self.gray, disabledforeground='black')
    dcd1_field.grid(row=4, column=1, sticky="EW")
    s4.configure(command=dcd1_field.xview, bg='white')

    tk.Label(self.dcd_compare_tab, text='Name 1: ', anchor='w', bg='white', fg='black').grid(row=6, column=0, sticky='W')
    # self.sim_name1 = tk.StringVar(value='sim1')
    field1 = tk.Frame(self.dcd_compare_tab, bg='white')
    field1.grid(row=6, column=1, sticky="EW", pady=(0,12))
    field1.columnconfigure(0, weight=1)
    self.compare_dcd1_name = tk.Entry(field1, bg='white', fg='black', highlightbackground='white', insertbackground='black')
    self.compare_dcd1_name.insert(0, 'sim1') # remove when test resolved
    self.compare_dcd1_name.grid(row=6, column=0, sticky="EW")
    self.color_dcd1 = '#1b3ede'
    color_dcd_field1 = tk.Label(field1, width=2, bg=self.color_dcd1, anchor="w")
    color_dcd_field1.grid(row=6, column=1, sticky='W')
    color_dcd_field1.bind("<Button-1>", lambda x=self.color_dcd1, y=color_dcd_field1:self._choose_dcd_color1(x, y))

    tk.Button(self.dcd_compare_tab, text='PSF 2', command=lambda:self._select_psf2_compare(psf2_field), bg='white', fg='black', highlightbackground='white').grid(row=7, column=0, sticky="EW")
    s5 = self._add_horisontal_scroll(self.dcd_compare_tab, row=8, column=1)
    psf2_field = tk.Entry(self.dcd_compare_tab, state='disabled', xscrollcommand=s5.set, bg='white', fg='black', highlightbackground='white', disabledbackground=self.gray, disabledforeground='black')
    psf2_field.grid(row=7, column=1, sticky="EW")
    s5.configure(command=psf2_field.xview, bg='white')

    tk.Button(self.dcd_compare_tab, text='DCDs 2', command=lambda:self._select_dcd2_compare(dcd2_field), bg='white', fg='black', highlightbackground='white').grid(row=9, column=0, sticky="EW")
    s6 = self._add_horisontal_scroll(self.dcd_compare_tab, row=10, column=1)
    dcd2_field = tk.Entry(self.dcd_compare_tab, state='disabled', xscrollcommand=s6.set, bg='white', fg='black', highlightbackground='white', disabledbackground=self.gray, disabledforeground='black')
    dcd2_field.grid(row=9, column=1, sticky="EW")
    s6.configure(command=dcd2_field.xview, bg='white')

    tk.Label(self.dcd_compare_tab, text='Name 2: ', anchor='w', bg='white', fg='black').grid(row=11, column=0, sticky='W')
    # self.sim_name1 = tk.StringVar(value='sim1')
    field2 = tk.Frame(self.dcd_compare_tab, bg='white')
    field2.grid(row=11, column=1, sticky="EW")
    field2.columnconfigure(0, weight=1)
    self.compare_dcd2_name = tk.Entry(field2, bg='white', fg='black', highlightbackground='white',  insertbackground='black')
    self.compare_dcd2_name.insert(0, 'sim2') # remove when test resolved
    self.compare_dcd2_name.grid(row=11, column=0, sticky="EW")
    l1 = tk.Label(self.dcd_compare_tab, text='Each simulation must be given a unique name, else the results will be overwritten.', anchor='w', justify='left', bg='white', fg='black')
    l1.grid(row=12, column=1, sticky='EW', pady=(0,24))
    l1.config(font=("Helvetica", 11))
    self.color_dcd2 = '#21c25f'
    color_dcd_field2 = tk.Label(field2, width=2, bg=self.color_dcd2, anchor="w")
    color_dcd_field2.grid(row=11, column=1, sticky='W')
    color_dcd_field2.bind("<Button-1>", lambda x=self.color_dcd2, y=color_dcd_field2:self._choose_dcd_color2(x, y))

    self.max_water_comp_dcd = tk.IntVar(value=3)
    tk.Label(self.dcd_compare_tab, text='Maximum number of waters in the bridge', anchor="w", bg='white', fg='black').grid(row=13, column=0, sticky='W')
    ttk.Combobox(self.dcd_compare_tab, textvariable=self.max_water_comp_dcd, values=[1,2,3,4,5], state='readonly').grid(row=13, column=1, sticky="EW")
    selsting_frame = tk.Frame(self.dcd_compare_tab, bg='white')
    selsting_frame.grid(row=14, column=0, columnspan=2, sticky="EW")
    selsting_frame.columnconfigure(1, weight=1)
    self.dcd_comp_selection_string, self.dcd_comp_selected_donors_pdb, self.dcd_comp_selected_acceptors_pdb = self.custom_selection_strin(selsting_frame, 1)

    tk.Button(self.dcd_compare_tab, text='Construct graph', command=lambda:self._construct_compare_graphs(psf1=self.psf_1, psf2=self.psf_2, dcd1=self.dcd_1, dcd2=self.dcd_2), bg='white', fg='black', highlightbackground='white').grid(self._create_big_button_grid(15), columnspan=2)

    self.compare_row = 15
