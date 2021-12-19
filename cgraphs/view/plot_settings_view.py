import tkinter as tk
from tkinter import ttk



def plot_settings(self):
    def _create_entries(self, frame, row, label, key):

        tk.Label(frame, text=label, anchor='w',  bg='white', fg='black').grid(row=row, column=0, sticky='W')
        var = tk.Entry(frame, bg='white', fg='black', highlightbackground='white', insertbackground='black',
            validate="key", validatecommand=(self.ifnum_cmd, '%S', '%P', 0, 20))
        var.grid(row=row, column=1, sticky='W')
        print(key, var)
        setattr(self, key, var)

    main_frame = ttk.LabelFrame(self.plotsettings, text='Plot settings')
    main_frame.grid(self._crate_frame_grid(0))
    main_frame.columnconfigure(0, weight=1)
    main_frame.columnconfigure(1, weight=1)


    self.edge_width = None
    self.node_label_size= None
    self.edge_label_size= None
    self.node_size= None
    self.node_color= None
    self.water_node_color= None
    self.edge_color= None
    self.plot_title_fontsize= None
    self.plot_label_fontsize= None
    self.plot_tick_fontsize= None
    self.plot_resolution= None

    settings_options = {
                'edge_width': {
                'var': self.edge_width,
                'label': 'Edge width',
                } ,
                'node_label_size': {
                'var': self.node_label_size,
                'label': 'Node label size',
                } ,
                'edge_label_size': {
                'var': self.edge_label_size,
                'label': 'Edge label size',
                } ,
                'node_size': {
                'var': self.node_size,
                'label': 'Node size',
                } ,
                'node_color': {
                'var': self.node_color,
                'label': 'Node color',
                } ,
                'water_node_color': {
                'var': self.water_node_color,
                'label':'Water node color' ,
                },
                'edge_color': {
                'var': self.edge_color,
                'label': 'Edge color',
                } ,
                'plot_title_fontsize': {
                'var': self.plot_title_fontsize,
                'label': 'Plot title font size',
                },
                'plot_label_fontsize': {
                'var': self.plot_label_fontsize,
                'label': 'Plot label font size',
                } ,
                'plot_tick_fontsize': {
                'var': self.plot_tick_fontsize,
                'label': 'Plot tick font size',
                },
                'plot_resolution': {
                'var': self.plot_resolution,
                'label': 'Plot resolution' ,
                } ,
                # 'figsize': {
                # 'var': self.figsize,
                # 'label': 'Figure size',
                # },
         }

    for i, (key, value) in enumerate(settings_options.items()):
        row=i
        label = value['label']
        var=value['var']
        _create_entries(self, main_frame, row, label, key)


    tk.Button(main_frame, text='Save', bg='white', fg='black', command=self._save_plot_settings, takefocus=False,  highlightbackground='white').grid(row=i+1, column=2, sticky="EW", columnspan=2)

