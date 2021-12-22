import tkinter as tk
from tkinter import ttk



def plot_settings(self):
    def _create_entries(self, frame, row, key, value):

        tk.Label(frame, text=value['label'], anchor='w',  bg='white', fg='black').grid(row=row, column=0, sticky='W')
        var = tk.Entry(frame, bg='white', fg='black', highlightbackground='white', insertbackground='black',
            validate="key", validatecommand=(self.ifnum_cmd, '%S', '%P', 0, 1000))
        var.insert(0, value['default'])
        var.grid(row=row, column=1, sticky='W')
        setattr(self, key, var)

    main_frame = ttk.LabelFrame(self.plotsettings, text='Plot settings')
    main_frame.grid(self._crate_frame_grid(0))
    main_frame.columnconfigure(0, weight=1)
    main_frame.columnconfigure(1, weight=1)

    settings_options = {
                'edge_width': {'label': 'Edge width',
                'default': '2'},
                'node_label_size': {'label': 'Node label size',
                'default': '12'},
                'edge_label_size': {'label': 'Edge label size',
                'default': '10'},
                'node_size': {'label': 'Node size',
                'default': '150'},
                'plot_title_fontsize': {'label': 'Plot title font size',
                'default': '20'},
                'plot_label_fontsize': {'label': 'Plot label font size',
                'default': '36'},
                'plot_tick_fontsize': {'label': 'Plot tick font size',
                'default': '33'},
                'plot_width': {'label': 'Plot width',
                'default': '15'},
                'plot_height': {'label': 'Plot height',
                'default': '16'},
                'plot_resolution': {'label': 'Plot resolution' ,
                'default': ''},
                # 'figsize': {
                # 'var': self.figsize,
                # 'label': 'Figure size',
                # },
         }

    for i, (key, value) in enumerate(settings_options.items()):
        row=i
        _create_entries(self, main_frame, row, key, value)


    tk.Label(main_frame, text='Graph color', anchor='w',  bg='white', fg='black').grid(row=i+1, column=0, sticky='W')
    self.graph_color = 'gray'
    graph_color_field = tk.Label(main_frame, width=2, bg=self.graph_color, anchor="w")
    graph_color_field.grid(row=i+1, column=1, sticky='W')
    graph_color_field.bind("<Button-1>", lambda x=self.graph_color, y=graph_color_field, var='graph_color':self._choose_color(x, y, var))

    tk.Label(main_frame, text='Water node color', anchor='w',  bg='white', fg='black').grid(row=i+2, column=0, sticky='W')
    self.water_node_color = '#db5c5c'
    water_node_color_field = tk.Label(main_frame, width=2, bg=self.water_node_color, anchor="w")
    water_node_color_field.grid(row=i+2, column=1, sticky='W')
    water_node_color_field.bind("<Button-1>", lambda x=self.water_node_color, y=water_node_color_field, var='water_node_color':self._choose_color(x, y, var))

    tk.Label(main_frame, text='Difference graph color', anchor='w',  bg='white', fg='black').grid(row=i+3, column=0, sticky='W')
    self.difference_graph_color = '#129fe6'
    difference_graph_color_field = tk.Label(main_frame, width=2, bg=self.difference_graph_color, anchor="w")
    difference_graph_color_field.grid(row=i+3, column=1, sticky='W')
    difference_graph_color_field.bind("<Button-1>", lambda x=self.difference_graph_color, y=difference_graph_color_field, var='difference_graph_color':self._choose_color(x, y, var))



    tk.Button(main_frame, text='Save', bg='white', fg='black', command=self._save_plot_settings, takefocus=False,  highlightbackground='white').grid(row=i+2, column=2, sticky="EW", columnspan=2)

