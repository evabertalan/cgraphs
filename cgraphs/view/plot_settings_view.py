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


                # 'graph_color': {'label': 'Node color',
                # 'default': 'gray'},
                # 'difference_graph_color': {'label': 'Edge color',
                # 'default': '#129fe6'},
                # 'water_node_color': {'label': 'Water node color' ,
                # 'default': '#db5c5c'},

    tk.Label(main_frame, text='Graph color', anchor='w',  bg='white', fg='black').grid(row=row, column=0, sticky='W')

    graph_color = 'gray'
    graph_color_field = tk.Label(main_frame, width=2, bg=graph_color, anchor="w")
    graph_color_field.grid(row=i+1, column=2, sticky='W')
    graph_color_field.bind("<Button-1>", lambda x=graph_color, y=graph_color_field, var='graph_color':self._choose_color(x, y, var))
    self.water_node_color= 'orange'
    self.difference_graph_color = 'green'


    tk.Button(main_frame, text='Save', bg='white', fg='black', command=self._save_plot_settings, takefocus=False,  highlightbackground='white').grid(row=i+2, column=2, sticky="EW", columnspan=2)

