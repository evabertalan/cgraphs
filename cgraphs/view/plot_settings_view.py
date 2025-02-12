import tkinter as tk
from tkinter import ttk


def plot_settings(self):
    def _create_entries(self, frame, row, key, value):

        tk.Label(frame, text=value["label"], anchor="w", bg="white", fg="black").grid(
            row=row, column=0, sticky="W"
        )
        var = tk.Entry(
            frame,
            bg="white",
            fg="black",
            highlightbackground="white",
            insertbackground="black",
            validate="key",
            validatecommand=(
                (self.ifnum_cmd, "%S", "%P", 0, 2000)
                if key not in ("xlabel", "ylabel")
                else None
            ),
        )
        var.insert(0, value["default"])
        var.grid(row=row, column=1, sticky="W")
        setattr(self, key, var)

    main_frame = ttk.LabelFrame(self.plotsettings, text="Plot settings")
    main_frame.grid(self._crate_frame_grid(0))
    main_frame.columnconfigure(0, weight=1)
    main_frame.columnconfigure(1, weight=1)

    settings_options = {
        "edge_width": {"label": "Edge width", "default": "2"},
        "edge_label_size": {"label": "Edge label size", "default": "10"},
        "node_size": {"label": "Node size", "default": "150"},
        "node_label_size": {"label": "Node label size", "default": "12"},
        "plot_title_fontsize": {"label": "Plot title font size", "default": "20"},
        "plot_label_fontsize": {"label": "Plot axis label font size", "default": "36"},
        "plot_tick_fontsize": {"label": "Plot tick font size", "default": "33"},
        "plot_width": {"label": "Plot width", "default": "15"},
        "plot_height": {"label": "Plot height", "default": "16"},
        "plot_resolution": {"label": "Plot resolution (dpi)", "default": "400"},
        "xlabel": {
            "label": "Plot X axis label",
            "default": "PCA projected membrane plane (Å)",
        },
        "ylabel": {"label": "Plot Y axis label", "default": "Membrane normal (Å)"},
    }

    for i, (key, value) in enumerate(settings_options.items()):
        row = i
        _create_entries(self, main_frame, row, key, value)

    tk.Label(main_frame, text="Graph color", anchor="w", bg="white", fg="black").grid(
        row=i + 1, column=0, sticky="W", pady=(2, 2), padx=(2, 2)
    )
    self.graph_color = "gray"
    graph_color_field = tk.Label(main_frame, width=2, bg=self.graph_color, anchor="w")
    graph_color_field.grid(row=i + 1, column=1, sticky="W", pady=(2, 2), padx=(2, 2))
    graph_color_field.bind(
        "<Button-1>",
        lambda x=self.graph_color, y=graph_color_field, var="graph_color": self._choose_color(
            x, y, var
        ),
    )

    tk.Label(
        main_frame, text="Water node color", anchor="w", bg="white", fg="black"
    ).grid(row=i + 2, column=0, sticky="W", pady=(2, 2), padx=(2, 2))
    self.water_node_color = "#db5c5c"
    water_node_color_field = tk.Label(
        main_frame, width=2, bg=self.water_node_color, anchor="w"
    )
    water_node_color_field.grid(
        row=i + 2, column=1, sticky="W", pady=(2, 2), padx=(2, 2)
    )
    water_node_color_field.bind(
        "<Button-1>",
        lambda x=self.water_node_color, y=water_node_color_field, var="water_node_color": self._choose_color(
            x, y, var
        ),
    )

    tk.Label(
        main_frame, text="Difference graph color", anchor="w", bg="white", fg="black"
    ).grid(row=i + 3, column=0, sticky="W", pady=(2, 2), padx=(2, 2))
    self.difference_graph_color = "#129fe6"
    difference_graph_color_field = tk.Label(
        main_frame, width=2, bg=self.difference_graph_color, anchor="w"
    )
    difference_graph_color_field.grid(
        row=i + 3, column=1, sticky="W", pady=(2, 2), padx=(2, 2)
    )
    difference_graph_color_field.bind(
        "<Button-1>",
        lambda x=self.difference_graph_color, y=difference_graph_color_field, var="difference_graph_color": self._choose_color(
            x, y, var
        ),
    )

    tk.Label(
        main_frame, text="Non-protein residue color", anchor="w", bg="white", fg="black"
    ).grid(row=i + 4, column=0, sticky="W", pady=(2, 2), padx=(2, 2))
    self.non_prot_color = "blue"
    non_prot_color_field = tk.Label(
        main_frame, width=2, bg=self.non_prot_color, anchor="w"
    )
    non_prot_color_field.grid(row=i + 4, column=1, sticky="W", pady=(2, 2), padx=(2, 2))
    non_prot_color_field.bind(
        "<Button-1>",
        lambda x=self.non_prot_color, y=non_prot_color_field, var="non_prot_color": self._choose_color(
            x, y, var
        ),
    )

    tk.Label(main_frame, text="Chain label", anchor="w", bg="white", fg="black").grid(
        row=i + 5, column=0, sticky="W"
    )
    self.show_chain_label = tk.BooleanVar()
    self.show_chain_label.set(False)
    tk.Checkbutton(
        main_frame,
        text="show",
        variable=self.show_chain_label,
        anchor="w",
        bg="white",
        fg="black",
    ).grid(row=i + 5, column=1, sticky="W")

    tk.Label(
        main_frame, text="Save figures in formats", anchor="w", bg="white", fg="black"
    ).grid(row=i + 6, column=0, sticky="W")
    self.png = tk.BooleanVar()
    self.png.set(True)
    tk.Checkbutton(
        main_frame,
        text="png",
        state="disabled",
        variable=self.png,
        anchor="w",
        bg="white",
        fg="black",
    ).grid(row=i + 6, column=1, sticky="W")

    self.eps_format = tk.BooleanVar()
    self.eps_format.set(False)
    tk.Checkbutton(
        main_frame,
        text="eps",
        variable=self.eps_format,
        anchor="w",
        bg="white",
        fg="black",
    ).grid(row=i + 7, column=1, sticky="W")

    self.svg = tk.BooleanVar()
    self.svg.set(False)
    tk.Checkbutton(
        main_frame, text="svg", variable=self.svg, anchor="w", bg="white", fg="black"
    ).grid(row=i + 8, column=1, sticky="W")

    tk.Button(
        main_frame,
        text="Save",
        bg="white",
        fg="black",
        command=self._save_plot_settings,
        takefocus=False,
        highlightbackground="white",
    ).grid(row=i + 9, column=0, sticky="EW", columnspan=2)
