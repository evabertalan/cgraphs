import os
import shutil
import logging
import json
import pickle
import warnings
import numpy as np
import networkx as nx
import MDAnalysis as _mda

# from Bio.SVDSuperimposer import SVDSuperimposer
# from Bio import pairwise2
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl


warnings.filterwarnings("ignore")

amino_d = {
    "CYS": "C",
    "ASP": "D",
    "SER": "S",
    "GLN": "Q",
    "LYS": "K",
    "ILE": "I",
    "PRO": "P",
    "THR": "T",
    "PHE": "F",
    "ASN": "N",
    "GLY": "G",
    "HIS": "H",
    "LEU": "L",
    "ARG": "R",
    "TRP": "W",
    "ALA": "A",
    "VAL": "V",
    "GLU": "E",
    "TYR": "Y",
    "MET": "M",
    "HSD": "H",
    "HSE": "H",
    "BWX": "X",
}
water_def = "(resname TIP3 and name OH2) or (resname HOH and name O) or (resname TIP4 and name OH2)"
water_types = ["HOH", "TIP3", "TIP4"]


def get_plot_parameters(plot_parameters):
    default_plot_parameters = {
        "edge_width": (
            plot_parameters["edge_width"]
            if "edge_width" in plot_parameters.keys()
            else 2
        ),
        "node_label_size": (
            plot_parameters["node_label_size"]
            if "node_label_size" in plot_parameters.keys()
            else 12
        ),
        "edge_label_size": (
            plot_parameters["edge_label_size"]
            if "edge_label_size" in plot_parameters.keys()
            else 10
        ),
        "node_size": (
            plot_parameters["node_size"]
            if "node_size" in plot_parameters.keys()
            else 150
        ),
        "graph_color": (
            plot_parameters["graph_color"]
            if "graph_color" in plot_parameters.keys()
            else "gray"
        ),
        "water_node_color": (
            plot_parameters["water_node_color"]
            if "water_node_color" in plot_parameters.keys()
            else "#db5c5c"
        ),
        "difference_graph_color": (
            plot_parameters["difference_graph_color"]
            if "difference_graph_color" in plot_parameters.keys()
            else "#129fe6"
        ),
        "non_prot_color": (
            plot_parameters["non_prot_color"]
            if "non_prot_color" in plot_parameters.keys()
            else "blue"
        ),
        "plot_title_fontsize": (
            plot_parameters["plot_title_fontsize"]
            if "plot_title_fontsize" in plot_parameters.keys()
            else 20
        ),
        "plot_label_fontsize": (
            plot_parameters["plot_label_fontsize"]
            if "plot_label_fontsize" in plot_parameters.keys()
            else 36
        ),
        "plot_tick_fontsize": (
            plot_parameters["plot_tick_fontsize"]
            if "plot_tick_fontsize" in plot_parameters.keys()
            else 33
        ),
        "plot_resolution": (
            plot_parameters["plot_resolution"]
            if "plot_resolution" in plot_parameters.keys()
            else 400
        ),
        "figsize": (
            plot_parameters["figsize"]
            if "figsize" in plot_parameters.keys()
            else (15, 16)
        ),
        "formats": (
            plot_parameters["formats"]
            if "formats" in plot_parameters.keys()
            else ["png"]
        ),
        "show_chain_label": (
            plot_parameters["show_chain_label"]
            if "show_chain_label" in plot_parameters.keys()
            else False
        ),
    }
    return default_plot_parameters


def create_logger(folder):
    logger = logging.getLogger("cgraphs")
    if not len(logger.handlers):
        logger.setLevel(logging.INFO)
        fh = logging.FileHandler(folder + "cgraphs_logs.log")
        with open(folder + "cgraphs_logs.log", "w"):
            pass
        fh.setLevel(logging.DEBUG)
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        fh_form = logging.Formatter(
            "%(asctime)s - %(levelname)s: %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
        )
        ch_form = logging.Formatter("%(levelname)s: %(message)s")
        fh.setFormatter(fh_form)
        ch.setFormatter(ch_form)
        logger.addHandler(fh)
        logger.addHandler(ch)
    return logger


def create_directory(directory):
    if not os.path.isdir(directory):
        os.makedirs(directory)
    return directory


def delete_directory(directory):
    if os.path.isdir(directory):
        shutil.rmtree(directory)


def get_pdb_files(folder):
    return [file for file in os.listdir(folder) if file.endswith(".pdb")]


def get_files(folder, endswith):
    return [file for file in os.listdir(folder) if file.endswith(endswith)]


def pickle_write_file(path, obj):
    with open(path, "wb") as fp:
        pickle.dump(obj, fp)


def pickle_load_file(path):
    with open(path, "rb") as fp:
        obj = pickle.load(fp)
    return obj


def json_write_file(path, obj):
    with open(path, "w", encoding="utf-8") as fp:
        json.dump(obj, fp)


def get_node_name(node):
    return node


def get_node_name_pats(node, with_group=False):
    if with_group:
        return (
            node.split("-")[0],
            node.split("-")[1],
            str(int(node.split("-")[2])),
            node.split("-")[3],
        )
    else:
        return node.split("-")[0], node.split("-")[1], str(int(node.split("-")[2]))


def concatenate_arrays(arrays):
    concatenated = []
    for arr in arrays:
        if arr.ndim > 1:
            for row in arr:
                concatenated.append(row)
        elif arr.size != 0:
            concatenated.append(arr)
    return np.array(concatenated)


def load_pdb_structure(pdb_file):
    u = _mda.Universe(pdb_file)
    return u


def water_in_pdb(pdb_file):
    structure = load_pdb_structure(pdb_file)
    waters = structure.select_atoms(water_def)
    return waters


def water_coordinates(pdb_file):
    waters = water_in_pdb(pdb_file)
    water_coord = waters.positions
    return np.array(water_coord)


def get_color_map(color_info, color_map="viridis", center=None, normalize=None):
    cmap = cm.get_cmap(color_map, len(color_info))
    _vals = np.array(list(color_info.values()), dtype=float)
    if len(color_info):
        if len(color_info) > 1:
            scaled_values = (_vals - _vals.min()) / (_vals.max() - _vals.min())
            value_colors = {
                key: cmap(scaled_values[i])
                for i, (key, values) in enumerate(color_info.items())
            }
        else:
            value_colors = {
                key: cmap(_vals[i])
                for i, (key, values) in enumerate(color_info.items())
            }
        if center:
            max_val = np.max([np.abs(_vals.min()), np.abs(_vals.max())])
            norm = mpl.colors.Normalize(vmin=-1 * max_val, vmax=max_val)
        elif normalize:
            norm = norm = mpl.colors.Normalize(vmin=normalize[0], vmax=normalize[1])
        else:
            norm = mpl.colors.Normalize(vmin=_vals.min(), vmax=_vals.max())
        return value_colors, cmap, norm
    else:
        return {}, cmap, None


def get_sequence(pdb_file, selection="protein and name CA"):
    structure = load_pdb_structure(pdb_file)
    protein = structure.select_atoms(selection)
    seq = ""
    for res in protein:
        seq += amino_d[res.resname]
    return seq


def get_best_alignment(alignments):
    best_score = 0
    best_i = 0
    for i, alignment in enumerate(alignments):
        if best_score <= alignment.score:
            best_score = alignment.score
            best_i = i

    best_alginment = alignments[best_i]
    return best_alginment, best_i


## INACTIVE option, but keep code in file, not jues in GitHub
# def get_residue_conservations(pdb_file, conservation_info):
#     conservation_info = np.loadtxt(conservation_info, skiprows=1, dtype=str)
#     conservation_sequence = ''.join(conservation_info[:,1])
#     structure = load_pdb_structure(pdb_file)
#     selection='protein and name CA and segid A'
#     protein = structure.select_atoms(selection)
#     seq = get_sequence(pdb_file, selection)

#     alignments =  pairwise2.align.globalxx(seq, conservation_sequence)
#     best_alginment, best_i = get_best_alignment(alignments)
#     conservation = {}

#     cons = []
#     for i, res in enumerate(alignments[best_i].seqB):
#         if res != '-':
#             index = np.where(conservation_info == res)[0][0]
#             cons.append(( res , conservation_info[index][2] ))
#             conservation_info = conservation_info[index+1:]
#         else: cons.append((res, 0))

#     res_index = 0
#     for i, res in enumerate(alignments[best_i].seqA):
#         if res != '-':
#             conservation.update({f'{protein[res_index].segid}-{protein[res_index].resname}-{protein[res_index].resid}': cons[i]})
#             res_index += 1
#     return conservation


def align_sequence(logger, pdb_ref, pdb_move, threshold=0.75):
    ref_sequence = get_sequence(pdb_ref, selection="protein and name CA")
    move_sequence = get_sequence(pdb_move, selection="protein and name CA")
    alignments = pairwise2.align.globalxx(ref_sequence, move_sequence)
    pdb_name = pdb_move.split("/")[-1]

    best_alginment, best_i = get_best_alignment(alignments)

    if len(best_alginment.seqA) != len(best_alginment.seqB):
        logger.warning("Aligned sequences have different lenght")
        logger.info("Thus " + pdb_name + " is excluded from further analysis.")
        return None, None
    if best_alginment.score / len(alignments[best_i].seqA) <= threshold:
        logger.warning(
            "Sequences of "
            + pdb_name
            + " has lower sequence identity than the threshold value ("
            + str(threshold * 100)
            + "%) compared to the reference structure"
        )
        logger.info("Thus " + pdb_name + " is excluded from further analysis.")
        return None, None
    if best_alginment.score / len(alignments[best_i].seqB) <= threshold:
        logger.warning(
            "Sequences of "
            + pdb_name
            + " has lower sequence identity than the threshold value ("
            + str(threshold * 100)
            + "%) compared to the reference structure"
        )
        logger.info("Thus " + pdb_name + " s excluded from further analysis.")
        return None, None
    return best_alginment.seqA, best_alginment.seqB


def superimpose_aligned_atoms(
    logger,
    seq_ref,
    pdb_ref,
    seq_move,
    pdb_move,
    save_file_to="",
    superimposition_threshold=5,
):
    if save_file_to == "":
        save_file_to = retrieve_pdb_code(pdb_move, ".pdb")
    else:
        save_file_to = save_file_to.split(".pdb")[0]
    # TODO: maybe creae regex or parameter to filnave OR retihnik this filename conscept
    pdb_name = pdb_move.split("/")[-1]

    ref_struct = load_pdb_structure(pdb_ref)
    move_struct = load_pdb_structure(pdb_move)
    move_struct.atoms.write(save_file_to + "_superimposed.pdb")
    return move_struct
    # ref_atoms = ref_struct.select_atoms('protein and name CA')
    # move_atoms = move_struct.select_atoms('protein and name CA')

    # unique_seg_move, unique_seg_ref = np.unique(move_atoms.segids), np.unique(ref_atoms.segids)
    # if len(unique_seg_move)==1 and len(unique_seg_ref)==1 and unique_seg_move != unique_seg_ref:
    #     logger.warning('Chains must have the same ID. To compare the conserved graph of multiple structures, same segments has to have the same chain ID.')
    #     logger.info('Chins of '+pdb_name+' has different chain ID than the reference structure. Thus excluded from further analysis.')
    #     return None

    # ref_atoms_pos = []
    # move_atoms_pos = []
    # i = -1
    # j = -1
    # for r, m in zip(seq_ref, seq_move):
    #     if r  != '-': i = i+1
    #     if m  != '-': j = j+1
    #     if (r  != '-' and m  != '-') and (r == m) and (r in amino_d.values()) and (ref_atoms[i].segid == move_atoms[j].segid):
    #         ref_atoms_pos.append(ref_atoms[i].position)
    #         move_atoms_pos.append(move_atoms[j].position)

    # move_atoms_pos = np.array(move_atoms_pos,  dtype='float64').reshape(-1, 3)
    # ref_atoms_pos = np.array(ref_atoms_pos,  dtype='float64').reshape(-1, 3)

    # sup = SVDSuperimposer()
    # sup.set(ref_atoms_pos, move_atoms_pos)
    # sup.run()
    # rot, tran = sup.get_rotran()
    # if sup.get_rms() > superimposition_threshold:
    #     logger.warning('Automatic superimposition of '+pdb_name+' was not sucessful. RMS '+str(round(sup.get_rms(),3))+' is too high. Please provide a PDB file superimposed to the reference structure. This structure is excluded from further analysis.')
    #     return
    # else:
    #     rot = rot.astype('f')
    #     tran = tran.astype('f')
    #     move_struct.atoms.positions = np.dot(move_struct.atoms.positions, rot) + tran
    #     move_struct.atoms.write(save_file_to+'_superimposed.pdb')

    #     logger.info('Superimposition RMS value of '+pdb_name+' to the reference structure is: '+str(round(sup.get_rms(),3)))
    #     logger.debug('Superimposed file is saved as: '+str(save_file_to+'_superimposed.pdb'))
    #     return move_struct


def get_connected_components(graph):
    return list(nx.connected_components(graph))


def get_water_coordinates(protein_chain, res_index):
    # FIX water id issue from mdhbond --> issue from MDAnalysis
    # if int(res_index) > 10000: res_index = int(res_index) - 10000
    # return protein_chain[('W', int(res_index), ' ')]['O'].get_coord()
    sel = protein_chain.select_atoms(water_def + " and resid " + res_index)
    if len(sel.positions):
        return sel.positions[0]
    else:
        print(
            f"Water {res_index} not found in the PDB file. INFO: https://github.com/evabertalan/cgraphs/blob/main/README.md."
        )
        return None


def calculate_connected_compontents_coordinates(
    connected_components, struct_object, option="pdb"
):
    all_chains = []
    for connected_chain in connected_components:
        chain_details = []
        for node in list(connected_chain):
            chain_id, res_name, res_id = (
                node.split("-")[0],
                node.split("-")[1],
                node.split("-")[2],
            )
            if option == "pdb":
                chain = struct_object.select_atoms("segid " + chain_id)
                if res_name in water_types:
                    coords = get_water_coordinates(chain, res_id)
                elif res_name in amino_d.keys():
                    coords = chain.select_atoms(
                        f"resname BWX or protein and name CA and resid {res_id}"
                    ).positions[0]
                else:
                    coords = chain.select_atoms(
                        f"resname {res_name} and resid {res_id}"
                    ).positions[0]

            else:
                coords = struct_object.select_atoms("resid " + res_id).positions[0]
            if coords is not None:
                chain_details.append(np.array([node, coords], dtype="object"))
        all_chains.append(chain_details)

    return [c for c in sorted(all_chains, key=len, reverse=True)]


def calculate_pca_positions(coordinates):
    pca_positions = {}
    XY = [i[0:2] for i in coordinates.values()]
    pca = PCA(n_components=1)
    xy = pca.fit_transform(XY)

    for i, (key, value) in enumerate(coordinates.items()):
        pca_positions[key] = [xy[i][0], value[2]]
    return pca_positions


def check_projection_sign(projection, reference):
    for i in range(len(reference.keys())):
        _k = list(reference.keys())[i]
        if _k in projection.keys() and reference[_k][0] > 5:
            if np.sign(reference[_k][0]) != np.sign(projection[_k][0]):
                t_projection = {}
                for k, v in projection.items():
                    t_projection.update({k: [v[0] * -1, v[1]]})
                projection = t_projection
                return projection
    return projection


def is_conserved_edge(other_graph_edges, e0, e1, with_group=False):
    conserved_edge = (
        len(np.where((other_graph_edges == [e0, e1]).all(axis=1))[0]) != 0
        or len(np.where((other_graph_edges == [e1, e0]).all(axis=1))[0]) != 0
    )
    conserved_edge_with_water = False
    for edge in other_graph_edges:
        if (
            e0 in edge
            and e1.split("-")[1] in water_types
            and (edge[0].startswith("X-w") or edge[1].startswith("X-w"))
        ) or (
            e1 in edge
            and e0.split("-")[1] in water_types
            and (edge[0].startswith("X-w") or edge[1].startswith("X-w"))
        ):
            conserved_edge_with_water = True
    return conserved_edge or conserved_edge_with_water


def retrieve_pdb_code(file_path, split_by):
    """split_by e.g.: '.pdb'"""
    return file_path.split("/")[-1].split(split_by)[0]


def _average_timeseries(hbond_dict):
    return {key: np.mean(hbond_dict[key]) for key in hbond_dict}


def get_edge_params(wba, edges):
    average_water_per_wire = wba.compute_average_water_per_wire()
    occupancy_per_wire = _average_timeseries(wba.filtered_results)

    keys = []
    waters = []
    occ_per_wire = []
    for edge in edges:
        key = str(edge[0]) + ":" + str(edge[1])
        keys.append(key)

        if key in average_water_per_wire:
            waters.append(average_water_per_wire[key])
        else:
            key = str(edge[1]) + ":" + str(edge[0])
            waters.append(average_water_per_wire[key])

        if key in occupancy_per_wire:
            occ_per_wire.append(occupancy_per_wire[key])
        else:
            key = str(edge[1]) + ":" + str(edge[0])
            occ_per_wire.append(occupancy_per_wire[key])

    return waters, occ_per_wire, keys


def edge_info(wba, edges):
    waters, occ_per_wire, keys = get_edge_params(wba, edges)
    edge_info = {}

    for w, o, k in zip(waters, occ_per_wire, keys):
        edge_info.update({k: {"waters": np.round(w, 1), "occupancy": o}})
    return edge_info


def write_text_file(file_path, text_content, logger=None):

    if file_path.endswith(".txt") and isinstance(text_content, list):
        f = open(file_path, "w")
        f.writelines(text_content)
        f.close()
    elif logger:
        logger.error(
            "The file name has to end to .txt and the text content has to be passed as a list."
        )


# TODO set back plot size from git
def create_plot(title="", xlabel="", ylabel="", plot_parameters={}):
    fig, ax = plt.subplots(figsize=plot_parameters["figsize"])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_title(title, fontsize=plot_parameters["plot_title_fontsize"])
    ax.set_xlabel(xlabel, fontsize=plot_parameters["plot_label_fontsize"])
    ax.set_ylabel(ylabel, fontsize=plot_parameters["plot_label_fontsize"])
    ax.tick_params(axis="x", labelsize=plot_parameters["plot_tick_fontsize"])
    ax.tick_params(axis="y", labelsize=plot_parameters["plot_tick_fontsize"])
    return fig, ax


def read_propka_file(file_path, selected_nodes):
    # this whole function should be refactored
    propka_info = {}
    with open(file_path, "r") as f:
        lines = f.readlines()
    is_summary = False
    for line in lines:
        # TODO: handle values higher than 1000 --> now they are not paresed
        line = line.lstrip()
        parts = str(f"{line[0:3]} {line[3:]}").split()
        if parts and parts[0] == "SUM":
            is_summary = True
        if parts and is_summary and parts[0] in amino_d.keys():
            res_name = parts[0]
            try:
                int(parts[1])
                res_id = parts[1]
                chain = parts[2]
                pka = parts[3]
            except:
                # rahter split on number and string
                res_id = parts[1][:-1]
                chain = parts[1][-1]
                pka = parts[2]
            selection = selected_nodes.select_atoms(
                f"resname {res_name} and resid {res_id} and segid {chain}"
            )
            if len(selection) and abs(float(pka)) < 50:
                propka_info.update({f"{chain}-{res_name}-{res_id}": pka})
        elif parts and is_summary and parts[0] in selected_nodes.resnames:
            res_name = parts[0]
            chain = parts[2]
            pka = parts[3]
            selection = selected_nodes.select_atoms(
                f"resname {res_name} and segid {chain}"
            )
            res_id = selection.resids[0]
            if len(selection) and abs(float(pka)) < 50:
                propka_info.update({f"{chain}-{res_name}-{res_id}": pka})
    return propka_info


def read_color_data_file(pdb_id, pdb_root_folder, selected_nodes):
    file_endings = ["_data.txt", "_color.txt", "data.txt", "color.txt"]
    # TODO add csv support

    for ending in file_endings:
        if os.path.isfile(f"{pdb_root_folder}/{pdb_id}{ending}"):
            color_file = f"{pdb_root_folder}/{pdb_id}{ending}"
            break
        else:
            color_file = ""

    if color_file.endswith("txt"):
        content = np.loadtxt(color_file, dtype=str)
    else:
        return {}

    color_info = {}
    for line in content:
        try:
            res_name = (
                list(amino_d.keys())[list(amino_d.values()).index(line[0])]
                if len(line[0]) == 1
                else line[0]
            )
            res_id, seg_id, value = line[1], line[2], line[3]
            selection = selected_nodes.select_atoms(
                f"resname {res_name} and resid {res_id} and segid {seg_id}"
            )
            if len(selection):
                color_info.update({f"{seg_id}-{res_name}-{res_id}": value})
        except:
            return {}

    return color_info


def read_edge_color_data(name, psf_file):
    folder = os.path.join(*[os.sep] + psf_file.split(os.sep)[:-1])
    if os.path.exists(folder):
        edge_color_file = [
            file for file in os.listdir(folder) if file.endswith("color_edges.txt")
        ]
        if len(edge_color_file):
            edge_colorc_data = np.loadtxt(
                os.path.join(folder, edge_color_file[0]), dtype=str
            )

            edge_value_dict = {
                tuple((edge[0], edge[1])): edge[2] for edge in edge_colorc_data
            }
            return edge_value_dict
    return {}
