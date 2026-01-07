# Installation guide
1. Download C-Graphs from the [GitHub repository](https://github.com/evabertalan/cgraphs)
2. Open a terminal window and navigate to the downloaded cgraphs folder
3. Make sure you have python3 and conda installed on your computer `conda --version`, if no version number is returned, please install conda.

**a) Linux, Mac:**

  4. `./install_cgraphs`
  5. `./start_cgraphs`
 
**b) Windows:**
  
  * optionally: create a conda environment before the following steps
  4. `conda install -c conda-forge mdanalysis`
  5. `py -m cgraphs start`

**c)** _Note:_ The installation script of C-Graphs creates a conda virtual environment. If you would like to use the tool without a virtual environment (e.g.: on Linux fonts can look different within the conda environment), please follow the following steps:

  4. `conda install -c conda-forge mdanalysis`
  5. `python3 -m cgraphs start`

# Main window
<img width="911" alt="Figure_1" src="https://github.com/user-attachments/assets/7bd18447-9224-476e-aed3-09fa677f25f2" />


# C-Graphs parameters and options
## I. Crystal structure analysis:
To calculate water clusters, individual H-bond graphs, conserved networks and conserved water sites in static crystal structures.

### Input locations:
* **Select PDB folder:** Select the folder containing the protein coordinate files (in standard Protein Data Bank (PDB) format), which are aimed to be analyzed. The tool will create a __workfolder__ in this location and save all the calculation results there.
* **Select reference file:** Select a single coordinate file as a reference. The reference structure will be used as basis to compute relative sequence identities for all other structures, for all structural overlaps, and for preparing the plots of conserved H-bond graphs. For membrane proteins, using structures oriented along the membrane normal is recommended, as this will facilitate interpretation of PCA-projected H-bond graphs.
* **Minimum sequence identity (%):** Structures from the selected PDB folder, whose sequence identity relative to the reference is below the minimum value set here will be discarded from the analyses.
* **Conservation of H-bonding groups across structures:** Percentage of structures that must include the graph nodes (H-bonding groups) to be consider the nodes as conserved and be included in the conserved H-bond graphs.
* **Superimpose PDB structures:** By default this option is selected and the tool will automatically superimpose the provided protein structures, however for more accurate results (especially by water clustering) it is recommended to pre-align the structures with a dedicated tool, provide the align structures as input, and deselect this options.

### Water cluster analysis with [DBSCAN](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html):
The default DBSCAN eps and RMS threshold parameters are recommended to use. Changing theses values influence the results of the water clustering algorithm and a parameter testing is recommended.
* **DBSCAN eps:** maximum distance parameter between two water molecules to be considered as in the neighborhood of each other. For the water clustering the default and recommended eps is 1.4 [[Bertalan, 2020]](https://www.sciencedirect.com/science/article/pii/S1047847720302070).
* **Superimposition RMS threshold:** Value to limit the deviation between the analyzed structures. DBSCAN is a distance based clustering algorithm, structures which have higher RMS deviation to the reference structures than the threshold value, are excluded from the clustering.

### Conserved network analysis with [Bridge](https://github.com/maltesie/bridge):
* **Selection:** By default, C-Graphs computes the H-bond network for protein groups and water molecules. 
Clicking on the “Selection” button opens a window in which a selection string may be used that supports the [MDAnalysis atom selection language](https://userguide.mdanalysis.org/stable/selections.html). 
  * *List of additional donors and acceptors:* To include non-protein donor and acceptor atoms in H-bond network calculations.
For example, to include a sodium ion in the H-bond graph computations, the selection string will read “protein or resname NA”, where NA is the residue name of the sodium ion in the PDB file; additionally, the name of sodium ion, e.g., Na, needs to be added to the list of acceptor and donor atoms. 

* **H-bond criteria:** All H-bond computations are performed according to the default geometric criteria of H-bonding.
  * *distance*: The distance between the heavy atoms of the H-bond. When the input protein structures are experimental structures that lack H atoms, H-bonds are computed with a single distance-based criterion. The default value is 3.5Å.
  * *angle*: Threshold value for the angle formed by the acceptor heavy atom, the H atom, and the donor heavy atom. When protein structures are read from an MD simulation trajectory, this additional H-bond angle criterion can be turned on. The default value is 60°.

* **Color nodes by:** An option to color nodes of the H-bond network by assigned values.
  * *Propka file*: Color titratable residues by their pKa values calculated with [PROPKA3](https://github.com/jensengroup/propka). Please the .pka files next to the PDB files in the PDB folder with the following naming convention: if the PDB file is named 6i9k.pdb place next to it the 6i9k.pka to color the nodes by their pKa value.
  * *B-factor*: Color the nodes of the H-bond network (based on the *Residues to color* selection, which defaults to protein) using the average B-factors of the atoms in each residue.
  * *User defined values* Nodes of the network can be colored by any value provided in an external data file. The file has to be put in the same folder where the input PDB files are located and follow the naming convention: if the PDB file is named 6i9k.pdb the file with the external values has to be named 6i9k_data.txt
  The .txt file has to have the following content RES_NAME  RES_ID  SEG_ID  VALUE, e.g:

                  ASN       1       A       7.9  
                  ASP       6       A       2.5  
  * **Residues to color...** Define which nodes to color in the H-bond by using the [MDAnalysis atom selection language](https://userguide.mdanalysis.org/stable/selections.html). By default only protein nodes are colored. In order to include waters add e.g: *protein or resname HOH*. By default the virdis color map of matplotlib is used, here there is an option to select other supported color maps.


* **Plot for each structure:** Perform the selected analysis types and create the plots for each of the structures included in the analysis from the PDB folder.
  * *Individual network:* Calculate the H-bond or water wire network with Bridge and plot the graph.
  * *Difference graph:* Computes the difference between the H-bond network of a given structure and the conserved graph. Graph nodes and edges present only in the given structure of interest for analysis are colored blue, conserved amino acid residues and their connections, gray.
  * *Linear length:* Creates a two-dimensional plot with the linear length of H-bond networks as vertical axis, and the linear projection of the 3rd, Cartesian Z coordinate of the amino acid residue node, as horizontal axis.
  * *H-bond distances:* Calculates the H-bond network of the individual PDB structures and writes the length in Å of the H-bond of the graph edges on the plot. Additionally creates a file (PDBID_H-bond_graph_distances.txt) with the length of each edges.


#### H-bond network:
* **Include crystallographic waters:** When the PBD files don't contain crystallographic water molecules relevant to the analysis, they can be excluded form the computations of the H-bond graphs.
* **Use water cluster coordinates:** Include the cluster centers obtained with DBSCAN as conserved water sites in the conserved H-bond graph calculation. These water clusters will appear with W1, W2, …, Wn labels on the plot,  where n is the number of water clusters.
* **Include sidechain-backbone interactions:** Include H-bonds between the sidechain of protein groups and the backbone of the protein in the H-bond network calculation.

#### Water wire network:
* **Maximum number of waters in bridge:** The maximum number of water molecules allowed in the connection between two protein sidechains.

## II. MD trajectory analysis:
To calculate and compare networks of one or more set of simulations of proteins of the same family. And calculate conserved networks when multiple set of simulations are selected.
* **Location of workfolder:** The location where the workfolder is going to be created. All results are saved in the created workfolder.

* **Select simulation:** C-Graphs can compute the water wire graph object for one simulation at a time. When the calculation is completed, a graph for an other simulation can be constructed or the already constructed graphs can be selected for comparisons.
  * *Select PSF:* Select the protein structure file.
  * *Select DCDs:* Select one or multiple trajectory files of the same simulations. Prior to this analysis, simulation atoms must be wrapped into a single unitcell (PBC wrap).
  * *Name as:* Each set of simulations must be given a unique name. The results of the calculations are going to be saved in a subfolder with this name.

* **Select graphs to compare:** Select graph objects which were constructed in the previous steps or at an other time. These graph objects are located in the *workfolder/graph_objects* folder.

* **Minimum H-bond occupancy:** The minimum occupancy of H-bonds, the percentage of the time of the trajectory segment used for analyses, during which the H-bond criteria are met. Connections with lower occupancy values are filtered out from the graph. 


See details for **Conservation of H-bonding groups across structures** and **Plot for each structure** options above.

## III. Compare 2 structures:
To perform direct comparison between the H-bond networks of two structures. The comparison can be performed between two PDB structures or between two separate set of MD simulations of proteins of the same family.
As an output it will create the conserved graph of the two structure and a comparison graph, where the common nodes and edges are colored with gray and the unique nodes and edges are colored by the selected color of the corresponding structure.

For detailed description of the options and parameter on this tab see the sections above.

* **Color common nodes by:** Input files have to follow the same naming convention as mentioned above in the *Conserved network analysis* section. On the resulted conserved graph plots, the common nodes of the two structures will be colored by the differences of the given value. E.g in case of Propka file, each common node will be colored by the pKa difference in the two structure.
* **Calculate H-bond distance difference:** If this option is selected an additional plot will be created, where the common edges (H-bonds) of the two networks will be colored by the length difference of the given H-bond in the two structure.


# How to cite:
* When using the tool please cite:

**C-Graphs Tool with Graphical User Interface to Dissect Conserved Hydrogen-Bond Networks: Applications to Visual Rhodopsins**
Éva Bertalan, Elena Lesca, Gebhard F. X. Schertler, and Ana-Nicoleta Bondar
Journal of Chemical Information and Modeling 2021 61 (11), 5692-5707
[DOI: 10.1021/acs.jcim.1c00827](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.1c00827)

* For version 2.0 cite the following paper as well:
  
**DNET: A Graph-Based Tool and Workflow for Dynamic Hydrogen-Bond Networks and Applications for Visual Rhodopsins**
Éva Bertalan, Matthew J. Rodrigues, Deborah Walter, Gebhard F. X. Schertler, and Ana-Nicoleta Bondar
Journal of Chemical Theory and Computation Article ASAP
[DOI: 10.1021/acs.jctc.5c01366](https://pubs.acs.org/doi/10.1021/acs.jctc.5c01366)


