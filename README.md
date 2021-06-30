# Installation guide
1. Download C-Graph from the [GitHub repository](https://github.com/evabertalan/cgraphs)
2. Open a terminal window and go to the cgraphs folder `cd <path_to_the_folder>/cgraphs`
3. Make sure you have python3 and conda installed on your computer `conda --version`, if no version number is returned, please install conda.

**a) Linux, Mac:**

  4. `./install_cgraphs`
  5. `./start_cgraphs`
 
**b) Windows:**
  
  * optionally: create a conda environment before the following steps
  4. `conda install -c conda-forge mdanalysis`
  5. `py -m cgraphs start`

**c)** _Note:_ The installation script of C-Graph creates a conda virtual environment. If you would like to use the tool without a virtual environment (e.g.: on Linux fonts can look different within the conda environment), please follow the following steps:

  4. `conda install -c conda-forge mdanalysis`
  5. `python3 -m cgraphs start`

# Main window
<img width="475" alt="Screenshot 2021-06-27 at 9 23 34" src="https://user-images.githubusercontent.com/15729207/123536226-6569de00-d729-11eb-9386-9535a8bd9bd1.png">


# C-Graph parameters and options
## I. Crystal structure analysis:
To calculate water clusters, individual H-bond graphs, conserved networks and conserved water sites in static crystal structures.

### Input locations:
* **Select PDB folder:** Select the folder containing the protein coordinate files (in standard Protein Data Bank (PDB) format), which are aimed to be analyzed.
* **Select reference file:** Select a single coordinate file as a reference. The reference structure will be used as basis to compute relative sequence identities for all other structures, for all structural overlaps, and for preparing the plots of conserved H-bond graphs. For membrane proteins, using structures oriented along the membrane normal is recommended, as this will facilitate interpretation of PCA-projected H-bond graphs.
* **Minimum sequence identity (%):** Structures from the selected PDB folder, whose sequence identity relative to the reference is below the minimum value set here will be discarded from the analyses.
Results are saved in the __workfolder__ created in the selected PDB folder.

### Water cluster analysis with DBSCAN:
The default DBSCAN eps and RMS threshold parameters are recommended to use. Changing theses values influence the results of the water clustering algorithm and a parameter testing is recommended.
* **DBSCAN eps:** maximum distance parameter between two water molecules to be considered as in the neighbourhood of each other. For the water clustering the default and recommended eps is 1.4 [[Bertalan, 2020]](https://www.sciencedirect.com/science/article/pii/S1047847720302070).
* **Superimposition RMS threshold:** Value to limit the deviation between the analyzed structures. DBSCAN is a distance based clustering algorithm, structures which have higher RMS deviation to the reference structures than the threshold value, are excluded from the clustering.

### Conserved network analysis with Bridge:
* **H-bond criteria:** All H-bond computations are performed according to the default geometric criteria of H-bonding.
  * *distance*: The distance between the heavy atoms of the H-bond. When the input protein structures are experimental structures that lack H atoms, H-bonds are computed with a single distance-based criterion. The default value is 3.5Å.
  * *angle*: Threshold value for the angle formed by the acceptor heavy atom, the H atom, and the donor heavy atom. When protein structures are read from an MD simulation trajectory, this additional H-bond angle criterion can be turned on. The default value is 60°.

* **Conservation of H-bonding groups across structures:** Percentage of structures that must include the graph nodes (H-bonding groups) to be consider the nodes as conserved and be included in the conserved H-bond graphs.
* **Plot for each structure:** Perform the selected analysis types and create the plots for each of the structures included in the analysis from the PDB folder.
  * *Individual network:* Calculate the H-bond or water wire network with Bridge and plot the graph.
  * *Difference graph:* Computes the difference between the H-bond network of a given structure and the conserved graph. Graph nodes and edges present only in the given structure of interest for analysis are colored blue, conserved amino acid residues and their connections, gray.
  * *Linear length:* Creates a two-dimensional plot with the linear length of H-bond networks as vertical axis, and the linear projection of the 3rd, Cartesian Z coordinate of the amino acid residue node, as horizontal axis.


#### H-bond network:
* **Include crystallographic waters:** When the PBD files don't contain crystallographic water molecules relevant to the analysis, they can be excluded form the computations of the H-bond graphs.
* **Use water cluster coordinates:** Include the cluster centers obtained with DBSCAN as conserved water sites in the conserved H-bond graph calculation. These water clusters will appear with W1, W2, …, Wn labels on the plot,  where n is the number of water clusters.
* **Include sidechain-backbone interactions:** Include H-bonds between the sidechain of protein groups and the backbone of the protein in the H-bond network calculation.

#### Water wire network:
* **Maximum number of waters in bridge:** The maximum number of water molecules allowed in the connection between two protein sidechains.

## II. MD trajectory analysis:
To calculate and compare networks of one or more set of simulations of proteins of the same family. And calculate conserved networks when multiple set of simulations are selected.
* **Location of workfolder:** The location where the workfolder is going to be created. All results are saved in the created workfolder.

* **Select simulation:** C-Graph can compute the water wire graph object for one simulation at a time. When the calculation is completed, a graph for an other simulation can be constructed or the already constructed graphs can be selected for comparisons.
  * *Select PSF:* Select the protein structure file.
  * *Select DCDs:* Select one or multiple trajectory files of the same simulations. Prior to this analysis, simulation atoms must be wrapped into a single unitcell (PBC wrap).
  * *Name as:* Each set of simulations must be given a unique name. The results of the calculations are going to be saved in a subfolder with this name.

* **Select graphs to compare:** Select graph objects which were constructed in the previous steps or at an other time. These graph objects are located in the *workfolder/graph_objects* folder.

* **Minimum H-bond occupancy:** The minimum occupancy of H-bonds, the percentage of the time of the trajectory segment used for analyses, during which the H-bond criteria are met. Connections with lower occupancy values are filtered out from the graph. 


See details for **Conservation of H-bonding groups across structures** and **Plot for each structure** options above.

## III. Compare 2 structures:
To perform direct comparison between the H-bond networks of two structures. The comparison can be performed between two PDB structures or between two separate set of MD simulations of proteins of the same family.

For detailed description of the options and parameter on this tab see the sections above.

# How to cite:
1.	Bertalan E, Lesnik S, Bren U, Bondar A-N. [Protein-water hydrogen-bond networks of G Protein-Coupled Receptors: Graph-based analyses of static structures and molecular dynamics.](https://www.sciencedirect.com/science/article/pii/S1047847720302070) Journal of Structural Biology 212, 107634 (2020)
2.	Bertalan E, Lesca E, Deupi X, Schertler GFX, Bondar A-N. C-Graphs Analyzer of Conserved Hydrogen-Bond Graphs: Applications to Hydrogen-Bond Networks of Visual Rhodopsins (2021, submitted)
