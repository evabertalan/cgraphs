## Installation guide

## C-Graph parameters and options
### I. Crystal structure analysis:

#### Input locations:
* **Select PDB folder:** Select the folder containing the protein coordinate files (in standard Protein Data Bank (PDB) format), which are aimed to be analized
* **Select reference file:** Select a single coordinate file as a reference. The reference structure will be used as basis to compute relative sequence identities for all other structures, for all structural overlaps, and for preparing the plots of conserved H-bond graphs. For membrane proteins, using structures oriented along the membrane normal is recomennded, as this will facilitate interpretation of PCA-projected H-bond graphs.
* **Minimum sequence identity (%):** Structures from the selected PDB folder, whose sequence identity relative to the reference is below the minimum value set here will be discarded from the analyses.
Results are saved in the __workfolder__ created in the selected PDB folder.

#### Water cluster analysis with DBSCAN:
Default DBSCAN eps and RMS thresholds parameters are recomenned to use, as changing theses values strongly influence the results of the water clustering algorithm.
* **DBSCAN eps:** maximum distance paramter between two water molecules to be considered as in the neighborhood of eahc other. For the water clustering the default and recommended eps is 1.4 [[Bertalan, 2020]](https://www.sciencedirect.com/science/article/pii/S1047847720302070).
* **Superimposiont RMS threshold:** Value to limit the deviation between the analyzed structures. DBSCAN is a distance based clustering algorithm, structres which have higher RMS deviation to the reference structres than the threshold value, are excluded from the clutering.

#### Conserved network analysisi with Bridge:
* **H-bond criteria:** All H-bond computations are performed according to the default geometric criteria of H-bonding.
  * __distance__: The distcance between the heavy atoms of the H-bond. When the input protein structures are experimental structures that lack H atoms, H-bonds are computed with a single distance-based criterion. The defualt value is 3.5Å.
  * __angle__: Threshold value for the angle formed by the acceptor heavy atom, the H atom, and the donor heavy atom. When protein structures are read from an MD simulation trajectory, this additional H-bond angle criterion can be turned on. The defualt value is 60°.

* **Conservation of H-bonfing groups across structures:** Percentage of structures that must include the graph nodes (H-bonding groups) to be consider the nodes as conserved and be included in the conserved H-bond graphs.
* **Plot for each structure:** Perform the selected analysis types and create the plots for each of the stuctured included in the analysis from the PDB folder.

##### H-bond network:
* **Include crystallographic waters:** When the PBD files don't contain crystallographic the water molecules from the PDB file can be excluded form the computions of the H-bond graphs.
* **Use water cluster coordinats:** To include the cluster centers obtained with DBSCAN as conserved water sites in the conserved H-bond graph calcualtion. These water clusters will appear with W1, W2, …, Wn labels on the plot,  where n is the number of water clusters.
* **Include sidechain-backbone interactions:** Include H-bonds between the sideshain of protein gtoups and the  backbone of the proetin in the H-bond network calculation.

##### Water wire network:
* **Maximum number of waters in bridge:** The maximum number of water molecules allowed in the connection between two protein sidechains.

### II. MD trajectory analysis:
* **Location of workfolder:**

### III. Compare 2 structures:

