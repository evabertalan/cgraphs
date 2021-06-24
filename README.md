## Installation guide

## C-Graph parameters and options
### I. Crystal structure analysis:

#### Input locations:
* **Select PDB folder:** Select the folder containing protein coordinate files (in standard Protein Data Bank (PDB) format), which are aimed to be analized
* **Select reference file:** Select a single coordinate file as a reference. The reference structure will be used as basis to compute relative sequence identities for all other structures, for all structural overlaps, and for preparing the plots of conserved H-bond graphs. For membrane proteins, using structures oriented along the membrane normal is recomennded, as this will facilitate interpretation of PCA-projected H-bond graphs.
* **Minimum sequence identity (%):** Structures from the selected PDB folder, whose sequence identity relative to the reference is below the minimum value set here will be discarded from the analyses.
Results are saved in the __workfolder__ created in the selected PDB folder.

#### Water cluster analysis with DBSCAN:
Default DBSCAN eps and RMS thresholds parameters are recomenned to use, as changing theses values strongly influence the results of the water clustering algorithm.
* **DBSCAN eps:** maximum distance paramter between two water molecules to be considered as in the neighborhood of eahc other. For the water clustering the default and recommended eps is 1.4 [[Bertalan, 2020]](https://www.sciencedirect.com/science/article/pii/S1047847720302070).
