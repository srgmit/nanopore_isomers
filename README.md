# Nanopore Isomers

Cataloging isomers of nanopores (extended vacancy defects) in 2D lattices

> Govind Rajan, A.; Silmore, K.; Swett, J.; Robertson, A. W.; Warner, J. H.; Blankschtein, D.; Strano, M. S. Nature Materials, Accepted.

## Contents

All bundled codes are written in MATLAB, and have been fully checked for compatibility with MATLAB R2016b. Post-processed results are stored in the directory `catalog` under two distinct folders: `without_edge_diffusion` and `with_edge_diffusion`. 

* `generate_isomers`: This program calls the `kmc_isomers` code repeatedly in order to stochastically generate a pre-determined number of isomers. To do so, the program requires:

  * `kmc_isomers`: This program carries out a kinetic Monte Carlo (KMC) simulation to stochastically generate one nanopore of a given size in the graphene lattice, and saves the XYZ file and the directed adjacency matrix of each generated nanopore.

* `analyze_isomers_directed_isomorphism`: This program analyzes all the nanopore isomers generated using the `generate_isomers` code to weed the duplicate isomers out of the list, and output the unique isomers. To do so, the program requires:

  * `add_weighing_nodes_in_between`: This program modifies the antimolecule adjacency matrix to add fictitious nodes to enable proper detection of C-C bonds with various orientations (see Nat. Mater. paper)

  * `compare_pores_with_weights`: This program compares two nanopore isomers, after fictitious nodes have been added into their adjacency matrices using `add_weighing_nodes_in_between`, and outputs whether the two nanopores are identical or not. 

