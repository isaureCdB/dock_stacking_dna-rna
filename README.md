# dock_stacking_dna-rna

###################################  
    Context  
###################################

We want to dock a ssDNA or ssRNA sequence on a protein structure/model, and we assume that some known aromatic amino-acids do stacking interactions with some nucleotide bases.

###################################  
  Approach  
###################################

We use :
- a fragment-based approach to tackle the huge conformational space of ssDNA/ssRNA (thereafyter called ssNA)
- home made libraries* of ssNA
- a coarse-grain representation of the protein and ssNA
- many random positions of the fragments around the protein 
- a gradient descent* minimisation of the interaction energy in a force-field
- distance restraints between the aromatic ring and base ring atoms

* see:
- tool to create libraries: https://doi.org/10.1093/bioinformatics/btac430
- dna library: https://zenodo.org/record/6517064#.YsV9tHjP1H4

** Monte-carlo approach is also available, with some script tuning

###################################  
  Procedure  
###################################

List the anchoring residue numbers (from your initial protein PDB file) in aromatics.list.
Create one list per cluster of neighbor anchors. For your tetramer, run the docking once for each monomer, by giving for each run the list of stacking aromatics in that monomer. 

!!! Each residue number must be present only once in your protein. If you have a multi-chains protein, renumber the residues using renum-at-res-nochain.py.

List the 3-nt motifs you want to dock in motifs.list, one motif per line. Add each motif only once in the list, even if it is present several times in the full sequence. ex: for AUUUUUG, the list is:
AUU
UUU
UUG

Run:
> dock-dna-anchors.sh protein.pdb aromatics.list [nb of CPU] [nb of pdb models*]

*this is just for visualisation of the few top-ranked models.

If you have a reference dna with which you want to compute RMSDs, call it dna.pdb. It must contain only one continuous chain.

----------------------------
create_rest.py
----------------------------
This script defines a maximal distance between each of the 3 coarse-grain beads of the base hexagonal cycle and the bead at the center of the amino-acid cycle.
Allow some wooble in the stacking: Use a max distance larger than the ~3.6A max distance observed in Xray structures.

The max distance restraints are writen in textfiles s.a. restraints/rest1.txt, and used by option --rest in the docking. See attract/restraints.txt for more details.
----------------------------

As we don't know which of the 3 nucleotides of a fragment stack to which aromatics, dock-dna-achors.sh dock with all possible combinations, in a for loop. The restraints files corresponding to each situation are listed in restraints.list.

With restraints, one doesn't need a very large number of starting positions around the receptor.
The restraints will orientate the ligand toward the stacking base at the start of the docking, in the "ghost only-rot" mode:
The protein-RNA force-field (ff) is disabled so that the ligand "feels" only the restraints, and only rotations of ligand and receptor are allowed.
Afterward, those constraints are released and both ff and restraints are on, with both rotations and translations allowed.

The resulting orientation of the ghost only-rot mode still depends on the starting orientation and position, so use a few of them.
The current script uses 100. To know if it is enough/overkill, use deredundant on the output of the ghost only-rot minimisation, and see if some/most poses get eliminated:

> $ATTRACTDIR/deredundant ${motif}-1.dat 2 --ens 0 $Nconf --lim 0.2 > test.dat
> tail -n 3  test.dat|grep "#"
> tail -n 3  ${motif}-0.dat|grep "#"
