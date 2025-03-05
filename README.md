# Design Script for Gx6G membrane protein building block

## Installation
First, install MASTER (Zhou and Grigoryan, 2015) from the link below. 
https://grigoryanlab.org/master/
The attached design script runs MASTER as a subprocess, so make sure script usage has same executable permissions as however you run MASTER. 

Follow their provided install instructions to create a database of input backbones to search over when doing design comparison.
```
./createPDS --type target target --pdbList [List of pdb file locations] --pdsList pdsList.txt
```

## 
Next, update the provided file designPARAMS to include the absolute file location of the ```createPDS``` and ```master``` executables, and the newly created pdslist. Simply add the path names after the respective locations as:
```
CREATE_PDS_LOC [fullpath]/createPDS
RUN_MASTER_LOC [fullpath]/master
PDSLIST_DATABSE [fullpath]/pdsList.txt
```

The designPARAMS file must always be in the same directory as the designByDimericInteractions.py script. It also includes option to adjust parameters used in running the design script. The default parameters reflect what was used in "Design prcinpiples of the common Gly-X6-Gly membrane protein building block" including the background amino acid frequency table and number of iterations per backbone.

Besides MASTER, the script also uses some python libraries, including: 

biopython v1.81

numpy v1.23.5

both of which can be installed with pip and should be fine at the versions listed or more recent versions. 

## Running

The design script can be run on either one individual backbone or a database of backbones. To enforce symmetry, residues with the same index in the PDB (regardless of chain) are treated as being forced to have the same amino acid assignment. Additionally, all residue assignments in the input structure are overwritten. To run on one backbone, you can use
```
python designByDimericInteractions.py -p designPDB -in 01322.d29bf65c881e.allbb_polyA.pdb -c AB -o output/
```

Or run on a directory containing many different backbones as
```
python designByDimericInteractions.py -p designPDBDatabase -in backbones/ -c AB -o output/
```

