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

