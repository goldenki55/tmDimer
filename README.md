## Installation
# Setting up the MASTER database - 
First, install MASTER (Zhou and Grigoryan, 2015) from the link below. 
https://grigoryanlab.org/master/
The attached design script runs MASTER as a subprocess, so make sure script usage has same executable permissions as however you run MASTER. 

Follow their provided install instructions to create a database of input backbones to search over when doing design comparison.
```
./createPDS --type target target --pdbList [List of pdb file locations] --pdsList [Output file location]
```

Next, update the provided file designPARAMS to include the absolute file location of the ```createPDS``` and ```master``` executables, and the newly created pdslist. Simply add the path names 
The designPARAMS file must always be in the same directory as the designByDimericInteractions.py script. 
biopython v1.81, numpy v1.23.5
