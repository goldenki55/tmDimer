import pandas as pd
import sys
dataIn = open('dimer_design_full/antipar_TMdimer_data.txt', 'r')
import os
from Bio.PDB.Polypeptide import *
from Bio.PDB import *
import numpy as np

columns = ['PDB ID','# Matches','Ro','w0','alpha','phi1','z_off']
data = []
dataMap = {}

# Preprocess the input file
for line in dataIn:
    line = line.strip().split()
    line = [int(line[0]), int(line[2]), float(line[4]), float(line[5]), float(line[6]), int(line[7]), float(line[8])]
    vals = tuple(line[2:])
    data.append(line)

df = pd.DataFrame(data, columns=columns)
uniqVals = {}
columnMax = []
for column in columns[2:]:
    a = list(df[column].unique())
    a.sort()
    uniqVals[column] = a
    columnMax.append(len(a)-1)


for index, row in df.iterrows():
    key = []
    for column in columns[2:]:
        key.append(uniqVals[column].index(row[column]))
    dataMap[tuple(key)] = (int(row['# Matches']),int(row['PDB ID']))

# Determine local maximum by iterating through adjacencies to see if an adjacent position has more matches, indicating it is not a local maxima
localMaxes = []
for paramInds in dataMap:
    localMax = True
    matchCount = dataMap[paramInds][0]
    paramList = list(paramInds)
    for i,param in enumerate(paramInds):
        paramTuple = paramList.copy()
        stepCount = 0
        while(tuple(paramTuple) not in dataMap and paramTuple[i] != columnMax[i]):
            paramTuple[i] += 1
        if(tuple(paramTuple) in dataMap and dataMap[tuple(paramTuple)][0] > matchCount):
            localMax = False
            break

        paramTuple = paramList.copy()
        while (tuple(paramTuple) not in dataMap and paramTuple[i] != 0):
            paramTuple[i] -= 1
        if (tuple(paramTuple) in dataMap and dataMap[tuple(paramTuple)][0] > matchCount):
            localMax = False
            break


    if(localMax):
        localMaxes.append(tuple(paramList))


smallBundleMaxes = {}
for entry in localMaxes:
    smallBundleMaxes[entry] = dataMap[entry][0]
count = 0

sorted = [k for k, v in sorted(smallBundleMaxes.items(), key=lambda item: item[1], reverse=True)]

import sys
outFile = open(sys.argv[1], 'w')
pdbList = []
pdbDir = 'Gx6G_dimers\AP_TM_dimer_database\AP_TM_dimer_database'
for localMax in sorted:
    if(dataMap[localMax][0] < 500):
        continue
    count += 1
    pdbId = str(dataMap[localMax][1])
    pdbId = (5-len(pdbId))*'0'+pdbId
    outString = 'LOCAL MAX: '+pdbId+', '+str(dataMap[localMax][0])+' matches.'+'\n'
    pdbList.append(os.path.join(pdbDir, pdbId+'.d29bf65c881e.allbb.pdb'))
    for param,column in zip(localMax,columns[2:]):
        outString += column+': '+str(uniqVals[column][param])+'\n'
    outFile.write(outString)


#print(sorted)
parser = PDBParser()
counter = 0
for pdbBB in pdbList:

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        structure = parser.get_structure("", pdbBB)

    chainA, chainB = list(structure.get_chains())[:2]  # Get first two chains of BioPDB structure object

    adjacenicesFull = {}

    for residueA in chainA:
        resA = residueA['CA'].get_coord()
        resA_id = residueA.get_id()[1]
        for residueB in chainB:
            dist = np.linalg.norm(np.array(resA) - np.array(residueB['CA'].get_coord()))
            adjacenicesFull[(resA_id, residueB.get_id()[1])] = dist
    minVal = min(adjacenicesFull.values())
    minCut = 1.2 * minVal
    adjCut = 2.2 * minVal

    adjacenices = {}

    for residueA in chainA:
        resA = residueA['CA'].get_coord()
        resA_id = residueA.get_id()[1]
        for residueB in chainB:
            dist = np.linalg.norm(np.array(resA) - np.array(residueB['CA'].get_coord()))
            if (dist < adjCut):
                adjacenices[(resA_id, residueB.get_id()[1])] = dist

    for keyPair in adjacenices:
        if(adjacenices[keyPair] != minVal):
            continue
        counter += 1
        n_term = (keyPair[0] - 7, keyPair[1] + 7)
        c_term = (keyPair[0] + 7, keyPair[1] - 7)
        if(n_term not in adjacenices or adjacenices[n_term] > minCut):
            continue
        if(c_term not in adjacenices or adjacenices[c_term] > minCut):
            continue
        print(pdbBB+'\t'+str(keyPair[0]))
        # print(n_term)
        # print(adjacenices[n_term])
        #
        # print(c_term)
        # print(adjacenices[c_term])
        # print(min
