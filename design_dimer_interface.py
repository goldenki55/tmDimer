from Bio.PDB.Polypeptide import *
from Bio.PDB import *
import os, sys, math, copy, random
import subprocess
import numpy as np

# July 5, 2023 version


parser = PDBParser()
io = PDBIO()

resFileHeader = '''# header - applied to all residues not specifies in resfile body

USE_INPUT_SC
NATAA

# body of resfile - apply residue specific commands

start

###  interacting chains & residues'''




def getAdjMap(inputStructFile):
    struct_temp = parser.get_structure("", inputStructFile)

    chainA, chainB = list(struct_temp.get_chains())

    adjacenicesFull = {}
    for residueA in chainA:
        resA = residueA['CA'].get_coord()
        resA_id = residueA.get_id()[1]
        for residueB in chainB:
            dist = np.linalg.norm(np.array(resA) - np.array(residueB['CA'].get_coord()))
            adjacenicesFull[(resA_id, residueB.get_id()[1])] = dist
    minVal = min(adjacenicesFull.values())
    maxDistCutOff = minVal + 2.6

    adjacenices = {}

    for residueA in chainA:
        resA = residueA['CA'].get_coord()
        resA_id = residueA.get_id()[1]
        for residueB in chainB:
            dist = np.linalg.norm(np.array(resA)-np.array(residueB['CA'].get_coord()))
            if(dist < maxDistCutOff):
                adjacenices[(resA_id, residueB.get_id()[1])] = dist

    minVal = min(adjacenices.values())
    resConnections = {}
    for keyPair in adjacenices:
        resConnections[keyPair] = (((adjacenices[keyPair]-minVal)/maxDistCutOff)-math.sqrt(2))**2

    return resConnections

# Run Master, Method of Accelerated Search for Tertiary Ensemble Representative, and add the output to a temporary location
def master(queryPDB, rmsdCut, workingDir, seqFile, targetlist='masterDB_by_chain/pds_loc.txt', exePathLocation='./'):
    queryPDS = os.path.join(workingDir, 'query.pds')
    subprocess.run([os.path.join(exePathLocation, 'createPDS'), '--type', 'query', '--pdb', queryPDB, '--pds', queryPDS,
                    '--cleanPDB'])
    subprocess.run(
        [os.path.join(exePathLocation, 'master'), '--query', queryPDS, '--targetList', targetlist, '--rmsdCut', str(rmsdCut),
         #'--matchOut', os.path.join(workingDir, 'match.match'),
         '--seqOut', seqFile,
         # '--structOut',os.path.join(workingDir,'struct/')
         ])

# Cuts a two helix interface into fragLen length fragments, each showcasing an interaction centered about an
# interaction, and returns the file of the cut fragment. Must be provided the BioPython residue object of the residue
# to center the cut about. This implementation for SYMMETRIC, PARALLEL, a-helix TM bundles.
# Internal method, uses instantiated structure, pdbio
def cutPDB(chainA, chainB, fragLen, central_res_A, central_res_B, outPath):
    flank_len = int(np.floor(fragLen/2))
    # Get the index in chain of the selected residue and surrounding residues and add to set, if insufficient available
    # only add available
    central_ind = list(chainA).index(central_res_A)
    a_set = list(chainA)[max(0, central_ind-flank_len) : min(len(list(chainA)), central_ind+flank_len+1)]
    central_ind = list(chainB).index(central_res_B)
    b_set = list(chainB)[max(0, central_ind-flank_len) : min(len(list(chainB)), central_ind+flank_len+1)]

    fragment_set = set(a_set).union(set(b_set))
    class StructSelect(Select):
        def accept_residue(self, residue):
            if residue in fragment_set:
                return True
            else:
                return False
    io.save(outPath, StructSelect())


# Takes a two chain interface
# return a dataframe which is initiated with positions labeled and has spots for storing position freq/prob,
# assumed adjacencies by a-g lettering w/ their union freq/prob and position-position covariance.
# Internal method, uses instantiated structure
def initStruct(chainA, chainB, resConnections, centralG, dir, fragLen=9, rmsd=1.0, pseudoCount=1):
    AsideCenters = [centralG-18, centralG-14, centralG-11, centralG-7, centralG-4, centralG, centralG+3, centralG+7, centralG+10, centralG+14, centralG+17]
    AsideCenters = [val for val in AsideCenters if (val > 0 and val <= 28)]
    BsideCenters = AsideCenters[::-1]

    fragmentPDBs = []
    ids = []
    for aside,bside in zip(AsideCenters,BsideCenters):

        id = 'A'+str(aside)+'_B'+str(bside)
        ids.append(id)
        output = os.path.join(dir, 'cutFrag_'+id+'.pdb')
        cutPDB(chainA, chainB, fragLen, chainA.__getitem__(aside), chainB.__getitem__(bside), output)
        fragmentPDBs.append(output)

    seqFiles = []
    # Run master
    for pdb,id in zip(fragmentPDBs,ids):
        seqFile = os.path.join(dir,id+'.seq')
        seqFiles.append(seqFile)
        master(pdb,rmsd,dir,seqFile)

    # Takes the chains and master sequence outputs and counts the frequency for each position and pairwise-joint frequency
    # of adjacent residues. Map has keys of (chain A/B, residue Index) which returns a map which uses residue types as key
    # and returns frequency. As well as a joint frequency map which takes keys (resInd(A), resInd(B)) and returns a map
    # of all possible residue pairs with keys (resA, resB).
    freqMap = {}
    totalPosCount = {}
    jointFreqMap = {}
    totalJointPosCount = {}
    aminoAcids = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','Y','W']
    for residue in chainA:
        posMap = {}
        for acid in aminoAcids:
            posMap[acid] = pseudoCount
        freqMap[('A',residue.get_id()[1])] = posMap
        totalPosCount[('A',residue.get_id()[1])] = 20*pseudoCount
    for residue in chainB:
        posMap = {}
        for acid in aminoAcids:
            posMap[acid] = pseudoCount
        freqMap[('B',residue.get_id()[1])] = posMap
        totalPosCount[('B',residue.get_id()[1])] = 20*pseudoCount

    pairwisePosMap = {}
    for Aacid in aminoAcids:
        for Bacid in aminoAcids:
            pairwisePosMap[(Aacid,Bacid)] = pseudoCount
    for keyPair in resConnections:
        jointFreqMap[keyPair] = copy.deepcopy(pairwisePosMap)
        totalJointPosCount[keyPair] = 400*pseudoCount


    for pdb,seqFile in zip(fragmentPDBs, seqFiles):
        fragment = parser.get_structure("",pdb)
        keys = []
        pairwiseKeys = {}
        chainA_endInd = 0
        for chain in fragment.get_chains():
            for i,residue in enumerate(chain):
                keys.append((chain.get_id(), residue.get_id()[1]))
            if(chainA_endInd == 0):
                chainA_endInd = i
        keys = keys[1:chainA_endInd-2]+keys[chainA_endInd+2:-1]             # Remove the most flanking residues from consideration
        for key in jointFreqMap:
            try:
                # If the two residues exist in the fragment, add to pairwise keys list
                chainA, chainB = fragment.get_chains()
                chainA.__getitem__(key[0])
                chainB.__getitem__(key[1])
                indA = 0
                indB = 0
                for ind,val in enumerate(keys):
                    if(val == ('A',key[0])):
                        indA = ind
                    if(val == ('B',key[1])):
                        indB = ind
                pairwiseKeys[key] = (indA,indB)
            except:
                pass
        scaling = math.exp(-1*rmsd)
        with open(seqFile) as seqFile:
            for line in seqFile:
                line = line.strip().split(' ')
                weight = math.exp(-1*float(line[0]))/scaling
                seq = []
                for res in line[1:]:
                    try:
                        seq.append(three_to_one(res))
                    except:
                        seq.append('X')
                seq = seq[1:chainA_endInd-2]+seq[chainA_endInd+2:-1]
                if(len(seq) != len(keys)):
                    print("Mismatch in fragment length with master found seq")

                for res,key in zip(seq,keys):
                    if(res == 'X'):
                        continue
                    freqMap[key][res] += weight
                    totalPosCount[key] += weight
                for key in pairwiseKeys:
                    resA = seq[pairwiseKeys[key][0]]
                    resB = seq[pairwiseKeys[key][1]]
                    if(resA == 'X' or resB == 'X'):
                        continue
                    jointFreqMap[key][(resA,resB)] += weight
                    totalJointPosCount[key] += weight


    for key in freqMap:
        for res in freqMap[key]:
            freqMap[key][res] = freqMap[key][res]/totalPosCount[key]


    # Update jointFreqMap to accommodate for homo-antiPar
    new_jointFreqMap = {}
    new_totalJointPosCount = {}
    for key in jointFreqMap:
        new_jointFreqMap[key] = {}
    for keyPair in jointFreqMap:
        reverseKey = (keyPair[1], keyPair[0])
        for resPair in jointFreqMap[keyPair]:
            new_jointFreqMap[keyPair][resPair] = jointFreqMap[reverseKey][resPair]+jointFreqMap[keyPair][resPair]
        new_totalJointPosCount[keyPair] = totalJointPosCount[reverseKey]+totalJointPosCount[key]
    jointFreqMap = new_jointFreqMap
    totalJointPosCount = new_totalJointPosCount

    for key in jointFreqMap:
        for resPair in jointFreqMap[key]:
            jointFreqMap[key][resPair] = jointFreqMap[key][resPair]/totalJointPosCount[key]


    return freqMap, jointFreqMap, AsideCenters, BsideCenters


# Initial conditions should be a list of pairs, (position, Letter)
def buildMap(inputSeq, posProbMap, jointProbMap, resConnections, fixedPos, centralG):
    TM_aaFrequency = {'A': 0.1102, 'R': 0.0256, 'N': 0.0225, 'D': 0.0135, 'C': 0.0141, 'Q': 0.0168, 'E': 0.0178,
                      'G': 0.0752, 'H': 0.0102, 'I': 0.1007, 'L': 0.1525, 'K': 0.0203, 'M': 0.0342, 'P': 0.0298,
                      'F': 0.0824, 'S': 0.0544, 'T': 0.0519, 'W': 0.0259, 'Y': 0.0409, 'V': 0.1005}
    aminoAcids = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','Y','W']

    adjacencies = {}

    for keyPair in jointProbMap:
        # if(keyPair[0] not in range(5,24) or keyPair[1] not in range(5,24)):
        #     continue
        if(keyPair[0] not in adjacencies):
            adjacencies[keyPair[0]] = [keyPair[1]]
        else:
            adjacencies[keyPair[0]].append(keyPair[1])
        if(keyPair[1] not in adjacencies):
            adjacencies[keyPair[1]] = [keyPair[0]]
        else:
            adjacencies[keyPair[1]].append(keyPair[0])

    inputSeq = [a for a in inputSeq]
    for val in adjacencies:
        adjacencies[val] = list(set(adjacencies[val]))
        inputSeq[val-1] = 'A'

    for position in fixedPos:
        inputSeq[centralG+position-1] = fixedPos[position]

    fixed = []
    for val in fixedPos:
        fixed.append(centralG+val)
    fixedPos = set(fixed)
    inputSeq = ''.join(inputSeq)

    edgeValues = {}
    for position in jointProbMap:
        temp = {}

        for aaPair in jointProbMap[position]:
            a = jointProbMap[position][aaPair]
            temp[aaPair] = resConnections[position]*a*math.log2(a/(TM_aaFrequency[aaPair[0]]*TM_aaFrequency[aaPair[1]]))

        edgeValues[position] = temp

    nodeValues = {}
    for position in posProbMap:
        temp = {}

        for aa in posProbMap[position]:
            a = posProbMap[position][aa]
            temp[aa] = a*math.log2(a/TM_aaFrequency[aa])
        nodeValues[position[1]] = temp

    states = {}
    initialSeq = {}
    for position in nodeValues:
        if(position not in adjacencies):
            continue
        initialSeq[position] = inputSeq[position-1]

    # Defines the numerical value of a specificed residue at a specified position in the sequence
    def posValue(position, positionRes):
        value = 0.5*nodeValues[position][positionRes]
        if(position in adjacencies):
            for adj in adjacencies[position]:
                if((position, adj) in jointProbMap):
                    value += edgeValues[(position, adj)][(positionRes, states[adj])]
                elif((adj, position) in jointProbMap):
                    value += edgeValues[(adj, position)][(states[adj], positionRes)]
                else:
                    print("ERROR, BAD POSITION")
        return value

    def evalSeq():
        seqValue = 0
        seq = list(inputSeq)
        for key in states:
            seq[key-1] = states[key]
            seqValue += posValue(key, states[key])
        return ''.join(seq), seqValue

    def genResFile(seq, delta):
        ppm = []
        skipRow = [0.05]*20
        # Undefined position
        undefPos = " A NATAA"

        resFile = resFileHeader + '\n'
        for i in range(1, len(seq) + 1):
            if (i not in adjacencies):
                resFile += str(i)
                resFile += undefPos + '\t\t\t\t\t# '  + '\n'
                ppm.append(skipRow)
                continue
            elif (i in fixedPos):
                resFile += str(i) + " A PIKAA " + seq[i - 1] + '\t\t\t\t\t#'  + ' fixed\n'
                ppm.append(skipRow)
                continue

            aaWeights = {posValue(i, a): a for a in aminoAcids}
            norm_weights = np.array(list(aaWeights.keys()))
            norm_weights[norm_weights < 0] = 0
            norm_weights = norm_weights/sum(norm_weights)
            ppm.append(np.array(list(norm_weights)))
            maxVal = max(aaWeights.keys())
            minAccepted = maxVal
            maxNotAccepted = 0
            acceptedAAs = ''
            for aaVal in aaWeights:
                #if (abs(maxVal - aaVal) < delta):
                if(aaVal > delta):
                    acceptedAAs += aaWeights[aaVal]
                    minAccepted = min(minAccepted, aaVal)
                else:
                    maxNotAccepted = max(maxNotAccepted, aaVal)
            resFile += str(i)
            tabs2add = '\t' * (5 - int(len(acceptedAAs) / 6))
            resFile += ' A PIKAA ' + acceptedAAs + tabs2add + '#'  + ' |'
            resFile += ' maxAccept: ' + aaWeights[maxVal] + ' - ' + str(round(maxVal, 2)) + ', minAccept: ' + \
                       aaWeights[minAccepted] + ' - ' + str(round(minAccepted, 2))
            if (maxNotAccepted != 0):
                resFile += ', maxReject: ' + aaWeights[maxNotAccepted] + ' - ' + str(
                    round(maxNotAccepted, 2)) + '\n'
            else:
                resFile += '\n'
        ppm = np.transpose(np.asarray(ppm))
        return resFile, ppm

    dataCollection = []
    seq2rosettaDesign = {}
    maxSeqVal = 0
    ppm = []
    # Iterations from input sequence
    for i in range(1000):
        states = initialSeq
        for j in range(100):
            pos2update = random.choice(list(states.keys()))
            if(pos2update not in fixedPos):
                # Choose based on position dependent weights of upregulated values
                newRes = random.choices(list(nodeValues[pos2update].keys()), weights=list(nodeValues[pos2update].values()))[0]
                states[pos2update] = newRes

            for adj in adjacencies[pos2update]:

                if(adj in fixedPos):
                    continue
                # Choose adjacency based on all edge dependent nodeValues
                aaWeights = [posValue(adj, a) for a in aminoAcids]
                aaWeights = [max(1e-5,a) for a in aaWeights]

                newRes = random.choices(aminoAcids, weights=aaWeights)[0]

                states[adj] = newRes

            if (pos2update not in fixedPos):
                # Choose best value based on edges
                # newRes = states[pos2update]
                # maxVal = posValue(pos2update, newRes)
                # for aa in aminoAcids:
                #     newVal = posValue(pos2update, aa)
                #     if(newVal > maxVal):
                #         maxVal = newVal
                #         newRes = aa

                # Choose adjaceny based on all edge dependent nodeValues
                aaWeights = [posValue(adj, a) for a in aminoAcids]
                aaWeights = [max(1e-5,a) for a in aaWeights]

                newRes = random.choices(aminoAcids, weights=aaWeights)[0]

                states[pos2update] = newRes


        seq, val = evalSeq()
        seq2rosettaDesign[(seq,val)], ppm_cand = genResFile(seq, 0.02)
        if(val > maxSeqVal):
            ppm = ppm_cand
            maxSeqVal = val
        if(val > 1):
            dataCollection.append((seq, val))



    dataCollection = list(set(dataCollection))
    dataCollection.sort(key = lambda x: x[1], reverse=True)



    return dataCollection, seq2rosettaDesign, ppm


if __name__ == '__main__':
    # Usage - python designDimer.py [pdb] [outputDir] [outputFileName w/ scores]
    struct = sys.argv[1]  # Input pdb
    dir = sys.argv[2]
    outPrefix = sys.argv[3]
    seqOut = os.path.join(dir,outPrefix+'.fasta')
    jointProbOut = os.path.join(dir,outPrefix+'_jointProb.txt')
    rosettaOut = os.path.join(dir,outPrefix+'_resFile.txt')
    aminoAcids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'Y', 'W']

    TM_aaFrequency = {'A': 0.1102, 'R': 0.0256, 'N': 0.0225, 'D': 0.0135, 'C': 0.0141, 'Q': 0.0168, 'E': 0.0178,
                      'G': 0.0752, 'H': 0.0102, 'I': 0.1007, 'L': 0.1525, 'K': 0.0203, 'M': 0.0342, 'P': 0.0298,
                      'F': 0.0824, 'S': 0.0544, 'T': 0.0519, 'W': 0.0259, 'Y': 0.0409, 'V': 0.1005}
    # RUN
    structure = parser.get_structure("", struct)
    io.set_structure(structure)
    chainA, chainB = list(structure.get_chains())[:2]  # Get first two chains of BioPDB structure object

    #inputSeq = "XXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    inputSeq = "XXXXAGXXAAXAGXXAAXAGXXAAXXXX"
    #inputSeq = "GGGGGGGGGGGGGGGGG"

    fixedPos = {0:'G', -7:'G', 7:'G'}
    #fixedPos = {}
    centralG = 13

    resConnections = getAdjMap(struct, dir)
    posProbMap,jointProbMap,acenters,bcenters = initStruct(chainA,chainB,resConnections,centralG,dir)
    out, rosetta, ppm = buildMap(inputSeq, posProbMap, jointProbMap, resConnections, fixedPos, centralG)

    with open(seqOut, 'w') as outFile:
        for entry in out:
            outFile.write('>'+str(round(entry[1],2))+'\n'+entry[0]+'\n')
    with open(rosettaOut, 'w') as outFile:
        for entry in out:
            outFile.write('>'+entry[0]+', '+str(round(entry[1],2))+'\n'+rosetta[entry]+'\n')


    indepJointProbs = {}
    for X in TM_aaFrequency:
        for Y in TM_aaFrequency:
            indepJointProbs[(X, Y)] = 2 * TM_aaFrequency[X] * TM_aaFrequency[Y]

    positions = []
    for key in jointProbMap:
        positions.append(key)
    positions.sort()
    with open(jointProbOut, 'w') as outFile:
        for key in positions:

            outFile.write("Position: " + str(key) + '\n')
            values = []

            for resPair in indepJointProbs:
                # if((key[0] in hexEnrichChangeMap and hexEnrichChangeMap[key[0]][resPair[0]] < -0.1) or (key[1] in hexEnrichChangeMap and hexEnrichChangeMap[key[1]][resPair[1]] < -0.1)):
                #     continue
                prob = jointProbMap[key][resPair]

                # logChange = math.log2(prob / pentjointProbMat[key][resPair])
                logChange = math.log2(prob / indepJointProbs[resPair])
                enrichment = prob * logChange
                values.append((prob, enrichment, logChange, str(resPair)))
            values.sort(reverse=True)
            for val in values:
                if (val[0] < 0.000001):
                    break



                outFile.write(val[3] + ": " + str(round(val[0], 3)) + ', ' + str(round(val[1], 3)) + ', ' + str(
                    round(val[2], 3)) + '\n')


