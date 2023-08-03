from Bio.PDB.Polypeptide import *
from Bio.PDB import *
import os, math, copy, random, seqlogo
import subprocess
import numpy as np
import warnings
# August 2023 version, J. Golden, Mravic Lab, Scripps Research
# ASSUMES THAT RESIDUES WITH THE SAME NUMBER ARE THE SAME RESIDUE INDEPENDENT OF CHAIN
# Please renumber residues if you want them to not be the same residue


def parsePARAMS(paramFile):
    masterCreatePDSPath = None
    masterPath = None
    masterDBPDSList = None
    masterRMSD = None
    outputMasterFiles = None
    fragLength = None
    deltaDistCut = None
    positionRelWeight = None
    sequenceCompRelWeight = None
    psuedoCount = None
    sequenceIterations = None
    sequenceMutations = None
    outSeqTables = None
    bg_AA_freq = {}
    with open(paramFile,'r') as f:
        for line in f:
            if(line[0] == '#' or line == '\n'): continue
            line = line.strip().split()
            if(len(line) < 2):
                raise ValueError("designPARAM input not understood: "+' '.join(line))
            if(line[0] == 'BG_AA_FREQ'):
                if(len(line) < 3):
                    raise ValueError("designPARAM input not understood: "+' '.join(line))
                bg_AA_freq[line[1]] = float(line[2])
            elif(line[0] == 'CREATE_PDS_LOC'):
                masterCreatePDSPath = os.path.normpath(line[1])
            elif(line[0] == 'RUN_MASTER_LOC'):
                masterPath = os.path.normpath(line[1])
            elif(line[0] == 'PDSLIST_DATABASE'):
                masterDBPDSList = os.path.normpath(line[1])
            elif(line[0] == 'RMSD_CUT'):
                masterRMSD = float(line[1])
            elif(line[0] == 'OUTPUT_MASTER_FILES'):
                outputMasterFiles = bool(line[1])
            elif(line[0] == 'CUT_FRAG_LEN'):
                fragLength = int(line[1])
            elif(line[0] == 'DEL_DIST_CUT'):
                deltaDistCut = float(line[1])
            elif(line[0] == 'PSU_COUNT'):
                psuedoCount = float(line[1])
            elif(line[0] == 'POSITION_SCORE_WEIGHT'):
                positionRelWeight = float(line[1])
            elif(line[0] == 'SEQUENCECOMP_SCORE_WEIGHT'):
                sequenceCompRelWeight = float(line[1])
            elif(line[0] == 'SEQ_ITERS'):
                sequenceIterations = int(line[1])
            elif(line[0] == 'SEQ_MUTATIONS'):
                sequenceMutations = int(line[1])
            elif(line[0] == 'OUTPUT_SEQ_TABLES'):
                outSeqTables = bool(line[1])
            else:
                raise ValueError("designPARAM input not understood: " + ' '.join(line))
    return masterCreatePDSPath, masterPath, masterDBPDSList, masterRMSD, fragLength, deltaDistCut, positionRelWeight, \
        sequenceCompRelWeight, psuedoCount, outputMasterFiles, sequenceIterations, sequenceMutations, outSeqTables, bg_AA_freq

_parent_dir = os.path.dirname(__file__)
if os.path.exists(os.path.join(_parent_dir, "designPARAMS")):
    masterCreatePDSPath, masterPath, masterDBPDSList, masterRMSD, fragLength, deltaDistCut, positionRelWeight, sequenceCompRelWeight,\
        psuedoCount, outputMasterFiles, sequenceIterations, sequenceMutations, outSeqTables, bg_AA_freq = parsePARAMS(os.path.join(_parent_dir, "designPARAMS"))
    # Throw exception if parameter is not provided for required parameter
    if(masterCreatePDSPath == None or os.path.exists(masterCreatePDSPath) == False):
        raise ValueError("Please specify master createPDS executable location in designPARAMS")
    if(masterPath == None or os.path.exists(masterPath) == False):
        raise ValueError("Please specify master executable location in designPARAMS")
    if(masterDBPDSList == None or os.path.exists(masterDBPDSList) == False):
        raise ValueError("Please specify master database PDS list location in designPARAMS")
    # If parameter is not provided for non-required parameters, set to default
    if(masterRMSD == None):
        masterRMSD = 1.0
    if(fragLength == None):
        fragLength = 9
    if(deltaDistCut == None):
        deltaDistCut = 2.6
    if(positionRelWeight == None):
        positionRelWeight = 0.5
    if(sequenceCompRelWeight == None):
        sequenceCompRelWeight = 0.3
    if(psuedoCount == None):
        psuedoCount = 1.0
    if(outputMasterFiles == None):
        outputMasterFiles = False
    if(sequenceIterations == None):
        sequenceIterations = 500
    if(sequenceMutations == None):
        sequenceMutations = 200
    if(outSeqTables == None):
        outSeqTables = False
    if(len(bg_AA_freq) == 0):
        warnings.warn("Background amino acid frequency not provided, using default TM background", UserWarning)
        bg_AA_freq = {'A': 0.1102, 'R': 0.0256, 'N': 0.0225, 'D': 0.0135, 'C': 0.0141, 'Q': 0.0168, 'E': 0.0178,
                      'G': 0.0752, 'H': 0.0102, 'I': 0.1007, 'L': 0.1525, 'K': 0.0203, 'M': 0.0342, 'P': 0.0298,
                      'F': 0.0824, 'S': 0.0544, 'T': 0.0519, 'W': 0.0259, 'Y': 0.0409, 'V': 0.1005}
    elif(bg_AA_freq.keys != set("ARNDCQEGHILKMPFSTWYV")):
        raise ValueError("Input background amino acid frequency is incomplete")
else:
    raise FileNotFoundError("designPARAMS file not found in script directory")

# Provide a list of Biopython.chain objects to use as basis for adj map
# deltaDistCut describes the difference between the minimum Ca-Ca distance across the helix interface
# from the input model and the allowed maximum. For a close helix-helix interface should be ~2.6, the distance of
# a helix radius. For other interactions, may be useful to raise this
def getAdjMap(inputChains, deltaDistCut = deltaDistCut):
    # Get list of Ca-Ca distances across helix-helix interface to find minimum distance
    minDist = np.inf
    for i in range(len(inputChains)):
        chainA = inputChains[i]
        for j in range(i+1, len(inputChains)):
            chainB = inputChains[j]
            for resA in chainA:
                resACa = np.array(resA['CA'].get_coord())
                for resB in chainB:
                    resBCa = np.array(resB['CA'].get_coord())
                    minDist = min(minDist, np.linalg.norm(resACa - resBCa))
    maxDistCutOff = minDist + deltaDistCut

    adjacenices = {}
    for i in range(len(inputChains)):
        chainA = inputChains[i]
        for j in range(i+1, len(inputChains)):
            chainB = inputChains[j]

            for residueA in chainA:
                resACa = residueA['CA'].get_coord()
                for residueB in chainB:
                    dist = np.linalg.norm(np.array(resACa)-np.array(residueB['CA'].get_coord()))
                    if(dist < maxDistCutOff):
                        adjacenices[(residueA, residueB)] = dist


    resConnections = {}     # List of interaction residue pairs that maps to normalized distance weight
                            # The smallest interaction is always weighted at 2, the seperation for larger distances becomes
                            # more significant with more separation
                            # This weighting (both the function itself and deltaCutoff) can definitely be optimized
    for keyPair in adjacenices:
        resConnections[keyPair] = (((adjacenices[keyPair]-minDist)/maxDistCutOff)-math.sqrt(2))**2

    return resConnections

# Run Master, Method of Accelerated Search for Tertiary Ensemble Representative, and add the output to a temporary location
def master(queryPDB, workingDir, seqFile, rmsdCut, targetlist=masterDBPDSList, createPDSLocation=masterCreatePDSPath, runMasterLocation=masterPath):
    queryPDS = os.path.join(workingDir, 'query.pds')
    subprocess.run([createPDSLocation, '--type', 'query', '--pdb', queryPDB, '--pds', queryPDS,
                    '--cleanPDB'], stdout=subprocess.DEVNULL)
    subprocess.run(
        [runMasterLocation, '--query', queryPDS, '--targetList', targetlist, '--rmsdCut', str(rmsdCut),
         #'--matchOut', os.path.join(workingDir, 'match.match'),
         '--seqOut', seqFile,
         # '--structOut',os.path.join(workingDir,'struct/')
         ], stdout=subprocess.DEVNULL)
    os.remove(queryPDS)

# Cuts a two helix interface into fragLen length fragments, each showcasing an interaction centered about an
# interaction, and returns the file of the cut fragment. Must be provided the BioPython residue object of the residue
# to center the cut about. This implementation for SYMMETRIC, PARALLEL, a-helix TM bundles.
# Internal method, uses instantiated structure, pdbio
def cutPDB(fragLen, central_res_A, central_res_B, outPath, io):
    flank_len = int(np.floor(fragLen/2))
    # Get the index in chain of the selected residue and surrounding residues and add to set, if insufficient available
    # only add available
    chainA = central_res_A.get_parent()
    chainB = central_res_B.get_parent()
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
# resConnections has the residues ORDERED, and chain specific
def initStruct(resConnections, workingDir, parser, rmsd=masterRMSD, fragLen=fragLength, pseudoCount=psuedoCount, silent=outputMasterFiles):
    AsideCenters = []
    BsideCenters = []
    for keyPair in resConnections:
        AsideCenters.append(keyPair[0])
        BsideCenters.append(keyPair[1])
    fragmentPDBs = []
    ids = []
    io = PDBIO()
    for aside,bside in zip(AsideCenters,BsideCenters):

        id = 'A'+str(aside.get_id()[1])+'_B'+str(bside.get_id()[1])
        ids.append(id)
        output = os.path.join(workingDir, 'cutFrag_' + id + '.pdb')
        cutPDB(fragLen, aside, bside, output, io)
        fragmentPDBs.append(output)

    seqFiles = []
    # Run master
    for pdb,id in zip(fragmentPDBs,ids):
        seqFile = os.path.join(workingDir, id + '.seq')
        seqFiles.append(seqFile)
        master(pdb, workingDir, seqFile, masterRMSD)

    # Takes the chains and master sequence outputs and counts the frequency for each position and pairwise-joint frequency
    # of adjacent residues. Map has keys of (chain A/B, residue Index) which returns a map which uses residue types as key
    # and returns frequency. As well as a joint frequency map which takes keys (resInd(A), resInd(B)) and returns a map
    # of all possible residue pairs with keys (resA, resB).
    freqMap = {}
    jointFreqMap = {}


    aminoAcids = list("ACDEFGHIKLMNPQRSTVYW")
    posMap = {'count': 20*pseudoCount}
    for acid in aminoAcids:
        posMap[acid] = pseudoCount
    for keyPair in resConnections:
        for key in keyPair:
            key = key.get_id()[1]
            if(key not in freqMap):
                freqMap[key] = copy.deepcopy(posMap)

    pairwisePosMap = {'count': 400*pseudoCount}
    for Aacid in aminoAcids:
        for Bacid in aminoAcids:
            pairwisePosMap[(Aacid,Bacid)] = pseudoCount
    for keyPair in resConnections:
        keyPair = (keyPair[0].get_id()[1], keyPair[1].get_id()[1])
        if(keyPair not in jointFreqMap):
            pairwiseMap = copy.deepcopy(pairwisePosMap)     # Both the key and inverseKey give the same dict object
            jointFreqMap[keyPair] = pairwiseMap
            inverseKey = (keyPair[1], keyPair[0])
            jointFreqMap[inverseKey] = pairwiseMap

    # Parse master generated sequence files and count amino acid occurences at each position and pairwise interaction
    for pdb,seqFile in zip(fragmentPDBs, seqFiles):
        fragment = parser.get_structure("",pdb)
        fragChainA, fragChainB = fragment.get_chains()
        if(silent):
            os.remove(pdb)

        chainA_endInd = 0
        keys = []
        for chain in [fragChainA, fragChainB]:
            for i, residue in enumerate(chain):
                keys.append((chain.get_id(), residue.get_id()[1]))
            if (chainA_endInd == 0):
                chainA_endInd = i
        keys = keys[1:chainA_endInd-2] + keys[chainA_endInd+2:-1]  # Remove the most flanking residues from consideration

        pairwiseKeys = {}
        for keyPair in jointFreqMap:
            try:
                # If the two residues exist in the fragment, add to pairwise keys list
                fragChainA.__getitem__(keyPair[0])
                fragChainB.__getitem__(keyPair[1])
                indA = 0
                indB = 0
                for ind, val in enumerate(keys):
                    if (val == keyPair[0]):
                        indA = ind
                    if (val == keyPair[1]):
                        indB = ind
                pairwiseKeys[keyPair] = (indA, indB)
            except:
                pass

        scaling = math.exp(-1*rmsd)
        with open(seqFile) as sFile:
            for line in sFile:
                line = line.strip().split(' ')      # Master sequence file, output is `RMSD AA1 AA2 AA3 ...`
                # RMSD deviation of found fragment weighting, this is a function that could be optimized
                weight = math.exp(-1*float(line[0]))/scaling
                seq = []
                for res in line[1:]:
                    try:
                        seq.append(three_to_one(res))
                    except:
                        seq.append('X')
                seq = seq[1:chainA_endInd-2]+seq[chainA_endInd+2:-1]        # Remove sequence elements from removed indices (N and C term residues of both fragments)
                if(len(seq) != len(keys)):
                    raise Exception("Mismatch in fragment length with master found seq")

                for res,key in zip(seq,keys):
                    if(res == 'X'):
                        continue
                    freqMap[key][res] += weight
                    freqMap[key]['count'] += weight
                for key in pairwiseKeys:
                    resA = seq[pairwiseKeys[key][0]]
                    resB = seq[pairwiseKeys[key][1]]
                    if(resA == 'X' or resB == 'X'):
                        continue
                    jointFreqMap[key][(resA,resB)] += weight
                    jointFreqMap[key]['count'] += weight
        if(silent):
            os.remove(seqFile)


    for key in freqMap:
        count = freqMap[key]['count']
        for res in freqMap[key]:
            freqMap[key][res] = freqMap[key][res]/count

    for keyPair in jointFreqMap:
        count = jointFreqMap[keyPair]['count']
        for resPair in jointFreqMap[keyPair]:
            jointFreqMap[keyPair][resPair] = jointFreqMap[keyPair][resPair]/count

    return freqMap, jointFreqMap


# Fixed position should be a map
# resConnections is ORDERED, with the keys as chain-specific Biopython.residue objects
# jointProbMap is UNORERED, with keys as indices, either ordering of keys (i, j) or (j, i) return the same probability map
# position realtive weight defines the weight of the single position sequence versus the pairwise information (weighted as 1)
def buildMap(seqLen, posProbMap, jointProbMap, resConnections, fixedPos, background_AA_frequency = bg_AA_freq,
             position_relative_weight = positionRelWeight, sequenceComp_relative_weight = sequenceCompRelWeight,
             sequenceIterations=sequenceIterations, sequenceMutations=sequenceMutations):
    aminoAcids = list("ACDEFGHIKLMNPQRSTVYW")

    adjacencies = {}       # Maps any index to all of its nearby (in 3D space) residue indicies
    for keyPair in jointProbMap:
        if(keyPair[0] not in adjacencies):
            adjacencies[keyPair[0]] = [keyPair[1]]
        else:
            adjacencies[keyPair[0]].append(keyPair[1])
        if(keyPair[1] not in adjacencies):
            adjacencies[keyPair[1]] = [keyPair[0]]
        else:
            adjacencies[keyPair[1]].append(keyPair[0])
    for keyPair in resConnections:
        indexKeyPair = (keyPair[0].get_id()[1], keyPair[1].get_id()[1])
        jointProbMap[indexKeyPair]['distance weight'] = resConnections[keyPair]

    inputSeq = list('X'*seqLen)
    for val in adjacencies:
        adjacencies[val] = list(set(adjacencies[val]))      # Remove duplicates
        inputSeq[val-1] = 'A'                               # Initiate all sequence to be designed to Ala

    for position in fixedPos:
        inputSeq[position-1] = fixedPos[position]
    fixedPos = set(fixedPos.keys())

    for position in jointProbMap.values():
        for aaPair in position:
            aaPairProb = position[aaPair]
            # Edge value function, weight*prob*log2(joint prob/independent background prob)
            position[aaPair] = position['distance weight'] \
                           * aaPairProb \
                           * math.log2(aaPairProb / (background_AA_frequency[aaPair[0]] * background_AA_frequency[aaPair[1]]))

    for position in posProbMap.values():
        for aa in position:
            aminoAcidProb = position[aa]
            # Node value function, probAA from sequence * log2(probAA from sequence/background aa prob)
            position[aa] = aminoAcidProb*math.log2(aminoAcidProb / background_AA_frequency[aa])

    initialSequenceComp = {aa: 0 for aa in aminoAcids}
    initialSequenceComp['total'] = 0
    initialSeq = {}
    for position in posProbMap:
        if(position not in adjacencies):
            continue
        res = inputSeq[position-1]
        initialSeq[position] = res
        initialSequenceComp[res] += 1
        initialSequenceComp['total'] += 1

    # Defines the numerical value of a specificed residue at a specified position in the sequence
    # This is the "score function" for choosing a residue at a specific position
    def posValue(position, positionRes):
        value = position_relative_weight*posProbMap[position][positionRes]

        residuePercent = sequenceComp[positionRes]/sequenceComp['total']
        value += sequenceComp_relative_weight*math.log2(residuePercent/background_AA_frequency[positionRes])

        if(position in adjacencies):
            for adj in adjacencies[position]:
                value += jointProbMap[(position, adj)][(positionRes, states[adj])]
        return value

    def evalSeq():
        seqValue = 0
        seq = list(inputSeq)
        for key in states:
            seq[key-1] = states[key]
            seqValue += posValue(key, states[key])
        return ''.join(seq), seqValue

    def genResFile(seq, delta):
        resFileHeader = '''# header - applied to all residues not specifies in resfile body

        USE_INPUT_SC
        NATAA

        # body of resfile - apply residue specific commands

        start

        ###  interacting chains & residues'''
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
    for i in range(sequenceIterations):
        states = copy.copy(initialSeq)
        sequenceComp = copy.copy(initialSequenceComp)
        for j in range(sequenceMutations):
            # At each sequence mutation, pick a location randomly, perform a "strong mutation" based only on the position
            # specific information, allowing for higher probability of sampling of different residues, update all
            # adjacent residues based on covariant interactions, performing "stabilizing mutations"
            # meant to benefit covariant interactions, then "stabilize" the initial position
            pos2update = random.choice(list(states.keys()))
            if (pos2update not in fixedPos):
                # Choose based on position dependent weights of upregulated values
                newRes = random.choices(list(posProbMap[pos2update].keys()), weights=list(posProbMap[pos2update].values()))[0]
                sequenceComp[states[pos2update]] -= 1
                sequenceComp[newRes] += 1
                states[pos2update] = newRes


            for adj in adjacencies[pos2update]:

                if(adj not in fixedPos):
                    # Choose adjacency based on all edge dependent nodeValues
                    aaWeights = [max(1e-5,posValue(adj, aa)) for aa in aminoAcids]
                    newRes = random.choices(aminoAcids, weights=aaWeights)[0]
                    sequenceComp[states[adj]] -= 1
                    sequenceComp[newRes] += 1
                    states[adj] = newRes

            if (pos2update not in fixedPos):
                # Choose adjaceny based on all edge dependent nodeValues
                aaWeights = [max(1e-5,posValue(pos2update, aa)) for aa in aminoAcids]
                newRes = random.choices(aminoAcids, weights=aaWeights)[0]
                sequenceComp[states[pos2update]] -= 1
                sequenceComp[newRes] += 1
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



def designPDB(pdbFile: os.path, chainSelection: list, sequenceFixMap: dict, outputLocation: os.path,
              sequenceLength: int, identifier:str = None, outSeqTables = outSeqTables, bg_AA_freq = bg_AA_freq):
    if(identifier == None):
        identifier = os.path.splitext(pdbFile)[0]
    parser = PDBParser()
    structure = parser.get_structure("", pdbFile)
    selectedChains = [chain for chain in structure.get_chains() if chain.get_id() in chainSelection]
    resConnections = getAdjMap(selectedChains)
    posProbMap, jointProbMap = initStruct(resConnections, outputLocation, parser)
    out, rosetta, ppm = buildMap(sequenceLength, posProbMap, jointProbMap, resConnections, sequenceFixMap)

    seqOut = os.path.join(outputLocation, identifier + '.fasta')
    with open(seqOut, 'w') as outFile:
        for entry in out:
            outFile.write('>' + str(round(entry[1], 2)) + '\n' + entry[0] + '\n')

    rosettaOut = os.path.join(outputLocation, identifier + '_resFile.txt')
    with open(rosettaOut, 'w') as outFile:
        for entry in out:
            outFile.write('>' + entry[0] + ', ' + str(round(entry[1], 2)) + '\n' + rosetta[entry] + '\n')

    ppm_obj = seqlogo.Ppm(ppm, alphabet_type='AA')
    seqLogoOut = os.path.join(outputLocation, identifier + '_seqLogo.png')
    seqlogo.seqlogo(ppm_obj, ic_scale=True, format='png', size='xlarge', filename=seqLogoOut)

    if(outSeqTables):
        jointProbOut = os.path.join(outputLocation, identifier + '_jointProb.txt')
        indepJointProbs = {}
        for X in bg_AA_freq:
            for Y in bg_AA_freq:
                indepJointProbs[(X, Y)] = 2 * bg_AA_freq[X] * bg_AA_freq[Y]
        positions = []
        for key in jointProbMap:
            positions.append(key)
        positions.sort()
        with open(jointProbOut, 'w') as outFile:
            for key in positions:
                outFile.write("Position: " + str(key) + '\n')
                values = []
                for resPair in indepJointProbs:
                    prob = jointProbMap[key][resPair]
                    logChange = math.log2(prob / indepJointProbs[resPair])
                    enrichment = prob * logChange
                    values.append((prob, enrichment, logChange, str(resPair)))
                values.sort(reverse=True)
                for val in values:
                    if (val[0] < 0.000001):
                        break
                    outFile.write(val[3] + ": " + str(round(val[0], 3)) + ', ' + str(round(val[1], 3)) + ', ' + str(
                        round(val[2], 3)) + '\n')



if __name__ == '__main__':
    import argparse

