import os
from designDimer import *
import seqlogo
generalDir = 'dimer_design_full'
topScore = {}
with open(sys.argv[1], encoding='utf16') as input:
    for line in input:
        line = line.strip().split()
        # with open(line[0]) as pdb:
        #     print(pdb)
        #     print(int(line[1]))
        # break
        # Usage - python designDimer.py [pdb] [outputDir] [outputFileName w/ scores]
        identifier = line[0][line[0].find('0'):line[0].find('.d')]
        file_suffix = line[0][line[0].find('.d'):]
        dimer_db = os.path.join('Gx6G_dimers','AP_TM_dimer_database','AP_TM_dimer_database')
        struct = os.path.join(dimer_db,identifier+file_suffix)
        outPrefix = identifier
        dir = os.path.join(generalDir,identifier+'_dimer')
        if not os.path.isdir(dir):
            os.makedirs(dir)
        seqOut = os.path.join(dir, outPrefix + '.fasta')
        jointProbOut = os.path.join(dir, outPrefix + '_jointProb.txt')
        rosettaOut = os.path.join(dir, outPrefix + '_resFile.txt')
        seqLogoOut = os.path.join(dir, outPrefix + '_seqLogo.png')
        aminoAcids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'Y',
                      'W']

        TM_aaFrequency = {'A': 0.1102, 'R': 0.0256, 'N': 0.0225, 'D': 0.0135, 'C': 0.0141, 'Q': 0.0168, 'E': 0.0178,
                          'G': 0.0752, 'H': 0.0102, 'I': 0.1007, 'L': 0.1525, 'K': 0.0203, 'M': 0.0342, 'P': 0.0298,
                          'F': 0.0824, 'S': 0.0544, 'T': 0.0519, 'W': 0.0259, 'Y': 0.0409, 'V': 0.1005}
        # RUN
        structure = parser.get_structure("", struct)
        io.set_structure(structure)
        chainA, chainB = list(structure.get_chains())[:2]  # Get first two chains of BioPDB structure object

        # inputSeq = "XXXXXXXXXXXXXXXXXXXXXXXXXXXX"
        inputSeq = "XXXXAGXXAAXAGXXAAXAGXXAAXXXX"
        # inputSeq = "GGGGGGGGGGGGGGGGG"

        #fixedPos = {0: 'G', -7: 'G', 7: 'G'}
        fixedPos = {}
        centralG = int(line[1])

        try:
            resConnections = getAdjMap(struct, dir)
            posProbMap, jointProbMap, acenters, bcenters = initStruct(chainA, chainB, resConnections, centralG, dir)

            out, rosetta, ppm = buildMap(inputSeq, posProbMap, jointProbMap, resConnections, fixedPos, centralG)


            ppm_obj = seqlogo.Ppm(ppm, alphabet_type='AA')
            seqlogo.seqlogo(ppm_obj, ic_scale=True, format='png', size='xlarge', filename=seqLogoOut)
            topScore[identifier] = round(out[0][1], 2)
            with open(seqOut, 'w') as outFile:
                for entry in out:
                    outFile.write('>' + str(round(entry[1], 2)) + '\n' + entry[0] + '\n')
            with open(rosettaOut, 'w') as outFile:
                for entry in out:
                    outFile.write('>' + entry[0] + ', ' + str(round(entry[1], 2)) + '\n' + rosetta[entry] + '\n')

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
        except:
            pass



topScore = {k: v for k, v in sorted(topScore.items(), key=lambda item: item[1], reverse=True)}
with open(os.path.join(generalDir,'topScores.txt'),'w') as outFile:
    for entry in topScore:
        outFile.write(entry+'\t'+str(topScore[entry])+'\n')