import numpy as np
import sys, os
import warnings
from Bio.PDB.Polypeptide import *
from Bio.PDB import *

parser = PDBParser()
io = PDBIO()

# Generates hydrogen positions in a list, uses biophyton object as input. Based on fundamental bond length and angles
def genHydroFromBBAtoms(biopyChain):
    hydrogenCoords = {}
    for residue in biopyChain:
        ca = np.asarray(residue['CA'].get_coord())
        c = np.asarray(residue['C'].get_coord())
        n = np.asarray(residue['N'].get_coord())
        c_vect = c - ca
        c_vect = c_vect / np.linalg.norm(c_vect)
        n_vect = n - ca
        n_vect = n_vect / np.linalg.norm(n_vect)

        bisect_vect = -1 * np.array([(a + b) / 2 for a, b in zip(c_vect, n_vect)])
        bisect_vect = bisect_vect / np.linalg.norm(bisect_vect)
        norm_vect = np.cross(c_vect, n_vect)
        norm_vect = norm_vect / np.linalg.norm(norm_vect)

        # vector sum of normal vector, bisecting vector, and position vector, adjusted to length assuming 109.5 deg H-Ca-H bond and 1.09A C-H
        h1 = ca + norm_vect * (0.892876) + bisect_vect * (0.625198)
        h2 = ca - norm_vect * (0.892876) + bisect_vect * (0.625198)
        resNum = str(residue.get_id()[1])
        resName = residue.get_resname()
        hydrogenCoords[(resNum,resName)] = (h1,h2)
    return hydrogenCoords

def countHBond(pdbFile):
    hBondCount = 0
    hBondList = "Residue 1\tResidue 2\tXi angle\tCa-H . . 0 Distance\n"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        structure = parser.get_structure("", pdbFile)
        io.set_structure(structure)
        chainA, chainB = list(structure.get_chains())[:2]  # Get first two chains of BioPDB structure object

    chainA_hydros = genHydroFromBBAtoms(chainA)
    chainB_hydros = genHydroFromBBAtoms(chainB)

    # Using definition from https://www.pnas.org/doi/full/10.1073/pnas.161280798, The Ca-H...O hydrogen bond..., Senes, A.
    # of h-bond in experimental structure (H...O d < 3.5 && Ca-H..O ang > 120 || H...O d < 3 && Ca-H..O ang > 90)
    for residueA in chainA:
        ca_pos = np.asarray(residueA['CA'].get_coord())
        o_pos = np.asarray(residueA['O'].get_coord())
        for residueB in chainB_hydros:
            for hydro in chainB_hydros[residueB]:
                vect1 = ca_pos - hydro
                vect2 = o_pos - hydro
                dist = np.linalg.norm(vect2)
                ang = np.arcsin(np.dot(vect1,vect2)/(np.linalg.norm(vect1)*np.linalg.norm(vect2)))
                if(dist < 3.5 and ang > np.pi/2 or dist < 3.0 and ang > np.pi/4):
                    resANum = str(residueA.get_id()[1])
                    resAName = residueA.get_resname()
                    hBondList += "A "+resAName+' '+resANum+' \tB '+residueB[1]+' '+residueB[0]+' \t'+str(round(ang,2))+' \t'+str(round(dist,2))+'\n'
                    hBondCount += 1

    for residueB in chainB:
        ca_pos = np.asarray(residueB['CA'].get_coord())
        o_pos = np.asarray(residueB['O'].get_coord())
        for residueA in chainA_hydros:
            for hydro in chainA_hydros[residueA]:
                vect1 = ca_pos - hydro
                vect2 = o_pos - hydro
                dist = np.linalg.norm(vect2)
                ang = np.arcsin(np.dot(vect1, vect2) / (np.linalg.norm(vect1) * np.linalg.norm(vect2)))
                if (dist < 3.5 and ang > np.pi / 2 or dist < 3.0 and ang > np.pi / 4):
                    resBNum = str(residueB.get_id()[1])
                    resBName = residueB.get_resname()
                    hBondList += "A "+residueA[1]+' '+residueA[0]+' \tB '+resBName+' '+resBNum+' \t'+str(round(ang,2))+' \t'+str(round(dist,2))+'\n'
                    hBondCount += 1

    return hBondCount, hBondList


if __name__ == '__main__':
    # Run as python .\plotHBondHeat.py .\dimer_design_full\antipar_TMdimer_data.txt .\dimer_design_full\[FILE OUT NAME]

    # antiPar_TMdimer_Data
    params = sys.argv[1]
    id_to_param = {}
    with open(params,'r') as file:
        for line in file:
            line = line.strip().split()
            id_to_param[str(line[0])] =[int(line[2]), float(line[4]), float(line[5]), float(line[6]), int(line[7]),
                    float(line[8])]

    id_to_hb_count = {}
    for key in id_to_param.keys():
        identifier = key
        file_suffix = '.d29bf65c881e.allbb.pdb'
        dimer_db = os.path.join('Gx6G_dimers','AP_TM_dimer_database','AP_TM_dimer_database')
        struct = os.path.join(dimer_db, identifier + file_suffix)

        hBondCount, hBondList = countHBond(struct)
        with open(os.path.join('Gx6G_dimers','hbonds',identifier+'_hbond.txt'),'x') as outFile:
            outFile.write(hBondList)
        id_to_hb_count[key] = hBondCount
        if(int(key)%100 == 0):
            print(key)

    with open(sys.argv[2],'w') as outputFile:
        for key in id_to_param:
            outputFile.write(key+'\t'+str(id_to_hb_count[key])+'\n')

