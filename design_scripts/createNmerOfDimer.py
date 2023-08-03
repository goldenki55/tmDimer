import numpy as np
import os
from Bio.PDB import *
from collections import deque

# Returns a string that can be written to a file in PDB format. Inputs are for atom number, atom name (C/H/N...), residue
# name, chain, residue number, and the x,y,z coordinates of the atom.
def writeAtom(atNum, atName, resName, chain, resNum, x,y,z):
    atom = "ATOM   "
    atom += " "*(4-len(str(atNum)))
    atom += str(atNum)
    atom += "  "
    atom += atName
    atom += " "*(4-len(str(atName)))
    atom += resName
    atom += " "
    atom += chain
    atom += " "*(4-len(str(resNum)))
    atom += str(resNum)

    x = round(x,3)
    atom += " "*(12-len(str(x)))
    atom += str(x)
    y = round(y,3)
    atom += " "*(8-len(str(y)))
    atom += str(y)
    z = round(z,3)
    atom += " "*(8-len(str(z)))
    atom += str(z)
    atom += "  1.00               "
    atom += chain
    atom += "\n"
    return atom

def angBetweenVects(vect1, vect2):
    # Rounding removes floating point error for identical/nearly identical vectors
    angle = np.degrees(np.arccos(np.round(np.dot(vect1, vect2) / (np.linalg.norm(vect1) * np.linalg.norm(vect2)), 6)))
    return angle


# Takes inputs as numpy arrays
def projectPointOnLine(point, linePoint, lineVector):
    pointVect = point - linePoint
    projVect = lineVector * (np.dot(pointVect, lineVector) / np.dot(lineVector, lineVector))
    projPoint = linePoint + projVect
    return projPoint


def projectVectOnPlane(vector, planeNormal):
    vectOnNorm = planeNormal * (np.dot(vector, planeNormal) / np.dot(planeNormal, planeNormal))
    return vector - vectOnNorm

def quaternionRotatePointAboutAxis(point, axis, rotationRadians):
    # Adapted from http://answers.google.com/answers/threadview/id/361441.html, for 3D point and axis.
    axis = axis/np.linalg.norm(axis)
    cosR = np.cos(rotationRadians/2)
    sinR = np.sin(rotationRadians/2)
    q = {0: cosR, 1: sinR*axis[0], 2: sinR*axis[1], 3: sinR*axis[2]}
    Qmat = np.matrix([
        [(q[0]**2+q[1]**2-q[2]**2-q[3]**2), (2*(q[1]*q[2]-q[0]*q[3])), (2*(q[1]*q[3]+q[0]*q[2]))],
        [(2*(q[2]*q[1]+q[0]*q[3])), (q[0]**2-q[1]**2+q[2]**2-q[3]**2), (2*(q[2]*q[3]-q[0]*q[1]))],
        [(2*(q[3]*q[1]-q[0]*q[2])), (2*(q[3]*q[2]+q[0]*q[1])), (q[0]**2-q[1]**2-q[2]**2+q[3]**2)]
    ])
    point.shape = (3,1)  # Explicitly make point 2D column vector
    rotatedPoint = Qmat*point
    return rotatedPoint


def rotatePointAboutAffineVector(point, vector, vectorOrigin, rotationRadians):
    if(rotationRadians == 0):
        return point
    point = point - vectorOrigin
    point = quaternionRotatePointAboutAxis(point, vector, rotationRadians)
    point.shape = (3,)  # Explicitly make 2D column vector 1D array

    point = np.squeeze(np.asarray(point)) + vectorOrigin
    return point


# Loosely based off the psico implementation https://github.com/speleo3/pymol-psico/blob/master/psico/orientation.py#L10
def getHelixOrientation(residue):
    chain = residue.get_parent()
    centralResInd = int(residue.get_id()[1])
    idealEnd = centralResInd + 3
    residueList = deque()

    # Add up to 11 residues to the list such that the 7 vectors can be formed from each residue i's C to residue i+4's
    for resi in chain:
        residueList.append(resi)
        if (len(residueList) > 11):
            if (residueList[0] == residue):
                residueList.pop()
                break
            residueList.popleft()
        if (len(residueList) == 11 and int(residueList[-1].get_id()[1]) == idealEnd):
            break

    residueList = list(residueList)
    residueList = [resi for resi in residueList if abs(int(resi.get_id()[1]) - centralResInd) < 12]
    if (len(residueList) < 5):
        raise KeyError("Chain used to get helix orientation too small")
    orientationVectors = []
    oxyList = []
    for resLow, resHigh in zip(residueList[:-4], residueList[4:]):
        try:
            aC = np.array(resLow['C'].get_coord())
            aO = np.array(resLow['O'].get_coord())
            bN = np.array(resHigh['N'].get_coord())
        except KeyError:
            continue
        oxyList.append(aO)

        orientVect = (bN - aC)
        orientationVectors.append(orientVect)

    if (len(orientationVectors) == 0):
        raise KeyError("Insufficient C/O/N Containing Residues to form orientation vector")

    orientationVector = np.array(orientationVectors).sum(0)
    orientationVector = orientationVector / np.linalg.norm(orientationVector)
    center = np.array(oxyList).sum(0) / len(oxyList)
    center = center.astype('float')
    return center, orientationVector

# All input distances are in angstroms, all input angles are in degrees
# n the number of dimers to repeat
# alpha: [-180,180]degrees describes the angle between the centers of chainA-chainB-chainC, the technical range of alpha is +/-0 to (360-(360/n))
# newInterfaceDist: [0,inf)Angstrom describes the distance between chainB and chainC at the midpoint
# superHelRot: [-90,0]degrees describes the rotation of the superhelix about the midpoint
# dimerZRot: [-180,180]degrees describes the rotation about the vector spanning the two helix interface
# zOff: (-lenHelix,lenHelix)Angstrom describes the offset about the midpoint to do rotations about
# output is the filename
def createNmer(inputDimer, n, alpha, newInterfaceDist, dimerZRot, superHelRot, zOff, output):
    chainA, chainB = list(inputDimer.get_chains())[:2]
    chainAResList = list(chainA.get_residues())
    chainAMidRes = chainAResList[int(len(chainAResList)/2)]
    chainBResList = list(chainB.get_residues())
    chainBMidRes = chainBResList[int(len(chainBResList)/2)]

    centChainA, orientA = getHelixOrientation(chainAMidRes)
    centChainB, orientB = getHelixOrientation(chainBMidRes)
    resACa = np.asarray(chainAMidRes['CA'].get_coord())
    resBCa = np.asarray(chainBMidRes['CA'].get_coord())
    projRes1Ca = projectPointOnLine(resACa, centChainA, orientA)

    projRes2Ca = projectPointOnLine(resBCa, centChainB, orientB)
    dimerSpan = projRes2Ca - projRes1Ca
    unSignedPackAng = angBetweenVects(orientA, orientB)

    if(unSignedPackAng > 90):
        # Antiparallel interaction
        orientB = -orientB
    orientDimer = (1/2)*(orientA+orientB)
    orientDimer = orientDimer/np.linalg.norm(orientDimer)
    dimerCenter = (1/2)*(centChainA+centChainB)
    dimerNormal = np.cross(dimerSpan, orientDimer)

    beta = np.sign(alpha)*(360-(360/n)-abs(alpha))
    alpha = np.radians(alpha)
    beta = np.radians(beta)
    dimerRotationAmount = -1*np.sign(alpha)*np.radians(360/n)

    # Get centroid
    helixCenters = [centChainA, centChainB]
    prevA = centChainA
    prevB = centChainB
    # Implement offset of dimer center for this rotation along orientDimer
    dimerCenter = dimerCenter + zOff*orientDimer
    superHelVector = rotatePointAboutAffineVector(orientDimer, dimerNormal, dimerCenter, np.radians(superHelRot))
    for i in range(1, n):
        prevA = rotatePointAboutAffineVector(prevA, superHelVector, prevB, alpha)
        prevA = prevB + newInterfaceDist*((prevA - prevB)/np.linalg.norm(prevA - prevB))
        prevB = rotatePointAboutAffineVector(prevB, superHelVector, prevA, beta)
        prevB = prevA + newInterfaceDist*((prevB - prevA)/np.linalg.norm(prevB - prevA))
        helixCenters = helixCenters + [prevA, prevB]
    centroid = np.asarray(helixCenters).mean(axis=0)

    # Print output for display
    unitSuperHel = superHelVector/np.linalg.norm(superHelVector)
    #print('cgo_arrow '+str(list(centroid))+', '+str(list(centroid+6*unitSuperHel)))

    outPDBString = ""
    chainNames = [*"ABCDEFGHIJKLMNOPQRSTUVWXYZ"]
    atomNum = 1

    for i in range(n):
        for chainI, chain in enumerate([chainA, chainB]):
            for residue in chain:
                resName = residue.get_resname()
                resNum = residue.get_id()[1]
                for atom in residue:
                    coords = np.array(atom.get_coord())
                    coords = rotatePointAboutAffineVector(coords, dimerSpan, dimerCenter, dimerZRot)
                    coords = rotatePointAboutAffineVector(coords, superHelVector, centroid, i*dimerRotationAmount)
                    outPDBString += writeAtom(atomNum, atom.get_name(), resName, chainNames[2*i+chainI], resNum, coords[0], coords[1], coords[2])
                    atomNum += 1
            outPDBString += "TER\n"
    open(output,'w').write(outPDBString)


def createNmerDatabase(inputStructure, n, outputDirectory):
    interFaceMin = 8
    interfaceMax = 12
    interfaceStep = 0.5
    alphaMin = 80
    alphaMax = 160
    alphaStep = 5
    superHelRotMin = -45
    superHelRotMax = 0
    superHelRotStep = 5
    zOffMin = -1
    zOffMax = 1
    zOffStep = 1
    dimerZRotMin = 0
    dimerZRotMax = 3
    dimerZRotStep = 3
    listOfCreated = ""
    #listOfCreated = "ID,Hel Dist,Int Ang,SuperHel Tilt,Zoff,dimerZ Rot\n"
    interfaceDist = interFaceMin
    alpha = alphaMin
    superHelRot = superHelRotMin
    zOff = zOffMin
    dimerZRot = dimerZRotMin
    id = 9181
    # createNmer(inputDimer, n, alpha, newInterfaceDist, dimerZRot, superHelRot, zOff, output):
    while(interfaceDist <= interfaceMax):
        while(alpha <= alphaMax):
            while(superHelRot <= superHelRotMax):
                while(zOff <= zOffMax):
                    while(dimerZRot <= dimerZRotMax):
                        createNmer(inputStructure, n, alpha, interfaceDist, dimerZRot, superHelRot, zOff, os.path.join(outputDirectory,str(id)+'_nmer.pdb'))
                        listOfCreated += str(id)+','+str(interfaceDist)+','+str(alpha)+','+str(superHelRot)+','+str(zOff)+','+str(dimerZRot)+'\n'
                        id += 1
                        createNmer(inputStructure, n, -1*alpha, interfaceDist, dimerZRot, superHelRot, zOff, os.path.join(outputDirectory,str(id)+'_nmer.pdb'))
                        listOfCreated += str(id) + ',' + str(interfaceDist) + ',' + str(-1*alpha) + ',' + str(superHelRot) + ',' + str(zOff) + ',' + str(dimerZRot) + '\n'
                        id += 1
                        dimerZRot += dimerZRotStep
                    zOff += zOffStep
                    dimerZRot = dimerZRotMin
                zOff = zOffMin
                superHelRot += superHelRotStep
            superHelRot = superHelRotMin
            alpha += alphaStep
        alpha = alphaMin
        interfaceDist += interfaceStep
    return listOfCreated

if __name__ == '__main__':
    import sys

    PDBparser = PDBParser()
    structure = PDBparser.get_structure("", sys.argv[1])
    #createNmer(structure, 3, 120, 10, 0, -30, 0, os.path.join('trimer_of_jgd42','test5.pdb'))

    list = createNmerDatabase(structure, 3, sys.argv[2])
    with open(os.path.join(sys.argv[2],'database_labels.csv'),'a') as outfile:
        outfile.write(list)
    # https: // pymolwiki.org / index.php / Cgo_arrow
    # Use to visualize axis