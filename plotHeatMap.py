import numpy as np
import sys, os
import warnings
import matplotlib.pyplot as plt

params = sys.argv[1]
id_to_param = {}
r0 = []
w0 = []
alp = []
phi1 = []
zoff = []
with open(params,'r') as file:
    for line in file:
        line = line.strip().split()
        id_to_param[str(line[0])] =[float(line[4]), float(line[5]), float(line[6]), int(line[7]),
                float(line[8])]
        if(float(line[4]) not in r0):
            r0.append(float(line[4]))
        if(float(line[5]) not in w0):
            w0.append(float(line[5]))
        if(float(line[6]) not in alp):
            alp.append(float(line[6]))
        if(int(line[7]) not in phi1):
            phi1.append(int(line[7]))
        if(float(line[8]) not in zoff):
            zoff.append(float(line[8]))
r0.sort()
w0.sort()
alp.sort()
phi1.sort()
zoff.sort()

hBond = sys.argv[2]
id_to_hbond = {}
with open(hBond,'r') as file:
    for line in file:
        line = line.strip().split()
        id_to_hbond[str(line[0])] = int(line[1])


def plot(ax1, ax2):
    params = [r0, w0, alp, phi1, zoff]
    labels = ['Inter-helix distance r0 (Å)', 'Superhelical frequency, ωo', 'Pitch Angle α (degrees)', 'phase φ1 (degrees)','Z off (Å)']
    ax1_options = params[ax1]
    ax2_options = params[ax2]
    data = np.zeros((len(ax1_options),len(ax2_options)))
    counts = np.ones((len(ax1_options),len(ax2_options)))
    print(phi1)
    for id in id_to_param:
        vals = id_to_param[id]
        ind1 = ax1_options.index(vals[ax1])
        ind2 = ax2_options.index(vals[ax2])
        hBond = id_to_hbond[id]
        if(vals[0] == 3.8):
            data[ind1,ind2] += hBond
            counts[ind1,ind2] += 1

    data = data/counts

    # plt.imshow(data)
    # plt.colorbar()
    # plt.show()
    plt.rcParams.update({'font.size': 20})
    plt.rcParams.update({'font.sans-serif': 'Arial'})
    fig, ax = plt.subplots(1,1, figsize=(10,8))
    cp = ax.contourf(ax1_options,ax2_options,np.transpose(data),cmap=plt.get_cmap('YlGnBu'))
    fig.colorbar(cp)
    ax.set_xlabel(labels[ax1])
    ax.set_ylabel(labels[ax2])
    ax.set_title('Antiparallel Dimer H-bond, r0 = 3.8Å')
    plt.show()
plot(2,1)