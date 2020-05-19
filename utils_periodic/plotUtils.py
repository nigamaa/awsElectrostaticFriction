import numpy as np
import scipy as sc
import sys
import matplotlib.pyplot as plt
import pickle


def read_tip_force(rfname,sl,fields,sampleHeight):
    data = dict()

    for f in fields:
        data[f] = []

    rfile = open(rfname,'r')

    for l in rfile:
        tmp = l.split()
        tmp2 = list(map(float,tmp))

        for i in range(len(tmp2)):
            data[fields[i]].append(tmp2[i])


    data['fnorm'] = []
    for i in range(len(data['fx'])):
        data['zmin'][i] = data['zmin'][i] - sampleHeight
        tmp = np.sqrt(data['fx'][i]**2 + data['fy'][i]**2 + data['fz'][i]**2)
        data['fnorm'].append(tmp)

    return data



def plot_tip_force(fields,sampleHeight,rfname,figname,title,approach,retract,xlim=(-3,16),ylim=(-5,5)):
    data = read_tip_force(rfname,0,fields,sampleHeight)
    print('min value of force approach: ', min(data['fz'][0:approach]))
    print('min value of force retract:  ', min(data['fz'][approach:]))
    plt.figure(1,figsize=(15,5))
    plt.subplot(121)
    plt.plot(data['zmin'][0:approach],1*np.array(data['fz'][0:approach]),label='approach')
    if retract == True :
        plt.plot(data['zmin'][approach:],1*np.array(data['fz'][approach:]), label='retract')

    plt.plot(np.linspace(xlim[0],xlim[1],100), np.zeros((100,)),'k--')
    plt.xlim(xlim)
    plt.xlabel('Tip-Sample seperation (Angstroms)',fontsize=14)
    plt.ylabel('Force on the tip (nN)',fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.gca().invert_xaxis()
    plt.legend()
    plt.subplot(122)
    plt.plot(data['zmin'][0:approach],1*np.array(data['fz'][0:approach]),label='approach')
    if retract == True :
        plt.plot(data['zmin'][approach:],1*np.array(data['fz'][approach:]), label='retract')

    plt.plot(np.linspace(xlim[0],xlim[1],100), np.zeros((100,)),'k--')
    plt.ylim(ylim)
    plt.xlim(xlim)
    plt.xlabel('Tip-Sample seperation (Angstroms)',fontsize=14)
    plt.ylabel('Force on the tip (nN)',fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.gca().invert_xaxis()
    plt.legend()
    plt.suptitle(title)
    plt.savefig(figname)
    # plt.show()

    return
