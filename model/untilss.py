# -*- coding: utf-8 -*-
"""
Created on Sat Jul 15 13:13:14 2023

@author: zahra.khodadadi
"""

import os

import pickle
from neuron  import h
import numpy as np
import neuron as nrn
from scipy.signal import argrelextrema, find_peaks
#import seaborn as sns
import os
import pickle
import matplotlib.pyplot as plt
import numpy as np
import sys
def get_file_path(fileNameList, index):
    current_dir = os.path.abspath(__file__)
    new_dir = os.path.join(current_dir, '..', '..', '..', 'data', 'Plasticity.data')
    parent = os.path.abspath(new_dir)
    
    keys_list = list(fileNameList.keys())
    directory = keys_list[index[0]]
    fileName = fileNameList[directory][index[1]]
    path = os.path.join(parent, directory)
    
    return path,fileName


def save_obj(obj, name ):
    '''
    functions used to save data in the pickle format. 
    Used to pickle dictionaries
    
    obj     = dictionary with data
    name    = name to be used, without file ending (.pkl will be added)
    '''
    
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
        
        

def load_obj(name ):
    '''
    functions used to load data in the pickle format. 
    Used to un-pickle data from dictionaries
    
    name    = file name of the pickled data, including ending
    '''
    
    with open(name, 'rb') as f:
        return pickle.load(f) 

def calculate_all_distances():
    # this calculates all distances twice. With some if statements only once would have been enough 
    d = {}
    for s1 in h.allsec():   
        if s1.name()=='soma[0]':
            soma=s1
            d[s1.name()] = {}
        elif s1.name()=='axon[0]':
            continue
        else:           
            d['soma[0]'][s1.name()] = h.distance(soma(0.5), s1(0.5))

    return d['soma[0]']
def calculate_distances(pos,dend):
    # this calculates all distances twice. With some if statements only once would have been enough 
    distanceSoma = 0
    for s1 in h.allsec():   
        if s1.name()=='soma[0]':
            soma=s1
    for s1 in h.allsec():   
        if s1.name()==dend:
             distanceSoma=h.distance(soma(0.5), s1(0))
             distanceSoma=distanceSoma+pos*s1.L

    return distanceSoma
def process_dendHist(sec_name, Rec_inputListSyn, output):
    dendHist = {}
    list_2 = {}
    for indx, i in enumerate(sec_name):
        listt = {}
          # for nnn in output[0][indx]:
        for hhhj in Rec_inputListSyn:
          # print('nn', hhhj)
           l = []

           for nnn_in, nnn in enumerate(output[indx]):

                 if nnn in Rec_inputListSyn[hhhj]:
                        # print('true',nnn_in)
                         l.append(nnn)
           listt[hhhj]=l
        list_2[i] = listt
    dendHist= list_2
    return dendHist

def define_peaks(Rec_MaxCal_NMDA,Rec_time,rang) :   
    
    synNumber =len(Rec_MaxCal_NMDA[0])
    a = Rec_MaxCal_NMDA[0][0]

    pe = {}
    ge = {}
    tf = {}
    for i in range(rang):
        peakss = []
        gs = []
        tffs = []
        for j in range(synNumber):
            a1 = np.random.rand(1, a.size())
            for iff, value in enumerate(Rec_MaxCal_NMDA[i][j]):
                a1[0][iff] = value
            peaks, g = find_peaks(a1[0], height=0.000004)
            peakss.append(peaks)
            gs.append(g)
            tff = []
            for itf in peaks:
                tff.append(Rec_time[0][itf])
            tffs.append(tff)
        tf[i] = tffs
        pe[i] = peakss
        ge[i] = gs
    return  tf,pe,ge     


def set_dend_hist(dendHist, dendpositions):
    fcolor1 = {'R':'r', 'Y':'y', 'S':'k', 'B':'b','E':'m'}

    for i in dendHist.keys():
        for ifd in dendHist[i]:
            if  len(dendHist[i][ifd]) != 0:
                for pat in dendHist[i][ifd]:
                    for f in dendpositions[i]:
                        if f[0]==pat:
                            distanceSoma=calculate_distances(f[1],i)
                            f[1] = distanceSoma  # Update f[1] with new distanceSoma
    return    dendpositions                   




def extract_data(dendpositions,dendHist):
    data = dendpositions
    # Extract dendrite names and y values
    dend_names = []
    dend_values = []
    for dend, points in data.items():
        dend_names.append(dend)
        if points:
            _, y_values = zip(*points)
            dend_values.append(y_values)
        else:
            dend_values.append(0)
    inputsdend=[]
    for i,j  in zip(dendpositions.keys(),dendHist.keys()):
        inputsdend_ind=[]
        for ind in range(len(dendpositions[i])):
    
        #  if dendpositions[i][ind][0] in dendHist[j]['B']:
         #      print('we are hereb')
          #     inputsdend_ind.append('B')
          if dendpositions[i][ind][0] in dendHist[j]['S']:
               inputsdend_ind.append('S')
          if dendpositions[i][ind][0] in dendHist[j]['Y']:
               inputsdend_ind.append('Y')
      #    if dendpositions[i][ind][0] in dendHist[j]['R']:
       #        print('we are herer')
        #       inputsdend_ind.append('R')
         # if dendpositions[i][ind][0] in dendHist[j]['E']:
          #     print('we are heree')
           #    inputsdend_ind.append('E')
            
        inputsdend.append(inputsdend_ind) 
    return dend_names, dend_values,inputsdend

def plot_swarm_plot(dend_names, dend_values):
    plt.figure(figsize=(12, 6))
    sns.swarmplot(data=dend_values, size=6)
    for i in range(len(dend_values)):
        plt.vlines(i, min(dend_values[i]), max(dend_values[i]), colors='gray', linewidth=0.5)
    plt.xlabel('Dendrite')
    plt.ylabel('Y')
    plt.title('Swarm Plot - Clustering within Dendrites')
    plt.xticks(range(len(dend_names)), dend_names, rotation=45)
    plt.grid(True, axis='y', linestyle='--', linewidth=0.5)
    plt.tight_layout()
    plt.show()  
def plot_swarm_plotWithInputs(dend_names, dend_values ,inputsdend):
    colors=[]
    for i in inputsdend:
        color=[]
        for ind in i:
            color.append(ind)
        colors.append(color)   

    # Reshape the data into suitable format and assign random groups
    reshaped_data = {
        'dend': [],
        'Value': [],
        'Group': []
    }

    group_names = ['Y', 'B', 'S', 'R','E']
    num_groups = len(group_names)

    for i, group in enumerate(dend_values):
        for ind_val,value in enumerate(group):
          #  print(colors[i][ind_val])
            reshaped_data['dend'].append(f'd{i+1}')
            reshaped_data['Value'].append(value)
            reshaped_data['Group'].append(colors[i][ind_val])

    # Create swarm plot with hue
    sns.set(style="whitegrid")
    plt.figure(figsize=(20, 6))
    sns.swarmplot(x='dend', y='Value', hue='Group', data=reshaped_data, palette='Set1')

    # Customize the plot

    plt.legend(title='Group')
    for i in range(len(dend_values)):
        plt.vlines(i, min(dend_values[i]), max(dend_values[i]), colors='gray', linewidth=0.5)
    plt.xlabel('Dendrite')
    plt.ylabel('Y')
    plt.title('Swarm Plot - Clustering within Dendrites')
    plt.xticks(range(len(dend_names)), dend_names, rotation=45)
    plt.grid(True, axis='y', linestyle='--', linewidth=0.5)
    plt.tight_layout()
    plt.show()  

    # Show the plot
    plt.show()
