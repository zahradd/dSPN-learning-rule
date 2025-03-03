import numpy
from neuron import h
import random
import pickle
import numpy as np
import itertools
import combinations_for_inputs as combi
from itertools import islice
from random import randint
from random import randrange

from collections import Counter

def randomListold(m, n,distribution,rank):
    random.seed(43+rank)

    # Create an array of size m where
    # every element is initialized to 0
    arr = [0] * m;
     
    # To make the sum of the final list as n
    for i in range(n) :
 
        # Increment any random element
        # from the array by 1
        if distribution=='random':
           arr[randint(0, n) % m] += 1;
        elif distribution=='uniform':
           arr[i % m] += 1;
  
    # Print the generated list
    return arr, m

def randomList(m, n, distribution, rank,length_list):
    np.random.seed(43 + rank)


    # Create an array of size m where every element is initialized to 0
    arr = [0] * m
    dend_list=list(range(m))
    # Calculate the total weight from the valuel list
    if distribution == 'random':
        # normalise length list (should sum to 1 in order to be able to be used with the numpy method)
        p = np.array(length_list) / sum(length_list)
        random_sections = np.random.choice(dend_list, size=n, p=p)
        # Count the occurrences of each dendrite
        dendrite_counts = Counter(random_sections)
        for dendrite, count in dendrite_counts.items():
          arr[dendrite] = count        
        
    elif distribution == 'uniform':
            # Distribute input uniformly
            for i in range(n):

                arr[i % m] += 1

    # Print the generated list
    return arr, m 


# Driver code
def multiples(value, length):
    return [value * i for i in range(1, length + 1)]

def getSpiketrain(train1,tstarts,synNumber1,inputDict,assignInputTask,timespan,NumberOfspike,distribution,rank):
        
    totalnumer=0
    inputprobebility=[]
    for i in inputDict:
        inputprobebility.append( int(inputDict[i]*synNumber1))
        totalnumer=int(inputDict[i]*synNumber1)+totalnumer
    Input =list(range(totalnumer))
    # Input1 =list(range(int(totalnumer/2)))
    # Input2 =list(range(int(totalnumer/2),totalnumer))

    if distribution=='random':
          random.seed(125+rank)

        
          random.shuffle(Input)


  

    
    Inputt = iter(Input)
    Output = [list(islice(Inputt, elem))
            for elem in inputprobebility]
    inputListSyn={}
    for ind,i in enumerate(Output):
        inputListSyn[list(inputDict)[ind]]=i
    assignTask={}
    for i in assignInputTask:
        l=[]
        for j in assignInputTask[i]:
             l.extend(inputListSyn[j])
        assignTask[i]=l
    dendsyntime={}

    for ist in range(totalnumer):
        dendsyntime[ist]=[]

    for i_iit,iit in  enumerate(train1): 
            # print(i_iit,iit)
            for i_inp in range(totalnumer):
                if i_inp in assignTask[iit]:
                  
                      spike1= random.sample(range(tstarts[i_iit], tstarts[i_iit]+timespan),NumberOfspike )
                      spike1.sort()
                      dendsyntime[i_inp].extend(spike1)
                else:
                    dendsyntime[i_inp].extend([])
    
    # print('outpkkkkut',dendsyntime)
  
    return dendsyntime,totalnumer,inputListSyn
    

    
          
def makeTrainingBatches(epoc, tasks, batch, expectedOutput,rank):
    expect = []
    train = []
    tr = tasks * (int(batch / len(tasks)))
    tr0 = tasks * (int(batch / len(tasks)))
    num = epoc * batch

    if num % batch == 0:

        
        tasks = tasks * int((batch / len(tasks)) * epoc)
        
        # Set a seed for reproducibility
        random.seed(42+rank)
        
        # Split the list into chunks of 12
        chunks = [tasks[i:i+batch] for i in range(0, len(tasks), batch)]
        
        # Shuffle only the chunks between the first and last one
        for i in range(1, len(chunks) - 1):
            random.shuffle(chunks[i])
        
        # Flatten the list of chunks to get the final list
        train = [item for chunk in chunks for item in chunk]
        
        print(train)

    else:
        print('the number of training should be multiple of len(task)')

    for i in train:
        expect.extend([expectedOutput[i], 0])
    return train, expect       

       

def makeRandomTrain(num):
    train=[]
    tr=['10','01','11','00','10','01','11','00','10','01','11','00']
    
    tr0=['10','01','11','00','10','01','11','00','10','01','11','00']

    # tr=['10','10','10','10','10','10','10','10','10','10','10','10']

    if num%12==0:
       n=int(num/12)
       for i in range(n):
           if i==0:
               train.extend(tr)
           elif i==79:
               
               train.extend(tr0)
           else:
                random.shuffle(tr)
                train.extend(tr)
    else:
        print('the number of training should be multiple of 12')
        
    return train
    
          
        
 
   
    
def makeSectionList(cell,dendl,clusterOrnot,num_syn):
    seclist=h.SectionList()
    secName=[] 
    pos={}
    for i ,sec in enumerate(dendl):
        dend_name   = 'dend[' + str(int(dendl[i])) + ']'
        for sec in cell.dendlist:
            if sec.name() == dend_name:
                seclist.append(sec=sec)
                secName.append(dend_name)
    ff=0
    for  i2_e,sec_e in enumerate(seclist):
        
        if clusterOrnot==1:
            sec_e.nseg=20
            pos[i2_e]=[]
            a2=np.random.uniform(0.85,0.95,size=num_syn)
            for x in a2:
                    pos[i2_e].append(x)

        else:
            pos[i2_e]=[]
            pos[i2_e].append([0.1,0.3,0.5,0.7,0.9,5,ff])
            ff=ff+5
          


    # print('pos',pos)
    return     seclist,secName,pos
def makeSectionList3(cell,dendl,clusterOrnot,num_syn,numberofsynapes):
    seclist=h.SectionList()
    secName=[] 
    pos={}
    for i ,sec in enumerate(dendl):
        dend_name   = 'dend[' + str(int(dendl[i])) + ']'
        for sec in cell.dendlist:
            if sec.name() == dend_name:
                seclist.append(sec=sec)
                secName.append(dend_name)
    ff=0
    for  i2_e,sec_e in enumerate(seclist):
        
        if clusterOrnot==1:
            sec_e.nseg=20
            pos[i2_e]=[]
            a2=np.random.uniform(0.85,0.95,size=num_syn)
            for x in a2:
                    pos[i2_e].append(x)

        else:
            pos[i2_e]=[]
            pos[i2_e].append([0.1,0.3,0.5,0.7,0.9,numberofsynapes,ff])
            ff=ff+numberofsynapes
          


    # print('pos',pos)
    return     seclist,secName,pos
def SetDendPosition(cell,dendl,num_syn,clusterOrDist,distribution,rank):
    seclist=h.SectionList()
    secName=[] 
    DendPosition={}
    DendPosition1={}
    for i ,sec in enumerate(dendl):
        dend_name   = 'dend[' + str(int(dendl[i])) + ']'
        for sec in cell.dendlist:
            if sec.name() == dend_name:
                seclist.append(sec=sec)
                secName.append(dend_name)
    # ff=0
    length_list=[149.95564708156397, 234.9497787721093, 212.87743881769075, 135.98179058422733, 109.16014441363774, 135.81257767517639, 122.34294330840919, 182.01664852383144, 102.44969358523703, 77.28090305713975, 53.51712542029221, 57.69529946180913, 119.0089598742341, 87.18195379742879,
 115.7911961383636, 57.819389668169215, 58.38736016719134, 186.15665551095813, 131.7771925283255, 148.21200520603827, 100.95230809946013, 78.42891647771215, 100.01810500760152, 11.473727688326658, 64.32599595106473, 47.12481434046398, 95.37066444687503, 81.01922773267117, 62.93034211825646,
 11.3560121815849]
    arr,m=randomList(len(dendl), num_syn,distribution,rank,length_list);
        
    #print('arr',arr)
    random.seed(546+rank)

    Input =list(range(num_syn))
    # random.shuffle(Input)
    
    Inputt = iter(Input)
    Output = [list(islice(Inputt, elem))
            for elem in arr]
    # random.randint(1, 9)*0.1
    d={}
    for  e,list_e in enumerate(Output):
        listma=[]
        for f in list_e:
            listma.append([f,random.uniform(clusterOrDist[0] * 10, clusterOrDist[1] * 10)*0.1])
        d[e]=listma    
        
    for  i2_e,sec_e in enumerate(seclist):
        sec_e.nseg=5
        listd=[]
        for i in  d[i2_e]:
            listd.append(i)
        DendPosition[sec_e]=listd
        DendPosition1[sec_e.name()]=listd


    return    DendPosition1,DendPosition,secName,Output,seclist
def SetDendPosition_differentdend(cell,dendl,num_syn,clusterOrDist,distribution,rank,length_list):
    seclist=h.SectionList()
    secName=[] 
    DendPosition={}
    DendPosition1={}
    for i ,sec in enumerate(dendl):
        dend_name   = 'dend[' + str(int(dendl[i])) + ']'
        for sec in cell.dendlist:
            if sec.name() == dend_name:
                seclist.append(sec=sec)
                secName.append(dend_name)
    # ff=0

    arr,m=randomList(len(dendl), num_syn,distribution,rank,length_list);
        
   # print('arr',arr)
    random.seed(546+rank)

    Input =list(range(num_syn))
    # random.shuffle(Input)
    
    Inputt = iter(Input)
    Output = [list(islice(Inputt, elem))
            for elem in arr]
    # random.randint(1, 9)*0.1
    d={}
    for  e,list_e in enumerate(Output):
        listma=[]
        for f in list_e:
            listma.append([f,random.uniform(clusterOrDist[0] * 10, clusterOrDist[1] * 10)*0.1])
        d[e]=listma    
        
    for  i2_e,sec_e in enumerate(seclist):
        sec_e.nseg=5
        listd=[]
        for i in  d[i2_e]:
            listd.append(i)
        DendPosition[sec_e]=listd
        DendPosition1[sec_e.name()]=listd


    return    DendPosition1,DendPosition,secName,Output,seclist
def set_bg_noise(cell,              \
                 ratio=0,    \
                 cell_type='D1',    \
                 syn_fact=False,    \
                 gabaMod=False,     \
                 delays=[]     ):
    
    ns      = {}
    nc      = {}
    Syn     = {}
    for s,sec in enumerate(cell.allseclist):
        
        # set bg noise----------------------------------
        
        if cell_type == 'D1':
            gbase = 0.3e-3 #previusly it was 0.3
        else:
            gbase = 0.2e-3
        
        if len(delays) == 0:
            delay = 0
        else:
            delay = delays[s]
            
        # create a glut synapse (glutamate)
        random_synapse(ns, nc, Syn, sec, 0.5,           \
                                NS_interval=80,    \
                                NC_conductance=gbase,       \
                                NS_start=delay )
        # create a gaba synapse (Exp2Syn)
        random_synapse(ns, nc, Syn, sec, 0.1,           \
                                Type='gaba',                \
                                NS_interval=50,       \
                                NC_conductance=1*1e-3 ,     \
                                NS_start=delay      )
        
        Syn[sec.name()+'_glut'].ratio = ratio  
        
        if syn_fact:
            Syn[sec.name()+'_glut'].ampa_scale_factor = syn_fact[0]
            Syn[sec.name()+'_glut'].nmda_scale_factor = syn_fact[1]
            
        
        if gabaMod:
            # scale gaba
            nc[sec.name()+'_gaba'].weight[0] = gbase * 3 * gabaMod
        
    
    return Syn, nc, ns
def set_bg_noise1(cell,              \
                 ratio=0,    \
                 cell_type='D1',    \
                 syn_fact=False,    \
                 gabaMod=False,     \
                 delays=[]     ):
    
    ns      = {}
    nc      = {}
    Syn     = {}
    for s,sec in enumerate(cell.allseclist):
        
        # set bg noise----------------------------------
        
        if cell_type == 'D1':
            gbase = 0.25e-3 #previusly it was 0.3
        else:
            gbase = 0.2e-3
        
        if len(delays) == 0:
            delay = 0
        else:
            delay = delays[s]
            
        # create a glut synapse (glutamate)
        random_synapse(ns, nc, Syn, sec, 0.5,           \
                                NS_interval=80,    \
                                NC_conductance=gbase,       \
                                NS_start=delay )
        # create a gaba synapse (Exp2Syn)
        random_synapse(ns, nc, Syn, sec, 0.1,           \
                                Type='gaba',                \
                                NS_interval=50,       \
                                NC_conductance=1*1e-3 ,     \
                                NS_start=delay      )
        
        Syn[sec.name()+'_glut'].ratio = ratio  
        
        if syn_fact:
            Syn[sec.name()+'_glut'].ampa_scale_factor = syn_fact[0]
            Syn[sec.name()+'_glut'].nmda_scale_factor = syn_fact[1]
            
        
        if gabaMod:
            # scale gaba
            nc[sec.name()+'_gaba'].weight[0] = gbase * 3 * gabaMod
        
    
    return Syn, nc, ns
def set_bg_noise_withconducancechanges(cell,  glut_conductance,gaba_conductance,ratio=0, cell_type='D1',syn_fact=False, gabaMod=False, delays=[]):
    
    ns      = {}
    nc      = {}
    Syn     = {}
    for s,sec in enumerate(cell.allseclist):
        
        # set bg noise----------------------------------
        
        if cell_type == 'D1':
            gbase = 0.3e-3 #previusly it was 0.3
        else:
            gbase = 0.2e-3
        
        if len(delays) == 0:
            delay = 0
        else:
            delay = delays[s]
            
        # create a glut synapse (glutamate)
        random_synapse(ns, nc, Syn, sec, 0.5,           \
                                NS_interval=80,    \
                                NC_conductance=glut_conductance,       \
                                NS_start=delay )
        # create a gaba synapse (Exp2Syn)
        random_synapse(ns, nc, Syn, sec, 0.1,           \
                                Type='gaba',                \
                                NS_interval=50,       \
                                NC_conductance=gaba_conductance ,     \
                                NS_start=delay      )
        
        Syn[sec.name()+'_glut'].ratio = ratio  
        
        if syn_fact:
            Syn[sec.name()+'_glut'].ampa_scale_factor = syn_fact[0]
            Syn[sec.name()+'_glut'].nmda_scale_factor = syn_fact[1]
            
        
        if gabaMod:
            # scale gaba
            nc[sec.name()+'_gaba'].weight[0] = gbase * 3 * gabaMod
        
    
    return Syn, nc, ns

def random_synapse(ns, nc, Syn, sec, x,         \
                Type='glut',                    \
                NS_start=0,                     \
                NS_interval=1000.0/18.0,        \
                NS_noise=1.0,                   \
                NS_number=100000,                 \
                S_AN_ratio=1.0,                 \
                S_tau_dep=100,                  \
                S_U=1,                          \
                S_e=-60,                        \
                S_tau1=0.25,                    \
                S_tau2=3.75,                    \
                NC_delay=0,                     \
                NC_conductance=0.6e-3,          \
                NC_threshold=0.1                ):
    '''
    random_synapse(argument, *, **)
    
    ---arg n removed. used for setting multiple gaba synapses in same segment
    
    creates a synapse in the segment closest to x in section sec, and updates the dict
    containing the synapses (as well as netStim and NetCon objects).
    
    Use the Type argument to specify synapse mechanism:
        Type        mechanism       description
        glut        tmglut          glutamatergic (ampa+nmda) with short term depression (default) 
        ampa        Exp2syn         NEURON native synapse with beta type dynamics
        gaba        Exp2syn         NEURON native synapse with beta type dynamics
        tmgabaa     tmgabaa         gabaergic with short term depression
    
    Any other Type than the above stated will result in an error.
        
    NS_arguments;       defines the NetStim object
    S_arguments;        defines the synapse mechanism
    NC_arguments;       defines the NetCon  object
    '''
    
    # create/set synapse in segment x of section
    if Type == 'tmglut':
        key                 = sec.name() + '_glut'
        Syn[key]            = h.tmGlut(x, sec=sec)
        Syn[key].nmda_ratio = S_AN_ratio
        Syn[key].tauR       = S_tau_dep
        Syn[key].U          = S_U
        
        
    elif Type == 'glut':
        key                 = sec.name() + '_glut'
        Syn[key]            = h.glutamate(x, sec=sec)
        Syn[key].ratio      = S_AN_ratio
        
    elif Type == 'gaba':
        key                 = sec.name() + '_gaba'  #+str(n)
        Syn[key]            = h.Exp2Syn(x, sec=sec)
        # what more?
        Syn[key].tau1       = S_tau1
        Syn[key].tau2       = S_tau2
        Syn[key].e          = -60
        
    elif Type == 'tmgabaa':
        
        key                 = sec.name() + '_gaba'
        Syn[key]            = h.tmGabaA(x, sec=sec)
        Syn[key].tauR       = S_tau_dep
        Syn[key].U          = S_U
        Syn[key].e          = S_e
        

        
    # create NetStim object
    ns[key]             = h.NetStim()
    ns[key].start       = NS_start
    ns[key].interval    = NS_interval # mean interval between two spikes in ms
    ns[key].noise       = NS_noise
    ns[key].number      = NS_number
    ns[key].seed( len(Syn) )
    
    # create NetCon object
    nc[key]             = h.NetCon(ns[key],Syn[key]) #  THIS IS WHERE THE ERROR WAS (Syn[sek] instead of Syn[key])
    nc[key].delay       = NC_delay
    nc[key].weight[0]   = NC_conductance
    nc[key].threshold   = NC_threshold
    

def putSynapse_justSteepANn(cell,spine_sec,secEx,pos,tra,w,tresh,stORsh,intfire1,intfire2,factor):
    stim  = h.VecStim()
    spike_times =tra
    spikes_vector = h.Vector(spike_times) # you put on vector just an spike time 
    stim.play(spikes_vector)
    synAmpa = h.adaptive_shom_AMPA(0.5,sec=spine_sec)
    synNmda = h.adaptive_shom_NMDA(0.5,sec=spine_sec)
    synNmdaEx = h.adaptive_shom_NMDAEX(pos,sec=secEx)
    synNmda.w0=w
    synNmda.treshf=tresh
    synNmda.eta=stORsh[0]
    synNmda.alpha=stORsh[1]
    synNmda.mg=stORsh[2]
    synNmda.rate_ltd=0.003
    synNmda.rate_ltp=0.000015
    synNmda.rate_ltd_thrsh=0.0000004
    synNmda.rate_ltp_tresh=0.0000001
    synNmdaEx.eta=stORsh[0]
    synNmdaEx.alpha=stORsh[1]
    synNmdaEx.mg=stORsh[2]
    synAmpa.gmax=1.5*1e-3
    synNmda.gmax=2.5*1e-3
    synNmdaEx.gmax=2.5*1e-3
    ncAmpa  = h.NetCon(stim , synAmpa)
    ncNmda  = h.NetCon(stim , synNmda)
    # ncNmdaEx  = h.NetCon(stim , synNmdaEx)
    # ncNmdaEx.delay=1.5
    ncAmpa.weight[0] = 1
    ncNmda.weight[0] = 1
    # ncNmdaEx.weight[0] = 1
    ncNmdaEx,netcon_s2t=generate_netcon(synNmdaEx,intfire1,intfire2,tra,w,stim,factor)


    return synAmpa,synNmda,ncAmpa,ncNmda,stim,synNmdaEx ,ncNmdaEx,netcon_s2t


def generate_netcon(synapse,intfire1,intfire2,t,w,stim,factor):
    name=str(synapse)
    my_list = []
    for i in range(20):
        string = f'adaptive_shom_NMDAEX[{i}]'
        my_list.append(string)
    if name in my_list:
        netcon = h.NetCon(intfire1, synapse)
        netcon_s2t=h.NetCon(stim, intfire1)

    else:
  
       netcon = h.NetCon(intfire2, synapse)
       netcon_s2t=h.NetCon(stim, intfire2)

    netcon.weight[0] = 1  # Weight of the connection (adjust as needed)
    netcon.delay=1.5
 
    netcon_s2t.weight[0] = w*factor
    netcon_s2t.delay = 0 
    return netcon,netcon_s2t
def putSynapse_task(cell,spine_sec,secEx,pos,tra,w,tresh,stORsh,ListSynapse ,Listintfire,factor):
    stim  = h.VecStim()
    spike_times =tra
    spikes_vector = h.Vector(spike_times) # you put on vector just an spike time 
    stim.play(spikes_vector)
    synAmpa = h.adaptive_shom_AMPA(0.5,sec=spine_sec)
    synNmda = h.adaptive_shom_NMDA(0.5,sec=spine_sec)
    synNmdaEx = h.adaptive_shom_NMDAEX(pos,sec=secEx)
    synNmda.w0=w
    synNmda.treshf=tresh
    synNmda.eta=stORsh[0]
    synNmda.alpha=stORsh[1]
    synNmda.mg=stORsh[2]
    synNmda.rate_ltd=0.003
    synNmda.rate_ltp=0.000015
    synNmda.rate_ltd_thrsh=0.0000004
    synNmda.rate_ltp_tresh=0.0000001
    synNmdaEx.eta=stORsh[0]
    synNmdaEx.alpha=stORsh[1]
    synNmdaEx.mg=stORsh[2]
    synAmpa.gmax=1.5*1e-3
    synNmda.gmax=2.5*1e-3
    synNmdaEx.gmax=2.5*1e-3
    ncAmpa  = h.NetCon(stim , synAmpa)
    ncNmda  = h.NetCon(stim , synNmda)
    # ncNmdaEx  = h.NetCon(stim , synNmdaEx)
    # ncNmdaEx.delay=1.5
    ncAmpa.weight[0] = 1
    ncNmda.weight[0] = 1
    # ncNmdaEx.weight[0] = 1
    ncNmdaEx,netcon_s2t=generate_netcon_task(synNmdaEx,tra,w,stim,ListSynapse ,Listintfire,factor)


    return synAmpa,synNmda,ncAmpa,ncNmda,stim,synNmdaEx ,ncNmdaEx,netcon_s2t

def putSynapse_task1(cell,spine_sec,secEx,pos,tra,w,tresh,stORsh):
    stim  = h.VecStim()
    spike_times =tra
    spikes_vector = h.Vector(spike_times) # you put on vector just an spike time 
    stim.play(spikes_vector)
    synAmpa = h.adaptive_shom_AMPA(0.5,sec=spine_sec)
    synNmda = h.adaptive_shom_NMDA(0.5,sec=spine_sec)
    synNmda.w0=w
    synNmda.treshf=tresh
    synNmda.eta=stORsh[0]
    synNmda.alpha=stORsh[1]
    synNmda.mg=stORsh[2]
    synNmda.rate_ltd=0.003
    synNmda.rate_ltp=0.000015
    synNmda.rate_ltd_thrsh=0.0000004
    synNmda.rate_ltp_tresh=0.0000001
    
    synAmpa.gmax=1.5*1e-3
    synNmda.gmax=2.5*1e-3
    ncAmpa  = h.NetCon(stim , synAmpa)
    ncNmda  = h.NetCon(stim , synNmda)
    # ncNmdaEx  = h.NetCon(stim , synNmdaEx)
    # ncNmdaEx.delay=1.5
    ncAmpa.weight[0] = 1
    ncNmda.weight[0] = 1
    # ncNmdaEx.weight[0] = 1
    return synAmpa,synNmda,ncAmpa,ncNmda,stim
def generate_netcon_task(synapse,t,w,stim,ListSynapse ,Listintfire,factor):
    name=str(synapse)
    for i in ListSynapse.keys():
        if name in ListSynapse[i]:
            netcon = h.NetCon(Listintfire[i], synapse)
            netcon_s2t=h.NetCon(stim, Listintfire[i])

    
    netcon.weight[0] = 1  # Weight of the connection (adjust as needed)
    netcon.delay=1.5
 
    netcon_s2t.weight[0] = w*factor
    netcon_s2t.delay = 0 
    return netcon,netcon_s2t
def putSynapse_justSteepANn_inh(cell,sec,pos,tra,w,stORsh):

    stimulator  = h.VecStim()
    spike_times =tra
    # print('tra',tra)
    syn = h.Gaba_mag(pos,sec=sec)
    syn.wmax=0.005
    syn.w0=w*stORsh[3]
    syn.e =-65
    spikes_vector = h.Vector(spike_times) # you put on vector just an spike time 
    stimulator.play(spikes_vector)       
    nc  = h.NetCon(stimulator , syn)
    nc.weight[0] = 1
    return syn,nc,stimulator
def putSynapse_justSteepANn_inhwithr(cell,sec,pos,tra,w,stORsh,wmax, tminrate , tmaxrate , wrate):

    stimulator  = h.VecStim()
    spike_times =tra
    # print('tra',tra)
    syn = h.Gaba_mag(pos,sec=sec)
    syn.wmax=wmax
    syn.w0=w*stORsh[3]
    syn.e =-65
    syn.tminrate=tminrate
    syn.tmaxrate=tmaxrate
    syn.wrate=wrate
    spikes_vector = h.Vector(spike_times) # you put on vector just an spike time 
    stimulator.play(spikes_vector)       
    nc  = h.NetCon(stimulator , syn)
    nc.weight[0] = 1
    return syn,nc,stimulator
def putSynapse_justSteepANncor(cell,spine_sec,secEx,pos,tra,w,tresh,stORsh):

    stim  = h.VecStim()
    spike_times =tra
    spikes_vector = h.Vector(spike_times) # you put on vector just an spike time 
    stim.play(spikes_vector)

    synAmpa = h.adaptive_shom_AMPA(0.5,sec=spine_sec)
    synNmda = h.adaptive_shom_NMDA(0.5,sec=spine_sec)
    synNmda.w0=w
    synNmda.treshf=tresh
    synNmda.eta=stORsh[0]
    synNmda.alpha=stORsh[1]
    synNmda.mg=stORsh[2]
    synNmda.rate_ltd=0.003
    synNmda.rate_ltp=0.000015
    synNmda.rate_ltd_thrsh=0.0000004
    synNmda.rate_ltp_tresh=0.0000001

    synAmpa.gmax=1.5*1e-3
    synNmda.gmax=2.5*1e-3
    ncAmpa  = h.NetCon(stim , synAmpa)
    ncNmda  = h.NetCon(stim , synNmda)
    ncAmpa.weight[0] = 1
    ncNmda.weight[0] = 1
    return synAmpa,synNmda,ncAmpa,ncNmda,stim

def putSynapse_justSteepANncorwratio(cell,spine_sec,secEx,pos,tra,w,tresh,stORsh,gmaxx,gmaxampa):

    stim  = h.VecStim()
    spike_times =tra
    spikes_vector = h.Vector(spike_times) # you put on vector just an spike time 
    stim.play(spikes_vector)

    synAmpa = h.adaptive_shom_AMPA(0.5,sec=spine_sec)
    synNmda = h.adaptive_shom_NMDA(0.5,sec=spine_sec)
    synNmda.w0=w
    synNmda.treshf=tresh
    synNmda.eta=stORsh[0]
    synNmda.alpha=stORsh[1]
    synNmda.mg=stORsh[2]
    synNmda.rate_ltd=0.003
    synNmda.rate_ltp=0.000015
    synNmda.rate_ltd_thrsh=0.0000004
    synNmda.rate_ltp_tresh=0.0000001

    synAmpa.gmax=gmaxampa
    synNmda.gmax=gmaxx
    ncAmpa  = h.NetCon(stim , synAmpa)
    ncNmda  = h.NetCon(stim , synNmda)
    ncAmpa.weight[0] = 1
    ncNmda.weight[0] = 1
    return synAmpa,synNmda,ncAmpa,ncNmda,stim


def putSynapse_justSteepANncorwratio1(cell,spine_sec,secEx,pos,tra,w,tresh,stORsh,gmaxx,gmaxampa):

    stim  = h.VecStim()
    spike_times =tra
    spikes_vector = h.Vector(spike_times) # you put on vector just an spike time 
    stim.play(spikes_vector)

    synAmpa = h.adaptive_shom_AMPA(0.5,sec=spine_sec)
    synNmda = h.adaptive_shom_NMDA(0.5,sec=spine_sec)
    synNmda.w0=w
    synNmda.treshf=tresh
    synNmda.eta=stORsh[0]
    synNmda.alpha=stORsh[1]
    synNmda.mg=stORsh[2]
    synNmda.rate_ltd=0.003
    synNmda.rate_ltp=0.000015
    synNmda.rate_ltd_thrsh=0.0000004
    synNmda.rate_ltp_tresh=0.0000001

    synAmpa.gmax=gmaxampa
    synNmda.gmax=gmaxx
    ncAmpa  = h.NetCon(stim , synAmpa)
    ncNmda  = h.NetCon(stim , synNmda)
    ncAmpa.weight[0] = 1
    ncNmda.weight[0] = 1

    
    return synAmpa,synNmda,ncAmpa,ncNmda,stim
def putSynapse_justSteepANn_acc(cell,spine_sec,secEx,pos,tra,w,tresh,stORsh):
    stim  = h.VecStim()
    spike_times =tra
    spikes_vector = h.Vector(spike_times) # you put on vector just an spike time 
    stim.play(spikes_vector)
    synAmpa = h.adaptive_shom_AMPA(0.5,sec=spine_sec)
    synNmda = h.adaptive_shom_NMDA(0.5,sec=spine_sec)
    synNmdaEx = h.adaptive_shom_NMDAEX(pos,sec=secEx)
    synNmda.w0=w
    synNmda.treshf=tresh
    synNmda.eta=stORsh[0]
    synNmda.alpha=stORsh[1]
    synNmda.mg=stORsh[2]
    synNmda.rate_ltd=0.003
    synNmda.rate_ltp=0.000015
    synNmda.rate_ltd_thrsh=0.0000004
    synNmda.rate_ltp_tresh=0.0000001
    synNmdaEx.eta=stORsh[0]
    synNmdaEx.alpha=stORsh[1]
    synNmdaEx.mg=stORsh[2]
    synAmpa.gmax=1.5*1e-3
    synNmda.gmax=2.5*1e-3
    synNmdaEx.gmax=5*1e-3
    ncAmpa  = h.NetCon(stim , synAmpa)
    ncNmda  = h.NetCon(stim , synNmda)
    ncNmdaEx  = h.NetCon(stim , synNmdaEx)
    ncNmdaEx.delay=5
    ncAmpa.weight[0] = 1
    ncNmda.weight[0] = 1
    ncNmdaEx.weight[0] = 1
    return synAmpa,synNmda,ncAmpa,ncNmda,stim,synNmdaEx ,ncNmdaEx

def putSynapse_withratio(cell,spine_sec,secEx,pos,tra,w,tresh,stORsh,ratio):

    stim  = h.VecStim()
    spike_times =tra
    spikes_vector = h.Vector(spike_times) # you put on vector just an spike time 
    stim.play(spikes_vector)

    synAmpa = h.adaptive_ratio_AMPA(0.5,sec=spine_sec)
    synNmda = h.adaptive_ratio_NMDA(0.5,sec=spine_sec)
    synNmda.w0=w
    synNmda.treshf=tresh
    synNmda.ratio_0=ratio
    synNmda.eta=stORsh[0]
    synNmda.alpha=stORsh[1]
    synNmda.mg=stORsh[2]
    synNmda.rate_ltd=0.003
    synNmda.rate_ltp=0.000015
    synNmda.rate_ltd_thrsh=0.0000004
    synNmda.rate_ltp_tresh=0.0000001

    synAmpa.gmax=1.5*1e-3
    synNmda.gmax=2.5*1e-3
    ncAmpa  = h.NetCon(stim , synAmpa)
    ncNmda  = h.NetCon(stim , synNmda)
    ncAmpa.weight[0] = 1
    ncNmda.weight[0] = 1


    
    return synAmpa,synNmda,ncAmpa,ncNmda,stim

def putSynapse_justSteep(cell,sec,pos,tra,w,tresh,stORsh):
    stim  = h.VecStim()
    spike_times =tra
    spikes_vector = h.Vector(spike_times) # you put on vector just an spike time 
    stim.play(spikes_vector) 
    syn = h.glutamate_steep(pos,sec=sec)
    syn.ratio=5
    syn.w0=w
    syn.treshf=tresh
    syn.wmax=1*1e-3
    syn.eta=stORsh[0]
    syn.alpha=stORsh[1]
    syn.mg=stORsh[2]
      
    nc  = h.NetCon(stim , syn)
    nc.weight[0] = 1
    return syn,nc,stim


import time

def custom_random():
    current_time = int(time.time() * 1000000)
    return current_time % 120

    



def putSynapse_in1(cell,sec,tra,w,stORsh,pos,magOrnot):
        #    difine syanpses
    tra2=[]   
    stimulator  = h.VecStim()
    #    tra = [x+10 for x in tra]
    #tra1 = [xtra for xtra in tra]
    if tra!=[]:
           for intra, value in enumerate(tra):
                    for i in range(0, 100,2):
                        tra2.append(value + i)

   
          # print('fffffffffffffff',tra2)
                   # tra2.append(value+280)
    spike_times =tra2

#    print('type tra',type(tra))
#    syn = h.glutamate_ica_nmda(0.95,sec)
    if magOrnot==1:
  #  syn = h.Gaba_p1(pos,sec=sec)
       syn = h.Gaba_mag(pos,sec=sec)
    else:
       syn = h.Gaba_p1(pos,sec=sec)

    syn.wmax=0.006
    syn.w0=w*stORsh[3]
    syn.e =-65
#    print('ratio',(3.5/w))
    spikes_vector = h.Vector(spike_times) # you put on vector just an spike time 
    stimulator.play(spikes_vector)       
    nc  = h.NetCon(stimulator , syn)
    nc.weight[0] = 1
    return syn,nc,stimulator



def putSynapse_in3(cell,sec,tra,w,stORsh,pos,magOrnot,wmax,interval):
        #    difine syanpses
    tra2=[]   
    stimulator  = h.VecStim()
    #    tra = [x+10 for x in tra]
    #tra1 = [xtra for xtra in tra]
    if tra!=[]:
           for intra, value in enumerate(tra):
                    for i in range(0, 100,interval):
                        tra2.append(value + i)

   
          # print('fffffffffffffff',tra2)
                   # tra2.append(value+280)
    spike_times =tra2

#    print('type tra',type(tra))
#    syn = h.glutamate_ica_nmda(0.95,sec)
    if magOrnot==1:
  #  syn = h.Gaba_p1(pos,sec=sec)
       syn = h.Gaba_mag(pos,sec=sec)
    else:
       syn = h.Gaba_p1(pos,sec=sec)

    syn.wmax=wmax
    syn.w0=w*stORsh[3]
    syn.e =-65
#    print('ratio',(3.5/w))
    spikes_vector = h.Vector(spike_times) # you put on vector just an spike time 
    stimulator.play(spikes_vector)       
    nc  = h.NetCon(stimulator , syn)
    nc.weight[0] = 1
    return syn,nc,stimulator

def putSynapse_in4(cell,sec,tra,w,stORsh,pos,magOrnot,wmax,interval, tminrate , tmaxrate , wrate):
        #    difine syanpses
    tra2=[]   
    stimulator  = h.VecStim()
    #    tra = [x+10 for x in tra]
    #tra1 = [xtra for xtra in tra]
    if tra!=[]:
           for intra, value in enumerate(tra):
                    for i in range(0, 100,interval):
                        tra2.append(value + i)

   
          # print('fffffffffffffff',tra2)
                   # tra2.append(value+280)
    spike_times =tra2


    syn = h.Gaba_mag(pos,sec=sec)
   

    syn.wmax=wmax
    syn.w0=w*stORsh[3]
    syn.e =-65
    syn.tminrate=tminrate
    syn.tmaxrate=tmaxrate
    syn.wrate=wrate
#    print('ratio',(3.5/w))
    spikes_vector = h.Vector(spike_times) # you put on vector just an spike time 
    stimulator.play(spikes_vector)       
    nc  = h.NetCon(stimulator , syn)
    nc.weight[0] = 1
    return syn,nc,stimulator
def putSynapse_in2(cell,sec,tra,w,stORsh,pos,magOrnot):
        #    difine syanpses
    stimulator  = h.VecStim()
#    tra = [x+10 for x in tra]
    tra1 = [xtra+random.randrange(150)-20 for xtra in tra]
    tra2=[]
    if tra1!=[]:
        for intra,itra in  enumerate(tra1):
            tra2.append(tra1[intra])
            tra2.append(tra1[intra]+5)
            tra2.append(tra1[intra]+10)
            tra2.append(tra1[intra]+15)
            tra2.append(tra1[intra]+20)
            tra2.append(tra1[intra]+25)
            tra2.append(tra1[intra]+30)
            tra2.append(tra1[intra]+35)
            tra2.append(tra1[intra]+40)
            tra2.append(tra1[intra]+45)
            tra2.append(tra1[intra]+50)

    spike_times =tra2
#    print('type tra',type(tra))
#    syn = h.glutamate_ica_nmda(0.95,sec)

    if magOrnot==1:
#  syn = h.Gaba_p1(pos,sec=sec)
         syn = h.Gaba_mag(pos,sec=sec)
    else:
         syn = h.Gaba_p1(pos,sec=sec)    
    syn.wmax=0.006
    syn.w0=w*stORsh[3]
    syn.e =-65
#    print('ratio',(3.5/w))
    spikes_vector = h.Vector(spike_times) # you put on vector just an spike time 
    stimulator.play(spikes_vector)       
    nc  = h.NetCon(stimulator , syn)
    nc.weight[0] = 1
    return syn,nc,stimulator

def finalHappyInput(flag,tstarts):
    if flag=='1':
      spike1=[]
      spike1.append(tstarts+random.randrange(21))
    else:
      spike1=[]
     
    return spike1
def getActiveInput(train):
    expec=[]
    activeInput=[]
    begu=[]
    for i in train:
        if i=='00':
            activeInput.append('1,3')
            expec.append(-1)
            begu.append(13)
        if i=='01':
            activeInput.append('1,4')
            expec.append(1)
            begu.append(14)
            
        if i=='10':
            activeInput.append('2,3')
            expec.append(1)
            begu.append(23)
        if i=='11':
            activeInput.append('2,4')
            expec.append(-1)
            begu.append(24)
            
    kexpect = [0] * len(expec)
    expect=[]
    for iexpe in range(len(expec)):
            expect.append(expec[iexpe])
            expect.append(kexpect[iexpe])       
    return begu,activeInput,expect

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

         

       
        

def taghsim(ich,num_syn):
    jcheck=ich%num_syn
    icheck=ich//num_syn
    return icheck,jcheck

def getDendinput(train1,tstarts,synNumber,num_syn):
    dendinput={}
    num_syn1=num_syn//4
    list1=list(range(num_syn1))
    list1=list1+[x+num_syn for x in list1]
    list2=list(range(num_syn1,num_syn1*2))
    list2=list2+[x+num_syn for x in list2]
    list3=list(range(num_syn1*2,num_syn1*3))
    list3=list3+[x+num_syn for x in list3]
    list4=list(range(num_syn1*3,num_syn1*4))
    list4=list4+[x+num_syn for x in list4]
  #  print('lissssssssssssssssssst','list1',list1,'list2',list2,'list3',list3,'list4',list4)
    for idendinput in range(synNumber):
        dendinput[idendinput]=[]
    st={}
    for ist in range(synNumber):
        st[ist]=[]
    
    for i_iit,iit in  enumerate(train1): 
        # inp=[]  
        if iit=='00':
            # inp=np.loadtxt('00.txt') 
            for i_inp in range(synNumber):
                if i_inp in list1+list3:
                   st[i_inp].append(finalHappyInput('1',tstarts[i_iit]))
                else:
                   st[i_inp].append(finalHappyInput('0',tstarts[i_iit]))
                

        if iit=='01':
            # inp=np.loadtxt('01.txt') 
            for i_inp in range(synNumber):
               if i_inp in list1+list4:
                   st[i_inp].append(finalHappyInput('1',tstarts[i_iit])) 
               else:
                   st[i_inp].append(finalHappyInput('0',tstarts[i_iit]))

        if iit=='10':
            # inp=np.loadtxt('10.txt') 
            for i_inp in range(synNumber):
                if i_inp in list2+list3:
                   st[i_inp].append(finalHappyInput('1',tstarts[i_iit])) 
                else:
                   st[i_inp].append(finalHappyInput('0',tstarts[i_iit]))
   
        if iit=='11':
            # inp=np.loadtxt('11.txt') 
            for i_inp in range(synNumber):
                if i_inp in list2+list4:
                   st[i_inp].append(finalHappyInput('1',tstarts[i_iit]))
                else:
                   st[i_inp].append(finalHappyInput('0',tstarts[i_iit]))
        if iit=='noise':
            for i_inp in range(synNumber):
                 st[i_inp].append(finalHappyInput_corticalnoise(tstarts[i_iit]))     



    for iFlat in range(synNumber):
         st[iFlat]= list(itertools.chain(*st[iFlat]))
               
    
#    print('st',st)    #................. 
    for idi in range(synNumber):
        dendinput[idi].append(st[idi])
    for dFlat in range(synNumber):
         dendinput[dFlat]= list(itertools.chain(*dendinput[dFlat])) 
         
    return dendinput       


def finalHappyInput_corticalnoise(tstarts):
      spike1=[]
      spike1.append(tstarts+random.randrange(51))
       
      return spike1
def makew_inex1(synNumber,ra,spillornot,num_syn,taskjob):
# w_in=np.ones((1,synNumber))*0*1e-3
       # w_exe=np.ones((1,synNumber))*0
       w_exe=np.zeros((1,synNumber))
       # w_exe=abs(np.random.normal(0.1, 0.0001, size=(1, synNumber))*1e-3)
       w_in=np.random.rand(1,synNumber)*0.01*1e-3
      
       ranra=[[[1,2,3],[1,2,4]]]*32
       inputs_linear=[]
       inputs_XOR=[]
       inputs_none=[]
       for i in combi.indices_linear:
            inputs_linear.append(combi.all_combinations[i])
       for i in combi.indices_XOR:
            inputs_XOR.append(combi.all_combinations[i])
       for i in combi.indices_none:
            inputs_none.append(combi.all_combinations[i])  
       if taskjob=='combination'   :  
       # input_combination=inputs_XOR
           input_combination=inputs_XOR
       elif taskjob=='linear':
           input_combination=inputs_linear
       elif taskjob=='non':
               input_combination=inputs_none
       elif taskjob=='3example':
                  input_combination=ranra
       #input_combination=ranra

       num_syn1=num_syn//4
       list1=[]
       list2=[]
       l1=list(range(num_syn1))
       l2=list(range(num_syn1,num_syn1*2))
       l3=list(range(num_syn1*2,num_syn1*3))
       l4=list(range(num_syn1*3,num_syn1*4))
       list3=list(range(num_syn))


       for i in input_combination[ra][0]:
           
           if i==1:
               for il in l1:
                  list1.append(il)
           if i==2:
               for il in l2:
                 list1.append(il)
           if i==3:
               for il in l3:
                 list1.append(il)
           if i==4:
               for il in l4:
                list1.append(il)
                
       for i in input_combination[ra][1]:
           if i==1:
               for il in l1:
                  list2.append(il)
           if i==2:
               for il in l2:
                 list2.append(il)
           if i==3:
               for il in l3:
                 list2.append(il)
           if i==4:
               for il in l4:
                list2.append(il)            
       # print(list1)
       listbenana=[]
       listStraberry=[]
       listbenana_in=[]
       listStraberry_in=[]
       for ib in [0]:
            for ittt in list1:
                #[5,6,7,8,9,10,11,12,13,14]
                listbenana.append(ib*num_syn+ittt)
       for ib in [1]:
            for ittt in list2:
                #[0,1,2,3,4,15,16,17,18,19]
                listStraberry.append(ib*num_syn+ittt)     
       for ib in [1]:
             for ittt in list3:
                 #[5,6,7,8,9,10,11,12,13,14]
                 listbenana_in.append(ib*num_syn+ittt)
       for ib in [0]:
             for ittt in list3:
                 #[0,1,2,3,4,15,16,17,18,19]
                 listStraberry_in.append(ib*num_syn+ittt)          
               
       for indexw in listbenana+listStraberry:
           np.random.seed(indexw+123456)
           if spillornot=='spill':
              # w_pla=np.random.rand(1,1)*0.6
               # w_pla= abs(np.random.normal(0.25, 0.01, size=(1, 1)))
               w_pla= abs(np.random.normal(0.25, 0.05, size=(1, 1)))
           elif spillornot=='spill_inh':
                  # w_pla=np.random.rand(1,1)*0.6
                   # w_pla= abs(np.random.normal(0.25, 0.01, size=(1, 1)))
                   w_pla= abs(np.random.normal(0.45, 0.05, size=(1, 1)))

           elif spillornot=='reward':
               w_pla= abs(np.random.normal(0.19, 0.1, size=(1, 1)))
           elif spillornot=='punish':
               w_pla= abs(np.random.normal(0.42, 0, size=(1, 1)))
           elif spillornot=='3_example':
               w_pla= abs(np.random.normal(0.25, 0, size=(1, 1)))

           # w_pla=np.ones((1,1))*1
           w_exe[0][indexw]=w_exe[0][indexw]+w_pla[0]
      # for indexw in listbenana_in+listStraberry_in:
       #    if spillornot=='3_example':
            # w_pla=np.random.rand(1,1)*0.5*1e-3
        #    w_pla=np.ones((1,1))*3*1e-3
         #   w_in[0][indexw]=w_in[0][indexw]+w_pla[0]
          # else:
       w_in=abs(np.random.normal(0.1, 0.001, size=(1, len(listbenana_in+listStraberry_in))))*1e-3
           # w_in[0][indexw]=w_in[0][indexw]+w_pla[0]
       return w_exe,w_in   
def set_position(minn,maxx,synNumber):
     b=[]
     a=np.random.uniform(minn,maxx,size=synNumber)
     for x in a:
        b.append(int(round(x*100))/100)
     return b
   
    