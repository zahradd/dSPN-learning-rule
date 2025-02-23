from __future__ import print_function, division
import sys
import os
cur_dir = os.path.abspath(__file__)
model_dir = os.path.join(cur_dir, '..', '..',  '..','model')
updatedmodel_dir = os.path.abspath(model_dir)
sys.path.append(updatedmodel_dir)  
from neuron import h
import  MSN_builder as build
import pickle
import dFunc as func
# Load model mechanisms
import neuron as nrn
import random
import matplotlib.pyplot as plt
import numpy as np
import efel
import copy
import operator 
import spine as sp
current_dir = os.path.abspath(__file__)
new_dir = os.path.join(current_dir, '..', '..', '..', 'mechanisms')
mechanisms = os.path.abspath(new_dir)



h.nrn_load_dll(mechanisms + 'x86_64/.libs/libnrnmech.so')
nrn.load_mechanisms(mechanisms)
h.load_file('stdlib.hoc')
h.load_file('import3d.hoc')
h.load_file('stdrun.hoc')


class justSteep:
       def __init__(self,cell,rank,v):
          self.cell = cell
          self.rank= rank
          self.v=v
       def create_dopamine():
            h('dopa = 0')        
            dopa = h.dopa
            h('expectt = 0')        
            expectt= h.expectt
        
       def set_exglu_weights(synlist, netcons_s2t):
             #SYNLIST IS THE LIST OF NMDA SYNAPSES; EXNC IS THE NETCONS FROM THESE SYNAPSES TO INFFIRE1
             for s, nc in zip(synlist, netcons_s2t):
                 scale_factor=0.5
                 nc.weight[0] = ( s.weight )*scale_factor
                 
       def reset_exglu(thresholds):
            # for s in exglu:
            #     s.m = 0
            for i in thresholds.keys():
                thresholds[i].m = 0

       def meanFreq(vm,tm,treward,start,expect,noghte1,train1,stORsh): 
#        check for perdictet result
            if treward%2==0:
                trace1 = {}
                trace1['T'] = tm
                trace1['V'] =vm
                trace1['stim_start'] = [start-290]
                trace1['stim_end'] = [start-20]
                traces = [trace1]
                traces_results = efel.getFeatureValues(traces,
                                                       ['mean_frequency'])
                for trace_results in traces_results:
                # trace_result is a dictionary, with as keys the requested features
                  for feature_name, feature_values in trace_results.items():
                      if feature_name=='mean_frequency':
                         feature_valuesf=feature_values
                         if feature_values==None:
                             feature_valuesf=0
            #end of prediction 
                error=78*12 
                dop=justSteep.giveDopamine(feature_valuesf,treward,expect,stORsh[5],error)
                nog=justSteep.setPerform(feature_valuesf,treward,expect)
                # print('nog',nog,'dop',dop,'treward+2',treward/2)
                h.dopa=dop
                noghte1.append(nog)
                h.expectt=1
               
                print('input=(',train1[int(treward/2)],')','dopa=(',h.dopa,')  ','expect=(',expect[treward],')', 'mean frequency',feature_valuesf)                                               
            else:
                h.dopa=0
                h.expectt=0

              
                       
            return  h.dopa ,h.expectt
       def  giveDopamine(feature_valuesf,treward,expect,subOrnon,error):
           if subOrnon=='non':
                if  feature_valuesf>0:
                    if expect[treward]==1:
                       dop=1
                    else:
                       dop=-1
                else:
                   if expect[treward]==1:
                       dop=0
                   else:
                       dop=0
           else : 
               if treward/2>error:
                    if expect[treward]==1:
                        if treward/2>936:
                           dop=0
                        else:
                           dop=1
                    else:
                      dop=0
               else:  
                    if expect[treward]==1:
                       dop=1
                    else:
                      dop=-1
                   
                                       
           return dop
       def  setPerform(feature_valuesf,treward,expect):
                if  feature_valuesf>0:
                    if expect[treward]==1:
                       nog=1
                    else:
                       nog=0
                else:
                   if expect[treward]==1:
                       nog=0
                   else:
                       nog=1
                       
                return nog
       def callbb(self,expect,tstops,w_exe,w_in,dend_exe,dendinput,tresh,noghte1,train1,stORsh,num_syn,dendinput_in,num_syn_in):
 
        
            justSteep.create_dopamine()
         
            
    
            list_syn_ampa=[]
            list_nc_ampa=[]
            list_stimulator=[]
            list_syn_nmda=[]
            list_nc_nmda=[]
            list_syn_nmdaEx=[]
            list_nc_nmdaEx=[]
            netcons_s2t=[]

           
            list_syn_i=[]
            list_nc_i=[]
            list_stimulator_i=[] 
            spines = []
            # spines_cor=[]
            clusterOrDist=[0.1,0.9]
            # dend_exe=[3,4,5,12,14,15,21,22,24,26,27,28,29,35,36,37,40,41,45,46,47,48]
            # dend_exe=[3,12,15,21,26,35,41,46,53,57]
            # dend_exe=[3,5,9,12,14,21,24,26,35,41,45,46,53,57]
            dend_exe=[4,8,12,15,21,22,24,26,28,35,36,37,41,46,47,51,52,53,57,3,14,17,18,29,40,45,48,56,27,44]


            # dend_exe=[4,5,8,12,15,21,22,24,26,28,35,36,37,41,46,47,51,52,53,57,3,14,17,18,29,38,40,45,48,56,2,13,20,27,43,44,9,10]

            DendPosition1,DendPosition,secName,Output,sec_list=func.SetDendPosition(self.cell,dend_exe,num_syn,clusterOrDist,'random',self.rank)
          
            clusterOrDist_in=[0.1,0.9]
            dend_exe_in=[4,8,12,15,21,22,24,26,28,35,36,37,41,46,47,51,52,53,57,3,14,17,18,29,40,45,48,56,27,44]



            DendPosition1_in,DendPosition_in,secName_in,Output_in,sec_list_in=func.SetDendPosition(self.cell,dend_exe_in,num_syn_in,clusterOrDist_in,'random',self.rank)


            ListSynapse={}
            Listintfire={}
            for j in DendPosition.keys():
                Listintfire[j]=h.IntFire1()
                Listintfire[j].tau = 1000  # Time constant (ms)
                Listintfire[j].refrac = 50  # Refractory period (ms)
                my_list = []
                for  e in DendPosition[j]: 
                    i=e[0]
                    string = f'adaptive_shom_NMDAEX[{i}]'
                    my_list.append(string)
                
                ListSynapse[j]=my_list
            factor=0.5

            for i in range(num_syn):
                for j in DendPosition.keys():
                     for  e in DendPosition[j]: 
                         if i==e[0]:
                             spine_name = 'spine_' + j.name() + '(' + str(e[1]) + ')'
                             spines.append(sp.Spine(j, spine_name))
                             spines[-1].attach(j, e[1], 0)
                             sec = spines[-1].head
                             syn_ampa,syn_nmda,nc_ampa,nc_nmda,stimulator_e,syn_nmdaEx,nc_nmdaEx,netcon_s2t=func.putSynapse_task(self.cell,sec,j,e[1],dendinput[e[0]],w_exe[0][e[0]],tresh[0][e[0]],stORsh, ListSynapse ,Listintfire,factor)
                             h.setpointer(syn_nmda._ref_weight, 'weight', syn_nmdaEx)
                             h.setpointer(syn_nmda._ref_weight, 'weight', syn_ampa)
                             list_syn_ampa.append(syn_ampa)
                             list_nc_ampa.append(nc_ampa)
                             list_stimulator.append(stimulator_e)
                             list_syn_nmda.append(syn_nmda)
                             list_nc_nmda.append(nc_nmda)
                             list_syn_nmdaEx.append(syn_nmdaEx)
                             list_nc_nmdaEx.append(nc_nmdaEx)
                             netcons_s2t.append(netcon_s2t)
            
            
           
                     
            for i in range(num_syn_in):
                for j in DendPosition_in.keys():

                    for  e in DendPosition_in[j]: 
                        if i==e[0]:

                            syn_i,nc_i,stimulator_i=func.putSynapse_justSteepANn_inh(self.cell,j,e[1],dendinput_in[e[0]],w_in[0][e[0]],stORsh)
                            list_syn_i.append(syn_i)
                            list_nc_i.append(nc_i)
                            list_stimulator_i.append(stimulator_i)
       
    #         for  i2_e,sec_e in enumerate(seclist2):
    #              b=func.set_position(0.85,0.95,num_syn)   
    #              for iff_e in range(num_syn):
    #                  f =(i2_e*num_syn)+iff_e
    # #                 print('dend',dendinput2[f])
    #                  syn_i,nc_i,stimulator_i=func.putSynapse_in1(self.cell,sec_e,dendinput[f],w_in[0][f],stORsh,b[iff_e])
    #                  list_syn_i.append(syn_i)
    #                  list_nc_i.append(nc_i)
    #                  list_stimulator_i.append(stimulator_i)
        
  
            
            for syn in list_syn_nmda:
                h.setpointer(h._ref_dopa, 'dopa', syn)
           
            
            for syn in list_syn_i:
                h.setpointer(h._ref_expectt, 'expectt', syn)
   


                  
             
             
                ###############
            Synb, ncb, nsb = func.set_bg_noise( self.cell)

            tm = h.Vector()
            tm.record(h._ref_t)
            tmfff = h.Vector()
            tmfff.record(h._ref_t,100)
           
            vm = h.Vector()
            vm.record(self.cell.soma(0.5)._ref_v) 
            Rec_Wg_NMDA={}
            Rec_Wg_GABA={}
            checkD={}
            checkDnaro={}
            checkD_i={}
            checkD_imin={}

            checkDnaro_i={}
            checkDLtype={}
            for isyni,syni in enumerate(list_syn_i):
                Rec_Wg_GABA[isyni]= h.Vector()
                Rec_Wg_GABA[isyni].record(syni._ref_weight,100)
                checkD_i[isyni]= h.Vector()
                checkD_i[isyni].record(syni._ref_tresh,100)
                checkD_imin[isyni]= h.Vector()
                checkD_imin[isyni].record(syni._ref_tresh_min,100)
                checkDnaro_i[isyni]= h.Vector()
                checkDnaro_i[isyni].record(syni._ref_mltype,100)
            for isyn,syn1 in enumerate(list_syn_nmda):
                Rec_Wg_NMDA[isyn]= h.Vector()
                Rec_Wg_NMDA[isyn].record(syn1._ref_weight,100)
                checkD[isyn]= h.Vector()
                checkD[isyn].record(syn1._ref_tresh,100)
                checkDnaro[isyn]= h.Vector()
                checkDnaro[isyn].record(syn1._ref_conc0,100)
                checkDLtype[isyn]= h.Vector()
                checkDLtype[isyn].record(syn1._ref_mltype,100)
          
                
                
           # #########ezafe 
            # Rec_Vol_Dend={}
       
          
            # for  i2c,secc in enumerate(sec_list_in): 
            #       Rec_Vol_Dend[i2c]= h.Vector()
            #       Rec_Vol_Dend[i2c].record(secc(0.8)._ref_v) 
       
            ####################    
      
            h.finitialize(-80)
            
            for treward in range(len(tstops)-1):
               
                    h.CVode().event(tstops[treward], (justSteep.set_exglu_weights, (list_syn_nmda, netcons_s2t)))
                    h.CVode().event(tstops[treward], (justSteep.reset_exglu, (Listintfire)))
            
                    h.CVode().event(tstops[treward], (justSteep.meanFreq,(vm,tm,treward,tstops[treward],expect,noghte1,train1,stORsh)))     
            h.continuerun(tstops[-1])
            return vm,tmfff,tm,Rec_Wg_NMDA,Rec_Wg_GABA,secName, checkD,checkDnaro, checkD_i,checkD_imin,checkDnaro_i,checkDLtype ,Output,Output_in,DendPosition1
        
       def timeStep(start_time,num_task):
                startToLearn=start_time
                tstarts=[startToLearn]
                for i in range(0, num_task-1):
                    tstarts.append(startToLearn+(1+i)*800)
                ktstarts = [300] * len(tstarts)
                ktdz = [350] * len(tstarts)
                tstops = list(map(operator.add, tstarts, ktstarts)) 
                tstop1= list(map(operator.add, tstarts, ktdz)) 
                tstops =tstops +tstop1
                tstops.append(tstarts[-1]+800)  
                tstops.sort()  
                return tstarts,tstops
       def beforMain(self):
            # stORsh=['eta','alpha','mg','inh or not','example reward punish regular','general or non dopamine']
            stORsh=[0.381679389,0.062,1,0,'else','non']
         
            batch=12
            epoc=120
            num=batch*epoc
            noghte1=[]
            # np.random.seed(1)
            # w_exe_cor=abs(np.random.normal(0.2, 0.05, size=(1, corticalnoise)))
            tasks=['YB','RS','YS','RB']
            # tasks=['YB']
            # expectedOutput={'YB':1}

            expectedOutput={'YB':1,'RS':1,'YS':-1,'RB':-1}
            train1,expect=func.makeTrainingBatches(epoc,tasks,batch,expectedOutput,self.rank)
            # print('trian2',train1)
            tstarts,tstops=justSteep.timeStep(100,num)
            # print(tstarts,tstops)
            # dendlst=self.v
            synNumber1=200
            timespan=20
            NumberOfspike=1
         
               
            inputDict={'B':0.2,'S':0.2,'Y':0.2,'R':0.2,'E':0.2}
            assignInputTask={'YB':['Y','B','E'],'RS':['R','S','E'],'YS':['Y','S','E'],'RB':['R','B','E']}

            # assignInputTask={'YB':['Y','B'],'RS':['R','S'],'YS':['Y','S'],'RB':['R','B']}
            dendinput,synNumber1,inputListSyn=func.getSpiketrain(train1,tstarts,synNumber1,inputDict,assignInputTask,timespan,NumberOfspike,'random',self.rank)
            # print('dendinput',dendinput,inputListSyn)
            synNumber1_in=20
            timespan_in=100
            NumberOfspike_in=50

            np.random.seed(123)

            inputDict_in={'B':0.25,'S':0.25,'Y':0.25,'R':0.25}
            assignInputTask_in={'YB':['Y','B'],'RS':['R','S'],'YS':['Y','S'],'RB':['R','B']}
            dendinput_in,synNumber1_in,inputListSyn_in=func.getSpiketrain(train1,tstarts,synNumber1_in,inputDict_in,assignInputTask_in,timespan_in,NumberOfspike_in,'random',self.rank)

            tresh=abs(np.random.normal(0.02, 0.001, size=(1, synNumber1)))
           # var = self.rank * (0.8 - 0.001) / 100 + 0.001
            #print('varrrr',var)
            w_exe = abs(np.random.normal(0.3, 0.1, size=(1, synNumber1)))
            w_in=abs(np.random.uniform(0.2, 1, size=(1, synNumber1_in)))*1e-3

            dend_exe=[4,8,12,15,21,22,24,26,28,35,36,37,41,46,47,51,52,53,57,3,14,17,18,29,40,45,48,56,27,44]




           
     
            vm,tmfff,tm,Rec_Wg_NMDA,Rec_Wg_GABA,secName,checkD,checkDnaro, checkD_i,checkD_imin,checkDnaro_i,checkDLtype,Output ,Output_in,DendPosition1=justSteep.callbb(self,expect,tstops,w_exe,w_in,dend_exe,dendinput,tresh,noghte1,train1,stORsh,synNumber1,dendinput_in,synNumber1_in)
            # ,checkDendVoltage,checkDendNmda ,checkDendCal,checkD,checkDnaro,checkD_cor,checkDnaro_cor
 
          
             ##############################################################
            return   w_exe ,Rec_Wg_NMDA,w_in,Rec_Wg_GABA,noghte1,secName,tmfff,checkD,checkDnaro,train1, checkD_i,checkD_imin,checkDnaro_i,checkDLtype,inputListSyn,Output, inputListSyn_in,Output_in,DendPosition1
