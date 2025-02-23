# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 23:40:02 2019

@author: zahra.khodadadi
"""
###################
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

####################

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
                 
       def reset_exglu(threshold,threshold2):
            # for s in exglu:
            #     s.m = 0 
           threshold.m = 0
           threshold2.m = 0


               
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
                error=35*12 
                dop=justSteep.giveDopamine(feature_valuesf,treward,expect,stORsh[5],error)
              
                nog=justSteep.setPerform(feature_valuesf,treward,expect)
                # print('nog',nog,'dop',dop,'treward+2',treward/2)
                h.dopa=dop
                noghte1.append(nog)
                h.expectt=1
               
                print('input=(',train1[int(treward/2)],')','dopa=(',h.dopa,') , ','expect=(',expect[treward],')', 'mean frequency',feature_valuesf)                                               
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
                        if treward/2>480:
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
       def callbb(self,expect,tstops,tstarts,w_exe,w_in,dend_exe,dendinput,dendinput_inh,tresh,noghte1,train1,stORsh,dendlst,w_exe_cor,dendinput_cor,tresh_cor,num_syn,num_syn_inh):
 
        
            justSteep.create_dopamine()
         
            
            seclist2,secName,pos=func.makeSectionList(self.cell,dendlst,1,num_syn)
            seclist2_inh,secName_inh,pos_inh=func.makeSectionList(self.cell,dendlst,1,num_syn_inh)

            dend_exe3=list(set(dend_exe) - set(dendlst)) + list(set(dendlst) - set(dend_exe))
            seclist20_cor,secName_cor,pos_cor=func.makeSectionList3(self.cell,dend_exe,0,num_syn,4)  
            list_syn_ampa=[]
            list_nc_ampa=[]
            list_stimulator=[]
            list_syn_nmda=[]
            list_nc_nmda=[]
            list_syn_nmdaEx=[]
            list_nc_nmdaEx=[]
            list_syn_e_corAmpa=[]
            list_nc_e_corAmpa=[]
            list_stimulator_e_cor=[]
            list_syn_e_corNmda=[]
            list_nc_e_corNmda=[]
            netcons_s2t=[]
            list_syn_i=[]
            list_nc_i=[]
            list_stimulator_i=[] 
            spines = []
            spines_cor=[]
            synlist=[]
           ############################################ # Glutamate thresholdjjjj
            intfire1 = h.IntFire1()
            intfire1.tau = 1000  # Time constant (ms)
            intfire1.refrac = 50  # Refractory period (ms)
            intfire2 = h.IntFire1()
            intfire2.tau = 1000  # Time constant (ms)
            intfire2.refrac = 50  # Refractory period (ms)
            factor=0.5
            ################################################################  jjjj
            
            for  i2_e,sec_e in enumerate(seclist2): 
                 
                 for iff_e in range(num_syn):
                     f =(i2_e*num_syn)+iff_e

                     spine_name = 'spine_' + sec_e.name() + '(' + str(pos[i2_e][iff_e]) + ')'
                     spines.append(sp.Spine(sec_e, spine_name))
                     spines[-1].attach(sec_e, pos[i2_e][iff_e], 0)
                     sec = spines[-1].head
                     syn_ampa,syn_nmda,nc_ampa,nc_nmda,stimulator_e,syn_nmdaEx,nc_nmdaEx,netcon_s2t=func.putSynapse_justSteepANn(self.cell,sec,sec_e,pos[i2_e][iff_e],dendinput[f],w_exe[0][f],tresh[0][f],stORsh,intfire1,intfire2,factor)

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

            synlist=list_syn_nmda  
                     
            for  i2_e,sec_e in enumerate(seclist20_cor): 
                 for iff_e in range(pos_cor[i2_e][0][-2]):
                      f =pos_cor[i2_e][0][-1]+iff_e
                      spine_name_cor = 'spine_' + sec_e.name() + '(' + str(pos_cor[i2_e][0][iff_e]) + ')'
                      spines_cor.append(sp.Spine(sec_e, spine_name_cor))
                      spines_cor[-1].attach(sec_e,pos_cor[i2_e][0][iff_e], 0)
                      sec = spines_cor[-1].head
                      syn_ampa_cor,syn_nmda_cor,nc_ampa_cor,nc_nmda_cor,stimulator_e_cor=func.putSynapse_justSteepANncor(self.cell,sec,sec_e,pos_cor[i2_e][0][iff_e],dendinput_cor[f],w_exe_cor[0][f],tresh_cor[0][f],stORsh)
                      h.setpointer(syn_nmda_cor._ref_weight, 'weight', syn_ampa_cor)
                      list_syn_e_corAmpa.append(syn_ampa_cor)
                      list_nc_e_corAmpa.append(nc_ampa_cor)
                      list_stimulator_e_cor.append(stimulator_e_cor)    
                      list_syn_e_corNmda.append(syn_nmda_cor)
                      list_nc_e_corNmda.append(nc_nmda_cor)
            wmax=0.005
            interval= 2 
            tminrate=0.00005
            tmaxrate=0.0005
            wrate=0.055
            for  i2_e,sec_e in enumerate(seclist2_inh):
                 b=func.set_position(0.55,0.6,num_syn_inh)   
                 for iff_e in range(num_syn_inh):
                     f =(i2_e*num_syn_inh)+iff_e
                     syn_i,nc_i,stimulator_i=func.putSynapse_in4(self.cell,sec_e,dendinput_inh[f],w_in[0][f],stORsh,b[iff_e],0,wmax,interval, tminrate , tmaxrate , wrate)
                     list_syn_i.append(syn_i)
                     list_nc_i.append(nc_i)
                     list_stimulator_i.append(stimulator_i)
  
            
            for syn in list_syn_nmda:
                h.setpointer(h._ref_dopa, 'dopa', syn)
           
            for syn in list_syn_e_corNmda:
                h.setpointer(h._ref_dopa, 'dopa', syn)
            
            for syn in list_syn_i:
                h.setpointer(h._ref_expectt, 'expectt', syn)

                  
             
             
                ###############
            Synb, ncb, nsb = func.set_bg_noise( self.cell)
            
            # m = h.Vector()
            # m.record(intfire1._ref_m)
            # m1 = h.Vector()
            # m1.record(intfire2._ref_m)
            tm = h.Vector()
            tm.record(h._ref_t)
            tmfff = h.Vector()
            tmfff.record(h._ref_t,100)
           
            vm = h.Vector()
            vm.record(self.cell.soma(0.5)._ref_v) 
            Rec_Wg_NMDA={}
            Rec_Wg_GABA={}
            Rec_Wg_NMDA_extraInput={}
            checkD={}
            checkDnaro={}
            checkD_cor={}
            checkDnaro_cor={}
            checkD_i={}
            checkDnaro_i={}
            checkDLtype_cor={}
            checkDLtype={}
            checkD_imin={}
            for isyni,syni in enumerate(list_syn_i):
                Rec_Wg_GABA[isyni]= h.Vector()
                Rec_Wg_GABA[isyni].record(syni._ref_weight,100)
                checkD_i[isyni]= h.Vector()
                checkD_i[isyni].record(syni._ref_tresh,100)
                checkDnaro_i[isyni]= h.Vector()
                checkDnaro_i[isyni].record(syni._ref_mltype,100)
                checkD_imin[isyni]= h.Vector()
                checkD_imin[isyni].record(syni._ref_tresh_min,100)
                
            for isyn,syn1 in enumerate(list_syn_nmda):
                Rec_Wg_NMDA[isyn]= h.Vector()
                Rec_Wg_NMDA[isyn].record(syn1._ref_weight,100)
                checkD[isyn]= h.Vector()
                checkD[isyn].record(syn1._ref_tresh,100)
                checkDnaro[isyn]= h.Vector()
                checkDnaro[isyn].record(syn1._ref_conc0,100)
                checkDLtype[isyn]= h.Vector()
                checkDLtype[isyn].record(syn1._ref_mltype,100)
                
                
            # for isyn,syn2 in enumerate(list_syn_nmdaEx):
            #        Rec_Wg_NMDA[isyn]= h.Vector()
            #        Rec_Wg_NMDA[isyn].record(syn2._ref_gluextra)
                   
                   
                   
                   
            for isyn,syn1 in enumerate(list_syn_e_corNmda):
                Rec_Wg_NMDA_extraInput[isyn]= h.Vector()
                Rec_Wg_NMDA_extraInput[isyn].record(syn1._ref_weight,100)
                checkD_cor[isyn]= h.Vector()
                checkD_cor[isyn].record(syn1._ref_tresh,100)
                checkDnaro_cor[isyn]= h.Vector()
                checkDnaro_cor[isyn].record(syn1._ref_conc0,100)
                checkDLtype_cor[isyn]= h.Vector()
                checkDLtype_cor[isyn].record(syn1._ref_mltype,100)
                
                
           # #########ezafe 
            # Rec_Vol_Dend={}
       
          
            # for  i2c,secc in enumerate(seclist2): 
            #       Rec_Vol_Dend[i2c]= h.Vector()
            #       Rec_Vol_Dend[i2c].record(secc(0.9)._ref_v) 
       
            ####################    
      
            h.finitialize(-80)
            # fih = h.FInitializeHandler((justSteep.seti_events,( tstarts,list_syn_nmdaEx , netcons_s2t, intfire1,intfire2)))

            for treward in range(len(tstops)-1):
             
                  h.CVode().event(tstops[treward], (justSteep.set_exglu_weights, (list_syn_nmda, netcons_s2t)))
                  h.CVode().event(tstops[treward], (justSteep.reset_exglu, (intfire1,intfire2)))
       
                  h.CVode().event(tstops[treward], (justSteep.meanFreq,(vm,tm,treward,tstops[treward],expect,noghte1,train1,stORsh)))
                  
            h.continuerun(tstops[-1])
            return vm,tmfff,tm,Rec_Wg_NMDA,Rec_Wg_GABA,Rec_Wg_NMDA_extraInput,secName, checkD,checkD_cor,checkDnaro,checkDnaro_cor, checkD_i,checkDnaro_i,checkDLtype_cor,checkDLtype,checkD_imin
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
            stORsh=[0.381679389,0.062,1,1,'else','non']
            nDend=2
            num_syn=4*5
            num_syn_inh=4*1
            synNumber=nDend*num_syn
            synNumber_inh=nDend*num_syn_inh
            corticalnoise=144
            batch=12
            epoc=80
            num=batch*epoc
            noghte1=[]
            np.random.seed(corticalnoise)

            taskjob='xor'
            w_exe,w_in=func.makew_inex1(synNumber,self.rank,'spill',num_syn,taskjob)
            w_exe_cor=abs(np.random.normal(0.25, 0.05, size=(1, corticalnoise)))
            # dend_exe=[2, 3, 4, 5, 8, 9, 12, 13, 15, 17, 18, 20, 21, 22, 24, 26, 27, 28, 29, 35, 36, 37, 38, 40, 43, 44, 45, 46, 47, 48, 50, 51, 52, 53, 56, 57]
            dend_exe=[4,5,8,12,15,21,22,24,26,28,35,36,37,41,46,47,51,52,53,57,3,14,17,18,29,38,40,45,48,56,2,13,20,27,43,44]
            tresh_cor=np.ones((1,corticalnoise))*0.02
            tresh=np.random.rand(1,synNumber)*0
            tres=np.ones((1,synNumber))*0.02
            my_new = [i_numx +0 for i_numx in tres]
            for iche in range(synNumber):
                tresh[0][iche]=my_new[0][iche]
           
           # tasks=['YB','RS','YS','RB']

            #expectedOutput={'YB':1,'RS':1,'YS':-1,'RB':-1}
            tasks=['10','01','11','00']

            expectedOutput={'10':1,'01':1,'11':-1,'00':-1}
            #tasks=['10']

            #expectedOutput={'10':1}
            train1,expect=func.makeTrainingBatches(epoc,tasks,batch,expectedOutput,self.rank)
            tstarts,tstops=justSteep.timeStep(100,num)
           # begu,activeIn,expect=func.getActiveInput(train1)
   
        
            dendlst=np.ones((1,2))
            dendlst=self.v
           
     
            dendinput=func.getDendinput(train1,tstarts,synNumber,num_syn)
            dendinput_inh=func.getDendinput(train1,tstarts,synNumber_inh,num_syn_inh)

            atrian_noise=['noise']*num
            dendinput_cor=func.getDendinput(atrian_noise,tstarts,corticalnoise,num_syn)


         
     
            vm,tmfff,tm,Rec_Wg_NMDA,Rec_Wg_GABA,Rec_Wg_NMDA_extraInput,secName,checkD,checkD_cor,checkDnaro,checkDnaro_cor, checkD_i,checkDnaro_i,checkDLtype_cor,checkDLtype,checkD_imin=justSteep.callbb(self,expect,tstops,tstarts,w_exe,w_in,dend_exe,dendinput,dendinput_inh,tresh,noghte1,train1,stORsh,dendlst,w_exe_cor,dendinput_cor,tresh_cor,num_syn,num_syn_inh)

          
             ##############################################################
            return   w_exe ,Rec_Wg_NMDA,w_in,Rec_Wg_GABA,w_exe_cor,Rec_Wg_NMDA_extraInput,noghte1,secName,tmfff,checkD,checkD_cor,checkDnaro,checkDnaro_cor,train1, checkD_i,checkDnaro_i,checkDLtype_cor,checkDLtype,checkD_imin
