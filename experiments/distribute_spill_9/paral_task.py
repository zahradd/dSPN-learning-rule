import task as pa
import sys
import os
cur_dir = os.path.abspath(__file__)
model_dir = os.path.join(cur_dir, '..', '..',  '..','model')
updatedmodel_dir =os.path.abspath(model_dir)
sys.path.append(updatedmodel_dir)  

from mpi4py import MPI
import time
import math
import numpy as np
import MSN_builder as build
import pickle
import matplotlib.pyplot as plt
import dFunc as func
import os
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()
sendbuf=[]
root=0
t0 = time.time()


# Specify the file name
bestfitmodel= 'D1_71bestFit_updRheob.pkl'
bestfitmodel_path = os.path.join(updatedmodel_dir, bestfitmodel)
with open(bestfitmodel_path, 'rb') as fpic:
        model_sets = pickle.load(fpic,encoding="latin-1")
    # basic parameters and morphology
par = 'params_dMSN.json'
morphology = 'MSN_morphology_D1.swc'  
par_path = os.path.join(updatedmodel_dir, par)
morphology_path = os.path.join(updatedmodel_dir, morphology)


      
for cell_index in [10]: 
    
    parameters = model_sets[cell_index]['variables']
    # initiate cell -
    cell = build.MSN(params=par_path, \
                     morphology=morphology_path, \
                     variables=parameters)
if rank == 0:
    
    dendlst=np.ones((nprocs,2))
    dendlst1=[[5., 41.]]*nprocs
  #  dendlst1=[[9., 41.],[12., 47.] ,[ 5., 41.] ,[22., 53.],[ 4., 36.],[24., 47.],[ 4., 36.],[24., 37.],[ 9., 57.],[ 9., 41.],[ 4., 57.],[22., 52.],[24., 37.],[ 4., 41.],[26., 51.],[ 9. ,41.],[22., 47.],[ 4., 47.],[15., 51.],[ 9., 52.],[24., 46.],[ 4., 37.],[24., 53.],[ 4. ,41.],[28., 46.],[21., 52.],[28., 46.],[15., 52.],[ 4., 41.],[28., 47.],[ 4., 47.],[24., 52.]]

    
    for i in range(nprocs):
        dendlst[i]=dendlst1[i]
    sendbuf=dendlst


v = comm.scatter(sendbuf, root)
results={}
newClass = pa.justSteep(cell,rank,v)
w_exe ,w_exef,w_in,w_inf,perform,secname,tm,tresh,calMax,train, checkTreshold_i,checkTreshold_imin,checkCalMax_i,checkDLtype,inputListSyn,Output,inputListSyn_in,Output_in,DendPosition =newClass.beforMain()   
wFirst= comm.gather(w_exe, root)
wEnd= comm.gather(w_exef, root)
wInFirst= comm.gather(w_in, root)
wInEnd= comm.gather(w_inf, root)
performs= comm.gather(perform, root)
secNames= comm.gather(secname, root)
tms= comm.gather(tm, root)
treshs= comm.gather(tresh, root)
calMaxs= comm.gather(calMax, root)
train= comm.gather(train, root)
treshs_in= comm.gather(checkTreshold_i, root)
treshsmin_in= comm.gather(checkTreshold_imin, root)
calMaxs_in= comm.gather(checkCalMax_i, root)
checkDLtypes= comm.gather(checkDLtype, root)
inputListSyns= comm.gather(inputListSyn, root)
Outputs= comm.gather(Output, root)
inputListSyns_in= comm.gather(inputListSyn_in, root)
Outputs_in= comm.gather(Output_in, root)
DendPositions= comm.gather(DendPosition, root)
# tmRs= comm.gather(tmR, root)
# vms= comm.gather(vm, root)
# Rec_Vol_Dends= comm.gather(Rec_Vol_Dend, root)
fileName='distribute_spill_9'
# NETWORK_DIR1 = '/scratch/snx3000/'+ os.environ['USER']+'/plasticity.data'
# NETWORK_DIR = 'Result/plasticity.data/'
NETWORK_DIR='/scratch/snx3000/bp000380/Plasticity.data'
INPUT_SAVE_DIR = os.path.join(NETWORK_DIR, str(fileName)+'_WFirst')
# NETWORK_DIR='/scratch/snx3000/$USER/STN_optimize.data'
if rank == 0:
    print('pi computed in {:.3f} sec'.format(time.time() - t0))
    func.save_obj( wFirst, os.path.join(NETWORK_DIR, str(fileName)+'_WFirst'))
    func.save_obj( wEnd,  os.path.join(NETWORK_DIR, str(fileName)+'_Wend') )
    func.save_obj( wInFirst, os.path.join(NETWORK_DIR, str(fileName)+'_WInFirst') )
    func.save_obj(wInEnd,  os.path.join(NETWORK_DIR, str(fileName)+'_WInEnd'))
    func.save_obj( performs, os.path.join(NETWORK_DIR, str(fileName)+'_performance'))
    func.save_obj( secNames, os.path.join(NETWORK_DIR, str(fileName)+'_secNames' ))
    func.save_obj( tms,  os.path.join(NETWORK_DIR, str(fileName)+'_tms'))
    func.save_obj( treshs, os.path.join(NETWORK_DIR, str(fileName)+'_tresh'))
    func.save_obj( calMaxs, os.path.join(NETWORK_DIR, str(fileName)+'_calMaxs'))
    func.save_obj( train,  os.path.join(NETWORK_DIR, str(fileName)+'_train'))
    func.save_obj( treshs_in, os.path.join(NETWORK_DIR, str(fileName)+'_treshs_in'))
    func.save_obj( treshsmin_in,os.path.join(NETWORK_DIR, str(fileName)+'_treshsmin_in'))
    func.save_obj( calMaxs_in, os.path.join(NETWORK_DIR, str(fileName)+'_calMaxs_in'))
    func.save_obj( checkDLtypes,os.path.join(NETWORK_DIR, str(fileName)+'_ltype'))
    func.save_obj( inputListSyns,os.path.join(NETWORK_DIR, str(fileName)+'_inputListSyn'))
    func.save_obj( Outputs,os.path.join(NETWORK_DIR, str(fileName)+'_output'))
    func.save_obj( inputListSyns_in, os.path.join(NETWORK_DIR, str(fileName)+'_inputListSyn_in'))
    func.save_obj( Outputs_in, os.path.join(NETWORK_DIR, str(fileName)+'_output_in') )
    func.save_obj( DendPositions, os.path.join(NETWORK_DIR, str(fileName)+'_DendPositions') )

    # func.save_obj( tmRs, './Results/SteepWithoutInhibitory/epochs1_tmRs' )

    # func.save_obj( vms, './Results/SteepWithoutInhibitory/epochs1_vms' )
    # func.save_obj( Rec_Vol_Dends, './Results/SteepWithoutInhibitory/epochs1_Rec_Vol_Dends' )



   


    