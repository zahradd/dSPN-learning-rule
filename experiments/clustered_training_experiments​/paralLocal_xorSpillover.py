
import sys
import os
cur_dir = os.path.abspath(__file__)
model_dir = os.path.join(cur_dir, '..', '..', '..', 'model')
updatedmodel_dir =os.path.abspath(model_dir)
sys.path.append(updatedmodel_dir)  

import time
import numpy as np
import MSN_builder as build
import Xor_spillOver as pa
import pickle
import matplotlib.pyplot as plt
import dFunc as func

t0 = time.time()


bestfitmodel= 'D1_71bestFit_updRheob.pkl'
bestfitmodel_path = os.path.join(updatedmodel_dir, bestfitmodel)
with open(bestfitmodel_path, 'rb') as fpic:
    model_sets = pickle.load(fpic,encoding="latin-1")

# basic parameters and morphology
par = 'params_dMSN.json'
morphology = 'MSN_morphology_D1.swc'  
par_path = os.path.join(updatedmodel_dir, par)
morphology_path = os.path.join(updatedmodel_dir, morphology)
      

# initiate cell
parameters = model_sets[10]['variables']
cell = build.MSN(params=par_path, \
                 morphology=morphology_path, \
                 variables=parameters)

dendlst = [14., 41.]

results={}
newClass = pa.justSteep(cell,0,dendlst)
w_exe ,w_exef,w_in,w_inf,w_exe_cor,w_cortical,perform,secname,tm,tresh,tresh_cor,calMax,calMax_cor,train, checkTreshold_i,checkCalMax_i,checkDLtype_cor,checkDLtype ,checkTreshold_imin,Rec_Vol_Dend,vm,tmR=newClass.beforMain()   
# exglusec=None
# exglu=None
# exnc=None

root=0
wFirst= comm.gather(w_exe, root)
wEnd= comm.gather(w_exef, root)
wInFirst= comm.gather(w_in, root)
wInEnd= comm.gather(w_inf, root)
wCorFirst= comm.gather(w_exe_cor, root)
wCorEnd= comm.gather(w_cortical, root)
performs= comm.gather(perform, root)
secNames= comm.gather(secname, root)
tms= comm.gather(tm, root)
treshs= comm.gather(tresh, root)
treshs_cor= comm.gather(tresh_cor, root)
calMaxs= comm.gather(calMax, root)
calMaxs_cor= comm.gather(calMax_cor, root)
train= comm.gather(train, root)
treshs_in= comm.gather(checkTreshold_i, root)
calMaxs_in= comm.gather(checkCalMax_i, root)
checkDLtype_cors= comm.gather(checkDLtype_cor, root)
checkDLtypes= comm.gather(checkDLtype, root)
treshsmin_in= comm.gather(checkTreshold_imin, root)

Rec_time= comm.gather(tmR, root)
Rec_Volt= comm.gather(vm, root)
Rec_Volt_Dend= comm.gather(Rec_Vol_Dend, root)
fileName='test'
NETWORK_DIR= r"C:\Users\zahra.khodadadi\Documents\GitHub\plasticity\experiments\cluster_0\data"

# if not os.path.exists(dir_name):
#     os.makedirs(dir_name)
# import os
# print("Current Directory:", os.getcwd())

print('pi computed in {:.3f} sec'.format(time.time() - t0))
func.save_obj( wFirst,  os.path.join(NETWORK_DIR, str(fileName)+'_WFirst' ))
func.save_obj( wEnd,  os.path.join(NETWORK_DIR, str(fileName)+'_Wend' ))
func.save_obj( wInFirst,  os.path.join(NETWORK_DIR, str(fileName)+'_WInFirst') )
func.save_obj(wInEnd,  os.path.join(NETWORK_DIR, str(fileName)+'_WInEnd' ))
func.save_obj( wCorFirst,  os.path.join(NETWORK_DIR, str(fileName)+'_WCorFirst' ))
func.save_obj(wCorEnd,  os.path.join(NETWORK_DIR, str(fileName)+'_WCorEnd' ))
func.save_obj( performs,  os.path.join(NETWORK_DIR, str(fileName)+'_performance' ))
func.save_obj( secNames,  os.path.join(NETWORK_DIR, str(fileName)+'_secNames' ))
func.save_obj( tms,  os.path.join(NETWORK_DIR, str(fileName)+'_tms' ))
func.save_obj( treshs,  os.path.join(NETWORK_DIR, str(fileName)+'_tresh' ))
func.save_obj( treshs_cor,  os.path.join(NETWORK_DIR, str(fileName)+'_treshs_cor' ))
func.save_obj( calMaxs,  os.path.join(NETWORK_DIR, str(fileName)+'_calMaxs' ))
func.save_obj( calMaxs_cor,  os.path.join(NETWORK_DIR, str(fileName)+'_calMaxs_cor' ))
func.save_obj( train,  os.path.join(NETWORK_DIR, str(fileName)+'_train' ))
func.save_obj( treshs_in,  os.path.join(NETWORK_DIR, str(fileName)+'_treshs_in' ))
func.save_obj( calMaxs_in,  os.path.join(NETWORK_DIR, str(fileName)+'_calMaxs_in' ))
func.save_obj( checkDLtype_cors,  os.path.join(NETWORK_DIR, str(fileName)+'_ltype_cor' ))
func.save_obj( checkDLtypes, os.path.join(NETWORK_DIR, str(fileName)+'_ltype' ))
func.save_obj( treshsmin_in, os.path.join(NETWORK_DIR, str(fileName)+'_treshsmin_in' ))
func.save_obj( Rec_time, os.path.join(NETWORK_DIR, str(fileName)+'_Rec_time' ))
func.save_obj( Rec_Volt, os.path.join(NETWORK_DIR, str(fileName)+'_Rec_Volt' ))
func.save_obj( Rec_Volt_Dend, os.path.join(NETWORK_DIR, str(fileName)+'_Rec_Volt_Dend' ))



# fig, axs = plt.subplots(3,figsize=(10,12))
# axs[0].plot(Rec_time[0],Rec_Volt_Dend[0][0],  linewidth=2)
# axs[1].plot(Rec_time[0],Rec_Volt_Dend[0][1],  linewidth=2)
# axs[2].plot(Rec_time[0], Rec_Volt[0], linewidth=2)

   


# def plot_data_in_subplots(data_list, titles):
#     fig, axs = plt.subplots(len(data_list), figsize=(10, 12))
#     for i, data in enumerate(data_list):
#         axs[i].plot(tm[0],data, linewidth=2)
#         axs[i].set_title(titles[i])

# # Second figure
# plot_data_in_subplots(
#     [treshs_cor[0][0], calMaxs[0][0], calMaxs_cor[0][0], treshs_in[0][0]],
#     ['treshs_cor']
# )
# plt.show()
# # Second figure
# # plot_data_in_subplots(
# #     [treshs_cor[0][0], calMaxs[0][0], calMaxs_cor[0][0], treshs_in[0][0]],
# #     ['treshs_cor', 'calMaxs', 'calMaxs_cor', 'treshs', 'treshs_in']
# # )

# # # Third figure
# # plot_data_in_subplots(
# #     [calMaxs_in[0][0], checkDLtype_cors[0][0], checkDLtypes[0][0], treshsmin_in[0][0], wFirst[0][0]],
# #     ['calMaxs_in', 'checkDLtype_cors', 'checkDLtypes', 'treshsmin_in', 'wFirst']
# # )

# # # Fourth figure
# # plot_data_in_subplots(
# #     [wEnd[0][0], wInFirst[0][0], wInEnd[0][0], wCorFirst[0][0], wCorEnd[0][0]],
# #     ['wEnd', 'wInFirst', 'wInEnd', 'wCorFirst', 'wCorEnd']
# # )

# # plt.show()

    
