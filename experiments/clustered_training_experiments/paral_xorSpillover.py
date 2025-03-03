experiment_index = 1  # Change this to run different experiments

import Xor_spillOver as pa
import sys
import os
import json
import numpy as np
import time
import pickle
import MSN_builder as build
import dFunc as func
from mpi4py import MPI
from neuron import h

# MPI initialization
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()
sendbuf = []
root = 0
t0 = time.time()

# Get paths
cur_dir = os.path.abspath(__file__)
model_dir = os.path.join(cur_dir, '..', '..', '..', 'model')
updatedmodel_dir = os.path.abspath(model_dir)
sys.path.append(updatedmodel_dir)  


# Get the absolute path of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))

mechanisms_dir = os.path.join(script_dir, '..', '..', "mechanisms")

# Windows: Look for nrnmech.dll
nrnmech_dll = os.path.join(mechanisms_dir, "nrnmech.dll")

# Linux/macOS: Look for libnrnmech.so
nrnmech_so = os.path.join(mechanisms_dir, "x86_64", ".libs", "libnrnmech.so")

if os.path.exists(nrnmech_dll):  # Windows
    h.nrn_load_dll(nrnmech_dll)
    print(f"Successfully loaded NEURON mechanisms from {nrnmech_dll}")
elif os.path.exists(nrnmech_so):  # Linux/macOS
    h.nrn_load_dll(nrnmech_so)
    print(f"Successfully loaded NEURON mechanisms from {nrnmech_dll}")
else:
    raise FileNotFoundError(f"❌ Mechanism files not found in {mechanisms_dir}. Please compile .mod files using `mknrndll` (Windows) or `nrnivmodl` (Linux/macOS).")

# Load best fit model
bestfitmodel = 'D1_71bestFit_updRheob.pkl'
bestfitmodel_path = os.path.join(updatedmodel_dir, bestfitmodel)
with open(bestfitmodel_path, 'rb') as fpic:
    model_sets = pickle.load(fpic, encoding="latin-1")

# Load morphology and parameters
par = 'params_dMSN.json'
morphology = 'MSN_morphology_D1.swc'  
par_path = os.path.join(updatedmodel_dir, par)
morphology_path = os.path.join(updatedmodel_dir, morphology)

# Load JSON file with experiment settings
json_path = os.path.join(cur_dir, '..', "cluster_experiments.json")
json_path = os.path.normpath(json_path)


with open(json_path, "r") as json_file:
    experiments = json.load(json_file)["experiments"]

# Select experiment (modify index as needed)
selected_experiment = experiments[experiment_index]

# Extract dendritic list
dendlst1 = selected_experiment["dendlst1"]

# Ensure there are enough dendrite entries for all processes
if rank == 0:
    while len(dendlst1) < nprocs:
        dendlst1.append(dendlst1[-1])  # Duplicate last entry if needed

    dendlst = np.array(dendlst1[:nprocs])  # Assign only as many as needed
    sendbuf = dendlst

# Scatter dendrite locations to all processes
v = comm.scatter(sendbuf, root)

# Run the main simulation
results = {}
for cell_index in [10]: 
    parameters = model_sets[cell_index]['variables']
    cell = build.MSN(params=par_path, morphology=morphology_path, variables=parameters)

newClass = pa.justSteep(cell, rank, v)
w_exe, w_exef, w_in, w_inf, w_exe_cor, w_cortical, perform, secname, tm, tresh, tresh_cor, calMax, calMax_cor, train, checkTreshold_i, checkCalMax_i, checkDLtype_cor, checkDLtype, checkTreshold_imin = newClass.beforMain()

# Gather results
wFirst = comm.gather(w_exe, root)
wEnd = comm.gather(w_exef, root)
wInFirst = comm.gather(w_in, root)
wInEnd = comm.gather(w_inf, root)
wCorFirst = comm.gather(w_exe_cor, root)
wCorEnd = comm.gather(w_cortical, root)
performs = comm.gather(perform, root)
secNames = comm.gather(secname, root)
tms = comm.gather(tm, root)
treshs = comm.gather(tresh, root)
treshs_cor = comm.gather(tresh_cor, root)
calMaxs = comm.gather(calMax, root)
calMaxs_cor = comm.gather(calMax_cor, root)
train = comm.gather(train, root)
treshs_in = comm.gather(checkTreshold_i, root)
calMaxs_in = comm.gather(checkCalMax_i, root)
checkDLtype_cors = comm.gather(checkDLtype_cor, root)
checkDLtypes = comm.gather(checkDLtype, root)
treshsmin_in = comm.gather(checkTreshold_imin, root)
# Rec_time= comm.gather(tmR, root)
# Rec_Volt= comm.gather(vm, root)
# Rec_Volt_Dend= comm.gather(Rec_Vol_Dend, root)

# Save results if rank 0
fileName = f"experiment_{experiment_index}"
# Automatically detect if running on an HPC system
if os.path.exists("/scratch"):
    NETWORK_DIR = '/scratch/snx3000/bp000380/Plasticity.data'  # HPC path
else:
    NETWORK_DIR = os.path.join(os.getcwd(), "PlasticityResults")  # Local directory

# Create the experiment-specific folder
experiment_dir = os.path.join(NETWORK_DIR, fileName)
os.makedirs(experiment_dir, exist_ok=True)  # Ensure the directory exists


if rank == 0:
    print(f'Simulation {selected_experiment["name"]} completed in {time.time() - t0:.3f} sec')
    
    func.save_obj(wFirst, os.path.join(experiment_dir, f"{fileName}_WFirst"))
    func.save_obj(wEnd, os.path.join(experiment_dir, f"{fileName}_Wend"))
    func.save_obj(wInFirst, os.path.join(experiment_dir, f"{fileName}_WInFirst"))
    func.save_obj(wInEnd, os.path.join(experiment_dir, f"{fileName}_WInEnd"))
    func.save_obj(wCorFirst, os.path.join(experiment_dir, f"{fileName}_WCorFirst"))
    func.save_obj(wCorEnd, os.path.join(experiment_dir, f"{fileName}_WCorEnd"))
    func.save_obj(performs, os.path.join(experiment_dir, f"{fileName}_performance"))
    func.save_obj(secNames, os.path.join(experiment_dir, f"{fileName}_secNames"))
    func.save_obj(tms, os.path.join(experiment_dir, f"{fileName}_tms"))
    func.save_obj(treshs, os.path.join(experiment_dir, f"{fileName}_tresh"))
    func.save_obj(treshs_cor, os.path.join(experiment_dir, f"{fileName}_treshs_cor"))
    func.save_obj(calMaxs, os.path.join(experiment_dir, f"{fileName}_calMaxs"))
    func.save_obj(calMaxs_cor, os.path.join(experiment_dir, f"{fileName}_calMaxs_cor"))
    func.save_obj(train, os.path.join(experiment_dir, f"{fileName}_train"))
    func.save_obj(treshs_in, os.path.join(experiment_dir, f"{fileName}_treshs_in"))
    func.save_obj(calMaxs_in, os.path.join(experiment_dir, f"{fileName}_calMaxs_in"))
    func.save_obj(checkDLtype_cors, os.path.join(experiment_dir, f"{fileName}_ltype_cor"))
    func.save_obj(checkDLtypes, os.path.join(experiment_dir, f"{fileName}_ltype"))
    func.save_obj(treshsmin_in, os.path.join(experiment_dir, f"{fileName}_treshsmin_in"))
    # func.save_obj( Rec_time, os.path.join(experiment_dir, str(fileName)+'_Rec_time' ))
    # func.save_obj( Rec_Volt, os.path.join(experiment_dir, str(fileName)+'_Rec_Volt' ))
    # func.save_obj( Rec_Volt_Dend, os.path.join(experiment_dir, str(fileName)+'_Rec_Volt_Dend' ))


   


    