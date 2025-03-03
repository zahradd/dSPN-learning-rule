import plasticity_model as pa
import sys
import os
import json
import time
import numpy as np
import pickle
import MSN_builder as build
import dFunc as func
import matplotlib.pyplot as plt
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()
root = 0
t0 = time.time()

# Define paths
cur_dir = os.path.dirname(os.path.abspath(__file__))
print('gggggggggggggggggggggggggggggggg',cur_dir)
model_dir = os.path.abspath(os.path.join(cur_dir, '..', '..', 'model'))
print('gggggggggggggggggggggggggggggggg',model_dir)

sys.path.append(model_dir)

# Load best fit model
bestfitmodel = 'D1_71bestFit_updRheob.pkl'
bestfitmodel_path = os.path.join(model_dir, bestfitmodel)
with open(bestfitmodel_path, 'rb') as fpic:
    model_sets = pickle.load(fpic, encoding="latin-1")

# Load morphology and parameters
par_path = os.path.join(model_dir, "params_dMSN.json")
morphology_path = os.path.join(model_dir, "MSN_morphology_D1.swc")

# Find JSON file dynamically
experiments_dir = os.path.join(cur_dir, "experiments_config")

# Use the first JSON file found (or modify to loop over multiple files)
json_filename =  "experiment1.json" # Picks the experiment
json_path = os.path.join(experiments_dir, json_filename)

# Extract the experiment name from the JSON filename (removing '.json')
experiment_name = os.path.splitext(json_filename)[0]

# Load experiment configuration
with open(json_path, "r") as json_file:
    selected_experiment = json.load(json_file)

# Initialize model for cell_index 10
for cell_index in [10]: 
    parameters = model_sets[cell_index]['variables']
    cell = build.MSN(params=par_path, morphology=morphology_path, variables=parameters)

if rank == 0:
    sendbuf = [f"Data for process {i}" for i in range(nprocs)]  # Example data
else:
    sendbuf = None  # Other processes don't need sendbuf


# MPI communication
v = comm.scatter(sendbuf, root)
results = {}
newClass = pa.justSteep(cell, rank, v)

# Run experiment
w_exe, w_exef, w_in, w_inf, perform, secname, tm, tresh, calMax, train, checkTreshold_i, checkTreshold_imin, \
calMax_i, checkDLtype, inputListSyn, Output, inputListSyn_in, Output_in, DendPosition = newClass.beforMain(selected_experiment)

# Gather results
wFirst = comm.gather(w_exe, root)
wEnd = comm.gather(w_exef, root)
wInFirst = comm.gather(w_in, root)
wInEnd = comm.gather(w_inf, root)
performs = comm.gather(perform, root)
secNames = comm.gather(secname, root)
tms = comm.gather(tm, root)
treshs = comm.gather(tresh, root)
calMaxs = comm.gather(calMax, root)
train = comm.gather(train, root)
treshs_in = comm.gather(checkTreshold_i, root)
treshsmin_in = comm.gather(checkTreshold_imin, root)
calMaxs_in = comm.gather(calMax_i, root)
checkDLtypes = comm.gather(checkDLtype, root)
inputListSyns = comm.gather(inputListSyn, root)
Outputs = comm.gather(Output, root)
inputListSyns_in = comm.gather(inputListSyn_in, root)
Outputs_in = comm.gather(Output_in, root)
DendPositions = comm.gather(DendPosition, root)

# Save results
fileName = f"{experiment_name}"  # Use JSON filename as the experiment name

# Automatically detect if running on an HPC system
if os.path.exists("/scratch"):
    NETWORK_DIR = '/scratch/snx3000/bp000380/Plasticity.data'  # HPC path
else:
    NETWORK_DIR = os.path.join(os.getcwd(), "PlasticityResults")  # Local directory

# Create the experiment-specific folder
experiment_dir = os.path.join(NETWORK_DIR, fileName)
os.makedirs(experiment_dir, exist_ok=True)  # Ensure the directory exists

if rank == 0:
    func.save_obj(wFirst, os.path.join(experiment_dir, f"{fileName}_WFirst"))
    func.save_obj(wEnd, os.path.join(experiment_dir, f"{fileName}_Wend"))
    func.save_obj(wInFirst, os.path.join(experiment_dir, f"{fileName}_WInFirst"))
    func.save_obj(wInEnd, os.path.join(experiment_dir, f"{fileName}_WInEnd"))
    func.save_obj(performs, os.path.join(experiment_dir, f"{fileName}_performance"))
    func.save_obj(secNames, os.path.join(experiment_dir, f"{fileName}_secNames"))
    func.save_obj(tms, os.path.join(experiment_dir, f"{fileName}_tms"))
    func.save_obj(treshs, os.path.join(experiment_dir, f"{fileName}_tresh"))
    func.save_obj(calMaxs, os.path.join(experiment_dir, f"{fileName}_calMaxs"))
    func.save_obj(train, os.path.join(experiment_dir, f"{fileName}_train"))
    func.save_obj(treshs_in, os.path.join(experiment_dir, f"{fileName}_treshs_in"))
    func.save_obj(treshsmin_in, os.path.join(experiment_dir, f"{fileName}_treshsmin_in"))
    func.save_obj(calMaxs_in, os.path.join(experiment_dir, f"{fileName}_calMaxs_in"))
    func.save_obj(checkDLtypes, os.path.join(experiment_dir, f"{fileName}_ltype"))
    func.save_obj(inputListSyns, os.path.join(experiment_dir, f"{fileName}_inputListSyn"))
    func.save_obj(Outputs, os.path.join(experiment_dir, f"{fileName}_output"))
    func.save_obj(inputListSyns_in, os.path.join(experiment_dir, f"{fileName}_inputListSyn_in"))
    func.save_obj(Outputs_in, os.path.join(experiment_dir, f"{fileName}_output_in"))
    func.save_obj(DendPositions, os.path.join(experiment_dir, f"{fileName}_DendPositions"))
    