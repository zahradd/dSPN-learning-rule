# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 20:39:04 2016

@author: daniel
"""
#from math import sqrt
results_directory = './results'

#-----------------------------------------------------------#
#      1. General recording and simulation parameters       #
#-----------------------------------------------------------#
step = 20.0
record_step = 0.1
g_record_step = 5
skip_first_x_ms = 100

nrn_dots_per_1ms = 1.0/record_step
time_to_avg_over = 20 # in seconds

simtime = 400
training_mode = 'supra'

num_trials = 20
NUMBER_OF_PROCESSES = 7

#-----------------------------------#
#      2. Synaptic parameters       #
#-----------------------------------#
esyn_tau = 6
isyn_tau = 6
isyn_plateau_tau = 6#87
e_esyn = 0
e_gaba = -60
erate = 0.8
irate = 0.8
pos = 0.55

# NMDA parameters
Mg =  1.4#1.0 #
alpha = 0.099# 0.062
eta = 1.0/18# 1.0/3.57

g_ramp_max = 0.000255
nmda_ampa_ratio = 1
gAMPAmax = 0.1e-3
gNMDAmax = gAMPAmax*nmda_ampa_ratio
gGABAmax = 1.5e-3
g_expsyn_max =  0.1e-3
g_inhexpsyn_max = gGABAmax

gAMPAmax_plateau = 0.1e-3 
gNMDAmax_plateau = 1.0e-3 
gGABAmax_plateau = 1.5e-3
nmda_ampa_ratio = gNMDAmax_plateau/gAMPAmax_plateau
ratio_glutamate_syn = 1.0

gAMPAmax_pf = 0.05e-3
gNMDAmax_pf = 0.75e-3 

glu_thresh1 = 0.1
glu_thresh2 = 0.1

tNMDAon = 2.76
tNMDAoff = 115.5
tau1_NMDA = 2.76
tau2_NMDA = 115.5
tau1_exp2syn = 1.9
tau2_exp2syn = 4.8
tau1_inhexp2syn = 1
tau2_inhexp2syn = 10

v_thresh = -50
glu_thresh = 0.06
#-----------------------------------------#
#      3. Synaptic input parameters       #
#-----------------------------------------#

plateau_syn_rate = 400
plateau_burst_start = 100
plateau_burst_end = 130
plateau_cluster_size = 20
plateau_cluster_size_max = 20
cluster_start_pos = 0.55
cluster_end_pos = 0.70

pf_input_rate = 20.0
pf_input_start = plateau_burst_start
pf_input_end = plateau_burst_end
pf_input_size = 4

inhibitory_syn_rate = 85.0
inhibitory_burst_start = 100
inhibitory_burst_end = 160
inhibitory_cluster_size = 1

distributed_input_rate = 1000.0/40
distributed_input_start = 330
distributed_input_end = 360
distributed_input_size = 62

deterministic_interval = 1.5

ramp_syn_rate = 100.0
ramp_slope = 0.5  # Hz/ms
ramp_burst_start = 200
ramp_burst_end = 1000

e_interval = 1.0/erate*(10**3)
i_interval = 1.0/irate*(10**3)
plateau_syn_interval = 1.0/plateau_syn_rate*(10**3)
distributed_input_interval = 1.0/distributed_input_rate*(10**3)
ramp_syn_interval = 1.0/ramp_syn_rate*(10**3)
inhibitory_syn_interval = 1.0/inhibitory_syn_rate*(10**3)

#--------------------------------#
#      4. Spine parameters       #
#--------------------------------#
head_L = 0.5
head_diam = 0.5
neck_L = 0.5
neck_diam = 0.125
neck_Ra = 1130.0
head_Ra = 150

#-----------------------------------------------------------#
#      5. XOR problem and adaptive synapse parameters       #
#-----------------------------------------------------------#

event_times = [200, 500, 900]

save_weights_file = 'syn_weights.dat'
load_weights_file = 'syn_weights_xor_hom.dat'
training_set_size_per_group = 20
training_set_size = training_set_size_per_group*4
training_input_length = 30
first_training_input_start = 300
time_to_reward = 100 - training_input_length
reward_length = 20
session_length = 500
test_set_size_per_group = 20
test_set_size = test_set_size_per_group*4

low_rate = 0.001
high_rate = 50
LTP_factor = 1.5
LTD_factor = 0.5
LTD_factor_NMDA = 0.5
l_thresh_LTP = 40
l_thresh_LTD = 40
learning_rate_w = 0.0125

random_weights = False
dense = False
read_input_config_from_file = True
input_config_file = 'xor_inputs_to_dends.dat'#'xor_dense.dat'
input_dends = [3, 22, 26, 35]#[3, 5, 8, 12, 15, 22, 26, 35, 41, 47, 53, 57]# 
pf_input_dends = [22, 26, 35]#[12, 22, 26, 35, 41, 53] # 

if dense == False and random_weights == False:
    thresh_LTP_max = 0.45
    thresh_LTP = 0.07#0.13
    thresh_LTP_min = 0.04
    learning_rate_thresh_LTP = 0.005

    thresh_LTD_max = 0.022
    thresh_LTD = 0.005#0.001
    thresh_LTD_min = 0.001#0.0005
    learning_rate_thresh_LTD = 0.005

    xor_input_size = 10
    p_conn = 0.6#1.0
    p_group = 0.5#1.0

elif dense == True and random_weights == True:
    thresh_LTP_max = 0.45
    thresh_LTP = 0.07#0.13
    thresh_LTP_min = 0.04
    learning_rate_thresh_LTP = 0.005
    
    thresh_LTD_max = 0.025
    thresh_LTD = 0.005#0.001
    thresh_LTD_min = 0.0005
    learning_rate_thresh_LTD = 0.005

    xor_input_size = 7
    input_dends = [3, 5, 8, 12, 15, 22, 26, 35, 41, 47, 53, 57]# 
    pf_input_dends = []#[12, 22, 26, 35, 41, 53] # 

    p_conn = 1.0
    p_group = 1.0
    
    
input_dends_dict = {'loc': [3, 22, 26, 35],
                    'pos': [0.55, 0.55, 0.55, 0.05],
                    'start': [plateau_burst_start]*4,
                    'end': [plateau_burst_end]*4 }

pf_input_dends_dict = {'loc': [22, 26, 35],
                    'pos': [0.55, 0.55, 0.05],
                    'start': [plateau_burst_start]*3,
                    'end': [plateau_burst_end]*3 }
                    
low_interval = 1.0/low_rate*(10**3)
high_interval = 1.0/high_rate*(10**3)

#simtime = first_training_input_start + training_set_size*session_length
#simtime = first_training_input_start + test_set_size*session_length
#simtime = first_training_input_start + training_set_size*(training_input_length+
#            time_to_reward + reward_length)

#----------------------------------#
#      6. Plotting parameters      #
#----------------------------------#

#-----------------------------------#
#      6.1. For XOR experiment      #
#-----------------------------------#

#-------------------------------------------------------#
#      7. Miscellaneous and parameter dictionaries      #
#-------------------------------------------------------#

dends_per_plot = 2
scale_conductance = 1000


#params = {
#        'erate' : erate,
#        'irate' : irate,
#        'e_esyn' : e_esyn,
#        'g_expsyn_max' : g_expsyn_max,
#        'g_inhexpsyn_max' : g_inhexpsyn_max,
#        'erate' : erate,
#        'irate' : irate,
#        'esyn_tau' : esyn_tau,
#        'isyn_tau' : isyn_tau,
#        'e_interval' : e_interval,
#        'i_interval' : i_interval
#}

par = {
        'gbar_naf_somatic': 9.0,
        'gbar_naf_axonal': 9.0
       }

       
nmda = {
        'gNMDAmax': gNMDAmax ,
        'gNMDAmax_plateau': gNMDAmax_plateau,
        'tcon': tNMDAon,
        'tcoff': tNMDAoff,
        'ratio': ratio_glutamate_syn
       }       
       
exp2syn = {
        'e': e_esyn ,
        'tau1': tau1_exp2syn ,
        'tau2': tau2_exp2syn
       }
