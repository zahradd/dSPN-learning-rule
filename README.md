# dSPN-learning-rule

**This repo is under re-construction for increased clarity**

## Overview

This repository contains the code and data used in the study **"Local, calcium- and reward-based synaptic learning rule that enhances dendritic nonlinearities can solve the nonlinear feature binding problem"** by **Zahra Khodadadi, Daniel Trpevski, Robert Lindroos, and Jeanette Hellgren Kotaleski**.

Using a **biophysically detailed model** of direct pathway Striatal Projection Neurons (dSPNs) from the striatal SPN library:

[ModelDB Repository (266775)](https://github.com/ModelDBRepository/266775)

we present two **synaptic learning rules**:

1. **An excitatory local synaptic learning rule** based on calcium dynamics and reward signals.
2. **An inhibitory plasticity rule** that enhances dendritic compartmentalization to improve feature discrimination.

The model demonstrates how a **single neuron can solve the Nonlinear Feature Binding Problem (NFBP)** by leveraging dendritic nonlinearities, metaplasticity, and inhibitory plasticity.

## Features

- **Multicompartmental dSPN Model:** Simulates dendritic plateau potentials and nonlinear synaptic integration.
- **Calcium-Based Synaptic Learning Rule:** Implements a biologically plausible LTP/LTD mechanism based on NMDA and L-type calcium channels.
- **Reward-Driven Synaptic Plasticity:** Uses dopamine signals to modify synaptic strengths dynamically.
- **Metaplasticity:** Regulates learning stability and maintains synaptic weights within a physiological range.
- **Inhibitory Plasticity Rule:** Enhances dendritic compartmentalization to improve feature discrimination.
- **Simulations for NFBP Learning:** Includes both clustered and randomly distributed synapses to test learning capabilities.

## Running Simulations

### Step 1: Compile Mechanisms

```bash
cd mechanisms/
nrnivmodl
cd ..
```

### Step 2: Navigate to the `experiments/` directory

```bash
cd experiments/
```

This directory contains all the different simulation setups used in the study.

The two subdirectories include:

- **`clustered_training_experiments`**  
  - Focuses on clustered synaptic inputs and their effects on plasticity.  
  - A **single JSON file (`experiment_config.json`)** defines multiple experiments, which can be selected in `run_cluster_experiment.py`.  
  - Includes scripts for running cluster-based plasticity experiments and analyzing results.

- **`distributed_training_experiments`**  
  - Investigates distributed synaptic inputs and their impact on plasticity.  
  - **Each experiment has its own JSON configuration file**, stored in the `experiments_config/` directory, which can be selected when running `run_plasticity_experiment.py`.  
  - These simulations introduce **more randomness** in synaptic placement and activation, allowing for greater flexibility in modifying the neuron’s behavior.

### Step 3: Run a Simulation

#### Clustered Training Experiments

To start a **clustered training** experiment, run:

```bash
cd clustered_training_experiments/
python run_cluster_experiment.py -e 1
```

(or use MPI for parallel execution)

```bash
mpiexec -n 1 python run_cluster_experiment.py -e 1
```

This will run **experiment 1** (as specified in `experiment_config.json`). 
You can select different experiments (1-8) as follows:

- **1:** Figure 3 example (no inhibition)
- **2:** Figure 4 combination (no inhibition)
- **3:** Figure 3 example (with inhibition)
- **4:** Figure 4 combination (with inhibition)
- **5:** Figure 4 distal example (with inhibition)
- **6:** Figure 4 proximal example (with inhibition)
- **7:** Figure 4 distal example (no inhibition)
- **8:** Figure 4 proximal example (no inhibition)

#### Distributed Training Experiments

To start a **distributed training** experiment, run:

```bash
cd distributed_training_experiments/
python run_plasticity_experiment.py -e 1
```

(or use MPI for parallel execution)

```bash
mpiexec -n 1 python run_plasticity_experiment.py -e 1
```

This will run **experiment 1** (as specified in the `experiment_config/` directory). 
You can select different experiments (1-2) as follows:

- **1:** Distributed training experiment with inhibition
- **2:** Distributed training experiment without inhibition

## Plotting Results

After running simulations, results can be visualized using the provided Jupyter notebooks.

For interactive analysis, open the plotting notebook:

```bash
jupyter-lab plot.ipynb
```

This notebook contains functions for visualizing synaptic plasticity changes, performance metrics, and other relevant analyses.

## Data Availability

- Raw data for large simulations can be requested from the first author (see contact information below).
- Model components adapted from:
  - **Lindroos & Hellgren Kotaleski (2021)**: ModelDB accession **266775**
  - **Trpevski et al. (2023)**: ModelDB accession **2017143**

## Citation

If you use this code, please cite the following paper:

```
@article{khodadadi2025dspn,
  author = {Khodadadi, Zahra and Trpevski, Daniel and Lindroos, Robert and Hellgren Kotaleski, Jeanette},
  title = {Local, calcium- and reward-based synaptic learning rule that enhances dendritic nonlinearities can solve the nonlinear feature binding problem},
  journal = {eLife},
  year = {2025}
}
```

## License

This project is licensed under the **GNU General Public License v3.0**. See the LICENSE file for details.

## Contact

For questions or issues, please reach out to **[zahra.khodadadi@scilifelab.se](mailto:zahra.khodadadi@scilifelab.se)** or open an issue in this repository.

Happy coding!
