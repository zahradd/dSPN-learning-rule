# dSPN-learning-rule

**This repo is under re-construction for increased clearity**

## Overview

This repository contains the code and data used in the study **"Local, calcium- and reward-based synaptic learning rule that enhances dendritic nonlinearities can solve the nonlinear feature binding problem"** by **Zahra Khodadadi, Daniel Trpevski, Robert Lindroos, and Jeanette Hellgren Kotaleski**.

Using a **biophysically detailed model** of direct pathway Striatal Projection Neurons (dSPNs), from the striatal SPN library:

https://github.com/ModelDBRepository/266775)

we here presents two **synaptic learning rules**:

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

### 1. Navigate to the `experiments/` directory

This directory contains all the different simulation setups used in the paper:

```bash
cd experiments/
```

The main subdirectories include:

- **`clustered_training_experiments/`**  
  - This folder contains experiments focused on clustered synaptic inputs and their effects on plasticity.
  - It includes scripts for running cluster-based plasticity experiments and analyzing results.

- **`distributed_training_experiments/`**  
  - This folder contains experiments investigating distributed synaptic inputs and their impact on plasticity.
  - The simulations in this directory distribute synaptic activity across a wider range of dendritic locations.

### 2. Running a Simulation

To start a simulation, navigate to the appropriate folder. For example, for the clustered training experiments:

```bash
cd clustered_training_experiments/
```

#### **Compiling Mechanisms**
Before running simulations, ensure that the necessary NEURON mechanisms are compiled:

```bash
nrnivmodl ../../mechanisms/
```

#### **Running the Simulation**
Once compiled, run the simulation using:

```bash
python run_cluster_experiment.py
```

For distributed training experiments, navigate to that directory and run:

```bash
cd ../distributed_training_experiments/
python run_plasticity_experiment.py
```

### 3. Plotting Results

After running the simulations, results can be visualized using the provided Jupyter notebooks.

For interactive analysis, open the plotting notebook:

```bash
jupyter-lab plot.ipynb
```

This notebook contains scripts for visualizing synaptic plasticity changes, performance metrics, and other relevant analyses.




## Data Availability

- The processed simulation data is available in the **data/** folder.
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

For questions or issues, please reach out to **[zahra.khodadadi@scilifelab.se](mailto\:zahra.khodadadi@scilifelab.se)** or open an issue in this repository.

Happy coding!

