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

First, start by **compiling mechanisms**

```bash
cd mechanisms/
nrnivmodl
cd ..
```  

then

### Navigate to the `experiments/` directory

```bash
cd experiments/
```

This directory contains all the different simulation setups used in the paper.

The two subdirectories include:

- **`clustered_training_experiments`**  
  - Contains experiments focused on clustered synaptic inputs and their effects on plasticity.  
  - A **single JSON file (`experiment_config.json`)** defines multiple experiments, which can be selected within `run_cluster_experiment.py`.  
  - Includes scripts for running cluster-based plasticity experiments and analyzing results.  

- **`distributed_training_experiments`**  
  - Contains experiments investigating distributed synaptic inputs and their impact on plasticity.  
  - **Each experiment has its own JSON configuration file**, stored in the `experiments_config/` directory, which can be chosen when running `run_plasticity_experiment.py`.  
  - These simulations introduce **more randomness** in synaptic placement and activation, allowing for greater flexibility in modifying the network’s behavior.  


### Start a simulation

To start a simulation, navigate to the appropriate folder:

- **For clustered training experiments**, execute:

  ```bash
  cd clustered_training_experiments/
  python run_cluster_experiment.py -i 1
  ```
  
  which will run experiment 1 (as specified in `experiment_config.json`). 
  
  You can select different experiments by changing the argument from 1-8:
  
  - 1: Figure 3 example, without (default)
  - 2: Figure 4 combination
  - 3: Figure 5, with inhibition
  - 4: ...
  - 5: ...
  - 6: ...
  - 7: ...
  - 8: ...

- **For distributed training experiments**, execute:

  ```bash
  cd distributed_training_experiments/
  python run_plasticity_experiment.py
  ```

  - The experiment configurations are stored as separate JSON files in `experiments_config/`.
  - The specific experiment is selected when executing `run_plasticity_experiment.py`.

## Plotting Results

After a simulations is run, results can be visualized using the provided Jupyter notebooks.

For interactive analysis, open the plotting notebook:

```bash
jupyter-lab plot.ipynb
```

This notebooks contains functions for visualizing synaptic plasticity changes, performance metrics, and other relevant analyses.



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

