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

### 1.

Navigate to the **experiments/** directory from were you can find all the different simulations run in the paper:

```bash
cd  experiments/
```

- **cluster_inh_spill_combine_2** comments describing the simulation here please
- **cluster_inh_spill_example_3**
- **distribute_spill_inh_10**
- **xor_dist2_inh_7**
- **xor_dist1_inh_6**

### 2. 

To start a simulation go to, for example, the first folder:

```bash
cd  cluster_inh_spill_combine_2
```

compile mechanisms in this folder:

```bash
nrnivmodl ../../mechanisms/
```

and then run:

```bash
python paralLocal_xorSpillover.py
```

### 3. plotting

Plotting is done in the notebook **plot.ipynb**. For interactive analysis open it using:

```bash
jupyter-lab plot.ipynb
```


## Data Availability

- The processed simulation data is available in the **data/** folder.
- Model components adapted from:
  - **Lindroos & Hellgren Kotaleski (2021)**: ModelDB accession **266775**
  - **Trpevski et al. (2023)**: ModelDB accession **2017143**

## Citation

If you use this code, please cite the following paper:

```
@article{khodadadi2024dspn,
  author = {Khodadadi, Zahra and Trpevski, Daniel and Lindroos, Robert and Hellgren Kotaleski, Jeanette},
  title = {Local, calcium- and reward-based synaptic learning rule that enhances dendritic nonlinearities can solve the nonlinear feature binding problem},
  journal = {Journal TBD},
  year = {2024}
}
```

## Contact

For questions or issues, please reach out to **[zahra.khodadadi@scilifelab.se](mailto\:zahra.khodadadi@scilifelab.se)** or open an issue in this repository.

## License

This project is licensed under the **GNU General Public License v3.0**. See the LICENSE file for details.

---

This repository enables full reproducibility of the results presented in the manuscript, providing detailed documentation for model execution and result interpretation. Happy coding!

