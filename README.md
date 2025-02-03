# Modelling the Phase Separation of Transcription Factors on DNA

## Overview  
This project employs a one-dimensional **Monte Carlo simulation** to investigate the **cluster size distribution** of a transcription factor (**protein A**) on a linear DNA strand. The primary objective is to evaluate whether introducing a second protein—an **intrinsically disordered protein (IDP, protein B)**—induces a **phase transition** in the clustering of protein A.

The simulation begins with **protein A binding exclusively to DNA**, followed by the gradual incorporation of **protein B** to analyze its impact. The study explores **phase separation** through **dimensionality reduction** of the parameter space.

### Key Outputs:
- **Cluster size distribution** of protein A.
- **Residence times** of transcription factors (TFs).
- Identification of **critical parameters** affecting cluster formation.
- **Optimization** of parameter ranges for computational efficiency.
- Evidence of **scale-invariance** and phase separation behavior.

---

## Methods  
This **Monte Carlo simulation**, implemented in **Python**, models discrete steps where proteins bind and unbind from specific sites along a **linear DNA chain**. The core functions include:

- **Monte Carlo steps** for protein binding/unbinding events.
- **Energy computations** for binding/unbinding processes.
- **State updates** for protein-DNA interactions over time.

To optimize performance, **numpy arrays** are used for efficient memory management. When arrays are unsuitable, **lists** are utilized and deleted after use to free memory. Simulations were conducted on the **Helvetios computing system**.

---

## Objectives  
The goal is to **predict the distribution of protein A clusters** on DNA using a **one-dimensional Monte Carlo model**. The simulation follows these steps:

1. **Initial state:** Only protein A binds to DNA.
2. **Progressive interactions:** A-A interactions are introduced.
3. **Incorporation of protein B:** Analyze its effect on A clustering.

### Main Outputs:
- **Mean cluster size** and **cluster size distribution** of proteins A and B.
- **Parameter identification** for efficient simulations.
- **Phase separation characteristics** within the system.

---

## Simulation Parameters  

| **Parameter** | **Description** | **Values Investigated** | **Comments** |
|--------------|----------------|------------------------|--------------|
| `N` | Total number of DNA binding sites | 3000 (short runs), 6000 (long runs) | Promoter regions are typically 100-1000 bp |
| `α` | Density of A proteins on DNA | 0.15, 0.3 (short), 0.15 (long) | - |
| `n_A` | Number of A proteins | Computed as `int(α * N)` | 450, 900 (short), 900 (long) |
| `n_B` | Number of B proteins | 100 (short), 200 (long) | - |
| `k` | Number of B binding sites | 2, 5, 10 (short), 2, 10 (long) | - |
| `L` | Distance between B binding sites | 10 (both short and long) | - |
| `stop_time` | Simulation end time | `2 × 10⁶` (short), `4 × 10⁶` (long) | - |
| `ignoring_time` | Time before sampling begins | `1 × 10⁶` (short), `2 × 10⁶` (long) | - |
| `m` | Sampling frequency | 50 (both short and long) | Must be a divisor of `stop_time - ignoring_steps` |
| `E_ad` | Energy of A binding to DNA | 0, 1, 3 (short), 1 (long) | In units of `k_BT` |
| `E_aa` | Interaction energy between A proteins | 0, 2 (short), 2 (long) | - |
| `E_ab` | Effect of B on A unbinding | Set to infinity | - |
| `E_ba` | Effect of A on B unbinding | 0, 1, 2, 4 (short), 0, 4 (long) | - |

---

## Analysis and Plots  

| **Plot** | **Description** | **X-Axis** | **Y-Axis** | **Purpose** |
|---------|---------------|------------|------------|------------|
| **First binding time of As that don’t unbind** | Tracks TFs that remain bound | TF index | Binding time | Identify persistent TFs |
| **Cluster size frequency** | Histogram of A cluster sizes | Cluster size | Frequency | Analyze cluster distributions |
| **Mean residence time per A** | Histogram of residence times | Mean residence time | Frequency | Assess binding dynamics |
| **Std Dev of residence time per A** | Scatter plot of standard deviation | TF index | Std Dev of residence time | Show variability in binding |
| **Fraction of B-occupied sites over time** | Tracks B binding stability | Time step | Occupied site fraction | Show B protein dynamics |
| **Bs bound to As at two time points** | Histograms of B binding at different times | DNA site index | Frequency | Analyze B clustering |
| **Number of binding events per A** | Histogram of A binding occurrences | Number of events | Frequency | Examine system fluidity |
| **Number of bound As over time** | Time series of A binding | Time step | Bound A count | Determine steady-state behavior |
| **Mean cluster size vs. Max cluster size** | Compares cluster size distributions | `E_aa` | Cluster size | Assess cluster formation |
| **Log(mean residence time) vs. `E_ad`** | Relationship between binding energy and stability | `E_ad` | ln(Mean Residence Time) | Understand TF binding energetics |

---

## Repository Structure  

📂 TF_DNA  
├── 📁 MC # Monte Carlo simulation 

│   ├── 📄 MC_simulation.py # Main file to run with the general structure of the simulation

│   ├── 📄 steps_MC_simulation.py #File with the steps of the MC simulation (step fpr the case with only protein A or both proteins)

│   ├── 📄 events_for_MC_steps.py # File with the specific events called in the steps

│   ├── 📄 general_functions.py # File with the functions used for plots or creating files 

│   ├── 📄 cluster_class.py # Started doing a cluster class- incomplete 

│   ├── 📁 Simulations_protein_A  #Results for just the A case 

│   │   ├── 📁 α_0.15 
 
│   │   ├── 📁 α_0.3  

│   ├── 📁 Simulations_proteins_A_B  #Results for both A and B case 

│   │   ├── 📁 α_0.15  

│   │   │   ├── 📁 N_3000

│   │   │   │   ├── 📁 nB_100_K_0_L_0

│   │   │   │   ├── 📁 nB_100_K_2_L_10

│   │   │   │   ├── 📁 nB_100_K_5_L_10

│   │   │   │   ├── 📁 nB_100_K_10_L_10

│   │   │   ├── 📁 N_6000

│   │   │   │   ├── 📁 nB_200_K_2_L_10

│   │   │   │   ├── 📁 nB_200_K_10_L_10

│   │   ├── 📁 α_0.3  

│   │   │   ├── 📁 nB_100_K_2_L_10

│   │   │   ├── 📁 nB_100_K_5_L_10

│   │   │   ├── 📁 nB_100_K_10_L_10

├── 📄 report.pdf                 # Detailed discussion of results  

├── 📄 README.md                  # Project documentation


