# Modelling the Phase Separation of Transcription Factors on DNA

## Overview
This project aims to simulate and understand the phase separation behavior of a transcription factor (TF), referred to as protein A, on a strand of DNA. Specifically, it seeks to predict the distribution of cluster sizes formed by protein A on the DNA, both independently and in the presence of a polymeric intrinsically disordered protein (IDP), designated as protein B. Protein A has non-specific binding sites on the DNA and can self-bind, while protein B serves as a platform with multiple binding sites for protein A. The interactions between these proteins and DNA are governed by three distinct energy scales:

1. **A-DNA Binding**
2. **A-B Binding**

In addition to these binding energies, the concentrations of proteins A and B, as well as the number of binding sites (k) on protein B, comprise the six parameters that define this system. Due to the high dimensionality of this parameter space, brute-force simulations would be computationally prohibitive, making parameter reduction and model simplification essential.

## Objectives
The objective is to implement a one-dimensional Monte Carlo model to predict the distribution of protein A clusters on DNA. The simulation will start with A-DNA binding only, then gradually introduce A-A interactions and the influence of protein B. The main outputs of this model will be the mean cluster size and the cluster size distribution for both proteins A and B on the DNA strand.
MC Simulations will be used to identify the most influential parameters in the cluster growth process and determine a way to limit the parameter ranges. This step will help make the problem more manageable for simulation.



## Project Structure

## Installation and Dependencies

## References

1. **Brackley, C. A., et al.** "Ephemeral protein binding to DNA shapes stable nuclear bodies and chromatin domains." *Biophys. J.* 112:1085-1093 (2017).
2. **Ancona, M., and Brackley, C. A.** "Simulating the chromatin-mediated phase separation of model proteins with multiple domains." *Biophys. J.* 121:2600 - 2612 (2022).

