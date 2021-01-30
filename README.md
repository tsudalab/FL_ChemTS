# FL_ChemTS

De novo fluorescent molecule generator based on ChemTS (https://github.com/tsudalab/ChemTS) and quantum chemical calculation.  
This program includes an extented version of ChemTS that is parallelizaed by MPI and a filter for Gaussian (https://gaussian.com) to compute absorption and fluorescence wavelength of a molecule.
The Gaussian filter automatically computes the optimized singlet ground (S0) state, vertical singlet excited (S1) state, and relaxed S1 state of a molecule at the density functional theory level.

# Requirements
1. [Gussian](https://gaussian.com)==16
2. [Python](https://www.anaconda.com/download/)>=2.7 
3. [Keras](https://github.com/fchollet/keras) (version 2.0.5) If you installed the newest version of keras, some errors will show up. Please change it back to keras 2.0.5 by pip install keras==2.0.5. 
4. [RDKit](https://anaconda.org/rdkit/rdkit)
5. Intel MPI environment

# Usage

The main python script is fl_chemts/mpi_thread_chemts_tree_vl.py. 

Paralleled search is performed based on Intel MPI environment using fl_chemts/job_sub.sh.

1. cd fl_chemts
2. qsub job_sub

Please setup your MPI and python environment, in fl_chemts/job_sub.sh.

# Dataset
We used 153,253 molecules that contain only H, O, N, and C elements obtained from the ZINC database for trainig of the RNN network.
The file is located at FL_ChemTS/data/.

# Trained RNN network
We trained a RNN network using the above SMILES dataset. The trained network files are located at FL_ChemTS/RNN_model/. 

# License
This package is distributed under the MIT License.
