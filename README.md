# _FPDock_
Code folder contains the source code of FPDock
 
**Prerequisites:**

Please install the following packages before running the program:

    biopython
    msms
    pandas
    emboss
    python-tk
    matplotlib
    scikit-learn
    pyquaternion
    pdbtools
    scipy

Alternatively, one can use the py3.yml to install all the dependencies.
conda env create -f py3.yml will serve the purpose.

FPDock uses  van der Waals and electrostatic potentials in meetdock scoring function, which can be found in https://github.com/maxibor/meetdock. 

Save the input files in Data folder. Sample inputs are uploaded in this folder. The folder name 1A2B_AB:C indicates that 1A2B is the target protein containing A,B, and C chains. Receptor PDB file is named as AB_model_st.pdb and contains A and B chains. Ligand PDB file is named as C_model_st.pdb and contains C chain. Final result will appear in the folder 'poses'
**To run:**
./fpdock.sh


**Reference**
1. https://github.com/maxibor/meetdock
