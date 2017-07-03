The regular weighted ensemble method simulation Python code was originally developed by H. Lee, Stanford. The bash 
script for multiple simultaneous simulations was originally developed by S. Ahn, Stanford. The resampling method 
implemented in the code is from E. Darve and J. Izaguirre (Chapter 7 in Schlick's "Innovations in Biomolecular Modeling 
and Simulations."). The Concurrent Adaptive Sampling (CAS) algorithm Python code was developed on top of the regular 
weighted ensemble method code by S. Ahn, Stanford. Improvements on the parameter input file had been made by J. 
Birgmeier, Stanford.

There are two folders provided in this Github repository: **2D_Toy_Model_Code** and **MD_Code**. Each folder provides 
all of the needed code with input files for a particular example, e.g., toymodel1, toymodel2, and penta-alanine. In the 
2D_Toy_Model_Code, **Code_With_Proper_Files** will make a separate folder with many output files for each walker and
**Code_Without_Proper_Files** will only output minimum number of files like flux.txt, total_weight_on_each_ball.txt, 
total_weight.txt, time_record.txt, etc. 

# List of input files

Here is the list of input files provided and the files that the user needs to edit are bolded. All of the files need to 
be downloaded and placed in the same directory.

1. **`check_state_function.py`**: Defines the reactant and product states. User needs to edit this file if the user 
wants to calculate the fluxes between these states.

2. **`clean_up.sh`** (only for MD_Code): Only used when simulation_flag in `parameters.py` is 2, or when the user is 
restarting a CAS simulation that didn't finish post-processing / getting the new collective variables' values from the
walkers that finished running. User needs to edit this file and it will be almost equivalent to the post-processing
part of `simulations.sh`. 

3. **`energy_function.py`** (only for 2D_Toy_Model_Code): Defines the free energy landscape of the system. User needs to 
edit this file if the user wants to run 2D toy model simulations (uses Metropolis algorithm to move walkers)

4. **`execute_pbs.sh`** or **`execute_slurm.sh`**: Shell script to submit to the user's cluster to run a CAS simulation. 
PBS and SLURM scripts are written separately as the name of the files imply.

5. **`functions.py`**: Has all of the functions of the CAS algorithm. User needs to edit a few parts in initialize 
function (highlighted in TODO) that has initial config and simulation files' names (only for MD_Code).

6. global_variables.py: Has all of the global variables used in the CAS algorithm. No need to edit.

7. **`initial_states.txt`**: Defines whether the initial states are in the reactant or product state. User needs to edit 
this file if the user wants to calculate the fluxes between these states. Note that if there is more than one initial 
state, then the user needs to write the initial state on the next line for the next initial state, i.e., if I have two 
initial states 0 and 1, then my initial_states.txt would look like the following:
```
    0
    1
```

8. **`initial_values.txt`**: Writes out the collective variables' values of the initial states. User needs to edit this 
file. Note that if there is more than one initial state, then the user needs to write the initial values on the next 
line for the next initial state, i.e., if I have two initial states x = 1.0, y = 2.0, and x = 3.0, y = 4.0, then my 
initial_values.txt would look like the following:
```
    1.0 2.0
    3.0 4.0
```

9. **`main.py`**: Main function of the CAS simulation. User needs to edit the main directory (highlighted in TODO), 
i.e. change it to the directory where all of the files are.

10. **`parameters.py`**: Parameters of the CAS simulation. User needs to **extensively** edit this file for the user's CAS 
simulation and how to do that is detailed below in the **"How to edit parameters.py"** section.

11. **`simulations.sh`** (only for MD_Code): Handles running walkers in parallel across the provided cores and 
post-processing / getting the new collective variables' values from the walkers that finished running. User needs to 
**extensively** edit this file for the user's CAS simulation, depending on the MD simulation program and what collective
variables are of interest. 

12. walker.py: Characterizes the walker object. No need to edit.

# How to edit parameters.py

Please report any bugs and questions/comments to sahn1@stanford.edu.
