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

Please cite the following paper when using the program: https://doi.org/10.1063/1.4999097

# List of input files

Here is the list of input files provided and the files that the user needs to edit are bolded. All of the files need to 
be downloaded and placed in the same directory.

1. **`check_state_function.py`**: Defines the reactant and product states. User needs to edit this file if the user 
wants to calculate the fluxes between these states. If the walker is not in one of the defined states, then the function 
should return -1 for that walker. Otherwise, the function should return some other integer value. Make sure each state 
returns an integer from 0 to number of states-1, i.e., if there are two defined states, then they are labeled by 0 and 1, 
respectively.

2. **`clean_up.sh`** (only for MD_Code): Only used when simulation_flag in `parameters.py` is 2, or when the user is 
restarting a CAS simulation that didn't finish post-processing / getting the new collective variables' values from the
walkers that finished running. User needs to edit this file and it will be almost equivalent to the post-processing
part of `simulations.sh`. 

3. **`energy_function.py`** (only for 2D_Toy_Model_Code): Defines the free energy landscape of the system. User needs to 
edit this file if the user wants to run 2D toy model simulations (uses Metropolis algorithm to move walkers)

4. **`execute_pbs.sh`** or **`execute_slurm.sh`**: Shell script to submit to the user's cluster to run a CAS simulation. 
PBS and SLURM scripts are written separately as the name of the files imply.

5. **`functions.py`**: Has all of the functions of the CAS algorithm. User needs to edit a few parts in initialize 
function (highlighted in TODO) that has initial configurations' and simulation files' names (only for MD_Code). Note 
that each initial configurations' files are assumed to be the same and are followed by a number, e.g., minim_0.gro, 
minim_1.gro, and so on. Hence, even if there is only one initial configuration, label it with a number, e.g., 
minim_0.gro.

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
variables are of interest. Note that the necessary simulation files (e.g., initial structure file, simulation setting 
file) need to be provided by the user for the simulations to run. 

12. walker.py: Characterizes the walker object. No need to edit.

# How to edit parameters.py

## 2D_Toy_Model_Code parameters.py

* main_directory: main directory of where the CAS simulation will take place, i.e., where all of the files are.

* balls_flag: either equal to 0 or 1. 

   0: create new macrostates at each step. 1: keep created macrostates (only allowed if enhanced sampling is off or 
   enhanced_sampling_flag = 0).

* flux_flag: either equal to 0 or 1.

   0: off. 1: on. fluxes between pre-defined states will be calculated. the walker's state is determined by 
   check_state_function.py.

* num_states: number of pre-defined states for flux calculation. only needed if flux_flag = 1, otherwise 0.

* enhanced_sampling_flag: either 0 or 1 or 2 or 3.

   0: off. 1: binning walkers if the walkers have some property less or greater than threshold. 2: spectral clustering.
   3: reweighting with equilibrium weights.

* num_balls_limit: maximum number of macrostates in the simulation.

* radius: radii of the Voronoi cells. make sure to have them inside brackets [], even if there's only one radius.

* num_walkers: number of walkers per macrostate.

* grid_dimensions: since this is for a 2D toy model, this is equal to 
[x_lower_bound, x_upper_bound, y_lower_bound, y_upper_bound].

* angle_cvs: list of 0's and/or 1's. make sure to have them inside brackets [], even if there's only one collective 
variable.

   0: if the collective variable is not an angle. 1: if the collective variable is an angle that ranges from 
   -180 deg to 180 deg.

* max_num_steps: maximum number of the CAS simulation steps.

* num_occupied_balls: starting number of initial states for the CAS simulation.

* first_walker: index of first walker in the list.

* last_walker: index of last walker in the list.

* m_steps_per_step: how many times the Metropolis algorithm should be executed per simulation step.

* step_size: how large the step size should be for each walker.

* beta: inverse temperature of the system.

* pbc: either equal to 0 or 1. 

   0: off. 1: periodic boundary conditions on.

##### The next four lines are needed if enhanced_sampling_flag = 1

* less_or_greater_flag: either equal to 0 or 1.  

   0: criteria for binning walkers if the walkers have some property LESS than the threshold.
   
   1: criteria for binning walkers if the walkers have some property GREATER than the threshold.
   
* threshold_values: list of threshold values. if at least one property of the walker has a value less or greater than 
the threshold value, then it is binned coarsely.

* properties_to_keep_track: list of properties of the walker that are compared against threshold values. This can be 
weight and/or some collective variables and multiple duplicates are allowed. If one of them is weight, then type -1. 
Otherwise type the indices of the collective variables, e.g., if there are 3 collective variables and you would like to 
keep track of the last one twice for different threshold values, then type [2, 2] (index starts from 0). If there is 
more than one property to keep track of, then type them sequentially.

* coarse_radius: radius for coarse binning.

##### The next four lines are needed if enhanced_sampling flag = 2

* num_balls_for_sc: minimum number of Voronoi cells present to start performing spectral clustering.

* num_clusters: number of clusters or macrostates (union of Voronoi cells )for k-means part of spectral clustering.

* num_walkers_for_sc: number of walkers for each cluster or macrostate (union of Voronoi cells).

* num_steps_for_sc: number of steps needed to calculate the transition matrix for spectral clustering.

##### The next four lines are needed if enhanced_sampling flag = 3

* initial_step_num_for_eq: starting CAS simulation step number for reweighting.

* num_steps_for_eq: number of steps needed to calculate the transition matrix for reweighting.

* eq_frequency: how many times to perform reweighting.

* num_steps_in_bw: number of steps in between reweightings.

## MD_Code parameters.py

* main_directory: main directory of where the CAS simulation will take place, i.e., where all of the files are.

* initial_configuration_directory: directory where initial configurations are.

* simulation_flag: either equal to 0 or 1 or 2 or 3 or 4.

   0: new simulation. 1: restarting simulation that didn't finish runnning all walkers. 
   2: restarting simulation that didn't finish post-processing walkers. 
   3: restarting simulation that didn't finish binning walkers to Voronoi cells.
   4: restarting simulation that didn't finish resampling.

* balls_flag: either equal to 0 or 1. 

   0: create new macrostates at each step. 1: keep created macrostates (only allowed if enhanced sampling is off or 
   enhanced_sampling_flag = 0).

* flux_flag: either equal to 0 or 1.

   0: off. 1: on. fluxes between pre-defined states will be calculated. the walker's state is determined by 
   check_state_function.py.

* num_states: number of pre-defined states for flux calculation. only needed if flux_flag = 1, otherwise 0.

* enhanced_sampling_flag: either 0 or 1 or 2 or 3.

   0: off. 1: binning walkers if the walkers have some property less or greater than threshold. 2: spectral clustering.
   3: reweighting with equilibrium weights.

* num_balls_limit: maximum number of macrostates in the simulation.

* radius: radii of the Voronoi cells.

* num_walkers: number of walkers per macrostate.

* num_cvs: number of collective variables.

* angle_cvs: list of 0's and/or 1's.

   0: if the collective variable is not an angle. 1: if the collective variable is an angle that ranges from 
   -180 deg to 180 deg.

* initial_step_num: starting CAS simulation step number. it should be changed from 0 when restarting a simulation.

* max_num_steps: maximum number of the CAS simulation steps.

* num_occupied_balls: starting number of initial states for the CAS simulation. 
it should be changed when restarting a simulation.
 
* first_walker: index of first walker in the list.

* last_walker: index of last walker in the list.

##### The next four lines are needed if enhanced_sampling_flag = 1

* less_or_greater_flag: either equal to 0 or 1.  

   0: criteria for binning walkers if the walkers have some property LESS than the threshold.
   
   1: criteria for binning walkers if the walkers have some property GREATER than the threshold.

* threshold_values: list of threshold values. if at least one property of the walker has a value less or greater than 
the threshold value, then it is binned coarsely.

* properties_to_keep_track: list of properties of the walker that are compared against threshold values. This can be 
weight and/or some collective variables and multiple duplicates are allowed. If one of them is weight, then type -1. 
Otherwise type the indices of the collective variables, e.g., if there are 3 collective variables and you would like to 
keep track of the last one twice for different threshold values, then type [2, 2] (index starts from 0). If there is 
more than one property to keep track of, then type them sequentially.

* coarse_radius: radius for coarse binning.

##### The next five lines are needed if enhanced_sampling flag = 2

* num_occupied_clusters: number of occupied clusters or macrostates (unions of Voronoi cells). it could be nonzero when 
restarting a simulation, otherwise 0.

* num_balls_for_sc: minimum number of Voronoi cells present to start performing spectral clustering.

* num_clusters: number of clusters or macrostates (unions of Voronoi cells )for k-means part of spectral clustering.

* num_walkers_for_sc: number of walkers for each cluster or macrostate (union of Voronoi cells).

* num_steps_for_sc: number of simulation steps needed to calculate the transition matrix for spectral clustering.

##### The next four lines are needed if enhanced_sampling flag = 3

* initial_step_num_for_eq: starting CAS simulation step number for reweighting.

* num_steps_for_eq: number of steps needed to calculate the transition matrix for reweighting.

* eq_frequency: how many times to perform reweighting.

* num_steps_in_bw: number of steps in between reweightings.
