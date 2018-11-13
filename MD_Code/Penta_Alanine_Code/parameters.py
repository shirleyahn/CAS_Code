main_directory='/scratch/users/sahn1/Penta_Alanine'  # main directory of where the CAS simulation will take place,
                                                    # i.e., where all of the files are.
initial_configuration_directory='/scratch/users/sahn1/Penta_Alanine/InitConfig'  # where initial configurations are.

simulation_flag=0  # 0: new simulation. 1: restarting simulation that didn't finish running all walkers.
                   # 2: restarting simulation that didn't finish post-processing walkers.
                   # 3: restarting simulation that didn't finish binning walkers to Voronoi cells.
                   # 4: restarting simulation that didn't finish resampling.
balls_flag=0  # 0: create new balls at each step. 1: keep created balls (only allowed if enhanced sampling is off).
flux_flag=1  # 0: off. 1: on. fluxes between pre-defined states will be calculated. the walker's state is determined by
             # check_state_function.py.
num_states=2  # number of pre-defined states for flux calculation. only needed if flux_flag = 1, otherwise 0.
enhanced_sampling_flag=0  # 0: off. 1: binning walkers if the walkers have some property less or greater than threshold.
                          # 2: spectral clustering. 3: reweighting with equilibrium weights.

num_balls_limit=600  # maximum number of macrostates in the simulation.
radius=[80.0, 80.0, 80.0, 80.0, 80.0, 80.0]  # radii of the Voronoi cells.
num_walkers=10  # number of walkers per macrostate.
num_cvs=6  # number of collective variables.
angle_cvs=[1, 1, 1, 1, 1, 1]  # 0: if the cv is not an angle. 1: if the cv is an angle.

initial_step_num=0  # starting CAS simulation step number. it should be changed from 0 when restarting a simulation.
max_num_steps=50  # maximum number of the CAS simulation steps.
num_occupied_balls=1  # starting number of initial states for the CAS simulation.
                      # it should be changed when restarting a simulation.
first_walker=0  # index of first walker in the list.
last_walker=0  # index of last walker in the list.

### the next four lines are needed if enhanced_sampling_flag = 1 ###
less_or_greater_flag=0  # 0: criteria for binning walkers if the walkers have some property LESS than the threshold.
                        # 1: criteria for binning walkers if the walkers have some property GREATER than the threshold.
threshold_values=[0.0, 10.0]  # list of threshold values. if at least one property of the walker has a value less or
                              # greater than the threshold value, then it is binned coarsely.
properties_to_keep_track=[0, 0]  # list of properties of the walker that are compared against threshold values. this can
                                 # be weight and/or some collective variables and multiple duplicates are allowed.
                                 # if one of them is weight, then type -1. otherwise type the indices of the collective
                                 # variables, e.g., if there are 3 collective variables and you would like to keep track
                                 # of the last one twice for different threshold values, then type [2, 2] (index starts
                                 # from 0). if there is more than one property to keep track of, then type them sequentially.
coarse_radius=[0.5, 0.5]  # radius for coarse binning.

### the next five lines are needed if enhanced_sampling flag = 2 ###
num_occupied_clusters=0  # number of occupied clusters or macrostates (union of Voronoi cells). it could be nonzero when
                         # restarting a simulation, otherwise 0.
num_balls_for_sc=100  # minimum number of Voronoi cells present to start performing spectral clustering.
num_clusters=5  # number of clusters or macrostates (union of Voronoi cells )for k-means part of spectral clustering.
num_walkers_for_sc=50  # number of walkers for each cluster or macrostate (union of Voronoi cells).
num_steps_for_sc=50  # number of steps needed to calculate the transition matrix for spectral clustering.

### the next four lines are needed if enhanced_sampling flag = 3 ###
initial_step_num_for_eq=1000  # starting CAS simulation step number for reweighting.
num_steps_for_eq=50  # number of steps needed to calculate the transition matrix for reweighting.
eq_frequency=2  # how many times to perform reweighting.
num_steps_in_bw=10  # number of steps in between reweightings.
