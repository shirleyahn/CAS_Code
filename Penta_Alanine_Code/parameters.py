main_directory='/scratch/users/sahn1/Penta_Alanine'
initial_configuration_directory='/scratch/users/sahn1/Penta_Alanine/InitConfig'

simulation_flag=0  # 0: new simulation. 1: restarting simulation that didn't run all walkers.
                   # 2: restarting simulation that didn't finish post-processing.
                   # 3: restarting simulation that didn't finish binning.
                   # 4: restarting simulation that didn't finish resampling.
balls_flag=0  # 0: create new balls at each step. 1: keep created balls.
rate_flag=1  # 0: off. 1: on. rates/fluxes between pre-defined states will be calculated. the walker's state is
             # determined by check_state_function.py.
num_states=2  # number of pre-defined states for rate/flux calculation. only needed if rate_flag = 1, otherwise 0.
enhanced_sampling_flag=2  # 0: off. 1: binning walkers if the walkers have some property less or greater than threshold.
                          # 2: spectral clustering.

num_balls_limit=600  # parameter needed in case the calculated max_num_balls is greater than the limit.
radius=80.0  # radius can be changed in the middle of the simulation.
num_walkers=10  # num_walkers should be fixed.
num_cvs=6  # number of collective variables (num_cvs) should be fixed.
lower_bound=-180.0  # lower bound value for the collective variables.
upper_bound=180.0  # upper bound value for the collective variables.
angle_cvs=[1, 1, 1, 1, 1, 1]  # 0: if the cv is not an angle. 1: if the cv is an angle.

initial_step_num=0  # initial_step_num should be changed from 0 when restarting a simulation.
max_num_steps=50  # maximum number of steps for the simulation.
num_occupied_balls=1  # num_occupied_balls should be changed when restarting a simulation.
first_walker=0  # only needed if simulation_flag is not equal to 0, otherwise 0.
last_walker=0  # only needed if simulation_flag is not equal to 0, otherwise 0.

### for the next four lines, if enhanced_sampling_flag = 1 ###
less_or_greater_flag=1  # 0: criteria for binning walkers if the walkers have some property LESS than the threshold.
                        # 1: criteria for binning walkers if the walkers have some property GREATER than the threshold.
static_threshold_flag=0  # 0: off, then the lowest (less_or_greater_flag = 0) or highest (less_or_greater_flag = 1)
                         # value is set to be the threshold for the next step and all walkers are replaced with the
                         # "reference walker." 1: on, initial threshold is kept throughout the simulation.
threshold_values=[20.0]  # if some properties of the walker have values less or greater than the threshold values,
                         # then it is binned to the one designated "leftover" macrostate.
properties_to_keep_track=[0]  # properties of the walker that are compared against threshold values. this can be weight
                              # and/or some cv(s). if one of them is weight, then type -1. otherwise type the
                              # indices of the collective variables, e.g. if there are 3 cvs and you would like to keep
                              # track of the last one, then type 2 (index starts from 0). if there is more than one
                              # property to keep track of, then type them sequentially.

### for the next three lines, if enhanced_sampling flag == 2 ###
num_occupied_big_clusters=0  # num_occupied_big_clusters could be changed when restarting a simulation, otherwise 0.
num_occupied_small_clusters=0  # num_occupied_small_clusters could be changed when restarting a simulation, otherwise 0.
num_balls_for_sc=500  # minimum number of balls present to perform spectral clustering for that step.
num_clusters=5  # number of clusters for k-means part of spectral clustering.
num_walkers_for_sc=500  # number of walkers for each cluster.
num_steps_for_sc=10  # number of steps of calculating transition matrix before spectral clustering.
