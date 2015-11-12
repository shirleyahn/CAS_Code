main_directory='/scratch/users/jbirgmei/CS229/penta_alanine'
initial_configuration_directory='/scratch/users/jbirgmei/CS229/penta_alanine/StarterFolderGROMACS_WE'

simulation_flag=1  # 0: new simulation. 1: restarting simulation from middle of simulation. 2: restarting simulation from
                   # middle of binning. 3: restarting simulation from middle of resampling.
balls_flag=0  # 0: create new balls at each step. 1: keep created balls.
sorting_flag=1  # 0: sort walkers' weights in descending order (most probable walkers first). 1: sort walkers' weights in
              # ascending order (rare walkers first).
rate_flag=1  # 0: off. 1: on. rates/fluxes between pre-defined states  will be calculated. the walker's state is
           #  determined by we_check_state_function.py.
num_states=2  # number of pre-defined states for rate/flux calculation. only needed if rate_flag = 1, otherwise 1.
enhanced_sampling_flag=3  # 0: off. 1: sub-binning balls by standard deviation distance from center of ball. 2: binning
                        # walkers if the walkers have some property less or greater than threshold. 3: spectral
                        # clustering.

num_balls_limit=100000  # limit is set depending on the available memory. parameter needed in case the calculated max_num_balls
                 # is greater or too much smaller than the limit.
radius=80  # radius can be changed in the middle of the simulation.
num_walkers=10  # num_walkers should be fixed.
num_cvs=6  # number of collective variables (num_cvs) should be fixed.
lower_bound=-180.0  # lower bound value for the collective variables. set it to arbitrary value if collective variables have
             # different units. this is only used to calculate the volume of the free energy landscape.
upper_bound=180.0  # upper bound value for the collective variables.

initial_step_num=60  # initial_step_num should be changed from 0 when restarting a simulation.
max_num_steps=5  # maximum number of steps for the simulation.
num_occupied_balls=1118  # num_occupied_balls should be changed when restarting a simulation.
first_walker=0  # only needed if simulation_flag is not equal to 0, otherwise put 0.
last_walker=11179  # only needed if simulation_flag is not equal to 0, otherwise put 0.

### for the next four lines, if enhanced_sampling_flag = 2 ###
less_or_greater_flag=1  # 0: criteria for binning walkers is if the walkers have some property LESS than the threshold.
                      # 1: criteria for binning walkers is if the walkers have some property GREATER than the threshold.
static_threshold_flag=0  # 0: off, then the lowest (less_or_greater_flag = 0) or highest (less_or_greater_flag = 1)
                       # current value is set as the threshold for the next step. 1: on, initial threshold is kept
                       # throughout the simulation.
threshold_values=[20.0, 20.0, 20.0]  # if some properties of the walker have values less or greater than the threshold
                                     # values, then it is binned to the nearest existing ball.
properties_to_keep_track=[20, 21, 22]  # properties of the 
                          # walker that are compared against the threshold values. this can be weight
                          # and/or some collective variable(s). if one of them is weight, then type -1. otherwise type
                          # the indices of the collective variable, e.g. if there are 3 collective variables and you
                          # would like to keep track of the last one, type 2
                          # (index starts from 0). if more than one property is kept track of, then type them sequentially

### for the next four lines, if enhanced_sampling flag == 3 ###
num_balls_for_sc=1000  # minimum number of balls present to perform spectral clustering for that step
num_clusters=5  # number of clusters for k-means part of spectral clustering
num_walkers_for_sc=2000  # number of walkers for each macrostate, usually set equal to the avg number of walkers per
                    # macrostate, which is (num_balls_for_sc/num_clusters)*num_walkers
timestep=10.0  # length of timestep in ps
