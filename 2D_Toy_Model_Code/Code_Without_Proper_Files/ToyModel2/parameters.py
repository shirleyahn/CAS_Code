main_directory='/scratch/users/sahn1/2D_Toy_Model'

balls_flag=0  # 0: create new balls at each step. 1: keep created balls.
rate_flag=1  # 0: off. 1: on. rates/fluxes between pre-defined states will be calculated. the walker's state is
             # determined by check_state_function.py.
num_states=2  # number of pre-defined states for rate/flux calculation. only needed if rate_flag = 1, otherwise 0.
num_pathways=2  # number of distinct pathways to keep track of. only needed if rate_flag = 1, otherwise 0.
enhanced_sampling_flag=0  # 0: off. 1: binning walkers if the walkers have some property less or greater than threshold.
                          # 2: spectral clustering.

num_balls_limit=1000  # parameter needed in case the calculated max_num_balls is greater than the limit.
radius=0.1
num_walkers=100
num_cvs=2  # number of collective variables (num_cvs) should be fixed to be 2 for this case.
grid_dimensions=[-1.5, 1.5, -0.5, 1.25]  # since num_cvs = 2, then type x_lower_bound x_upper_bound
                                         # y_lower_bound y_upper_bound
angle_cvs=[0, 0]  # 0: if the cv is not an angle. 1: if the cv is an angle.

max_num_steps=200  # maximum number of steps for the simulation.
num_occupied_balls=2

m_steps_per_step=5  # how many times the metropolis algorithm should be executed per step
step_size=0.05  # how large the step size should be for each walker
beta=10.0  # inverse temperature
pbc=0  # 0: off. 1: periodic boundary conditions on.

### for the next four lines, if enhanced_sampling_flag = 1 ###
less_or_greater_flag=0  # 0: criteria for binning walkers is if the walkers have some property LESS than the threshold.
                        # 1: criteria for binning walkers is if the walkers have some property GREATER than the
                        # threshold.
static_threshold_flag=1  # 0: off, then the lowest (less_or_greater_flag = 0) or highest (less_or_greater_flag = 1)
                         # current value is set as the threshold for the next step. 1: on, initial threshold is kept
                         # throughout the simulation.
threshold_values=[1.0e-100]  # if some properties of the walker have values less or greater than the threshold values,
                             # then it is binned to the nearest existing ball.
properties_to_keep_track=[-1]  # properties of the walker that are compared against the threshold values. this can be
                               # weight and/or some cv(s). if one of them is weight, then type -1.  otherwise type the
                               # indices of the cv, e.g. if there are 3 cvs and you would like to keep track of the
                               # last one, type 2 (index starts from 0). if more than one property is kept track of,
                               # then type them sequentially.

### for the next three lines, if enhanced_sampling flag == 2 ###
num_balls_for_sc=1000  # minimum number of balls present to perform spectral clustering for that step
num_clusters=5  # number of clusters for k-means part of spectral clustering
num_walkers_for_sc=1000  # number of walkers for each macrostate, usually set equal to the avg number of walkers per
                         # macrostate, which is (num_balls_for_sc/num_clusters)*num_walkers
