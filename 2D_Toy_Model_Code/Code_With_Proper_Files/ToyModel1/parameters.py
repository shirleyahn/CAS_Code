main_directory='/scratch/users/sahn1/2D_Toy_Model'  # main directory of where the CAS simulation will take place,
                                                    # i.e., where all of the files are.

balls_flag=0  # 0: create new macrostates at each step. 1: keep created macrostates (only allowed if enhanced sampling
              # is off or enhanced_sampling_flag = 0).
flux_flag=1  # 0: off. 1: on. fluxes between pre-defined states will be calculated. the walker's state is determined by
             # check_state_function.py.
num_states=2  # number of pre-defined states for flux calculation. only needed if flux_flag = 1, otherwise 0.
enhanced_sampling_flag=0  # 0: off. 1: binning walkers if the walkers have some property less or greater than threshold.
                          # 2: spectral clustering. 3: reweighting with equilibrium weights.

num_balls_limit=1000  # maximum number of macrostates in the simulation.
radius=[0.1, 0.1]  # radii of the Voronoi cells.
num_walkers=100  # number of walkers per macrostate.
grid_dimensions=[-1.0, 1.0, -1.0, 1.0]  # since this is for a 2D toy model, this is equal to
                                        # [x_lower_bound, x_upper_bound, y_lower_bound, y_upper_bound].
angle_cvs=[0, 0]  # list of 0's and/or 1's. 0: if the collective variable is not an angle.
                  # 1: if the collective variable is an angle that ranges from -180 deg to 180 deg.

max_num_steps=100000  # maximum number of the CAS simulation steps.
num_occupied_balls=1  # starting number of initial states for the CAS simulation.

m_steps_per_step=1  # how many times the Metropolis algorithm should be executed per simulation step.
step_size=0.05  # how large the step size should be for each walker.
beta=10.0  # inverse temperature of the system.
pbc=0  # 0: off. 1: periodic boundary conditions on.

### the next four lines are needed if enhanced_sampling_flag = 1 ###
less_or_greater_flag=0  # 0: criteria for binning walkers if the walkers have some property LESS than the threshold.
                        # 1: criteria for binning walkers if the walkers have some property GREATER than the threshold.
threshold_values=[-0.5, 0.5]  # list of threshold values. if at least one property of the walker has a value less or
                              # greater than the threshold value, then it is binned coarsely.
properties_to_keep_track=[0, 0]  # list of properties of the walker that are compared against threshold values. this can
                                 # be weight and/or some collective variables and multiple duplicates are allowed.
                                 # if one of them is weight, then type -1. otherwise type the indices of the collective
                                 # variables, e.g., if there are 3 collective variables and you would like to keep track
                                 # of the last one twice for different threshold values, then type [2, 2] (index starts
                                 # from 0). if there is more than one property to keep track of, then type them sequentially.
coarse_radius=[0.5, 0.5]  # radius for coarse binning.

### the next four lines are needed if enhanced_sampling flag = 2 ###
num_balls_for_sc=100  # minimum number of Voronoi cells present to start performing spectral clustering.
num_clusters=5  # number of clusters or macrostates (union of Voronoi cells )for k-means part of spectral clustering.
num_walkers_for_sc=50  # number of walkers for each cluster or macrostate (union of Voronoi cells).
num_steps_for_sc=50  # number of steps needed to calculate the transition matrix for spectral clustering.

### the next four lines are needed if enhanced_sampling flag = 3 ###
initial_step_num_for_eq=1000  # starting CAS simulation step number for reweighting.
num_steps_for_eq=50  # number of steps needed to calculate the transition matrix for reweighting.
eq_frequency=2  # how many times to perform reweighting.
num_steps_in_bw=10  # number of steps in between reweightings.
