The regular weighted ensemble simulation Python code was originally developed by H. Lee, Stanford. The bash script for 
multiple simultaneous simulations was originally developed by S. Ahn, Stanford. The resampling method implemented in the 
code is from E. Darve, Stanford and J. Izaguirre, University of Notre Dame. The Concurrent Adaptive Sampling (CAS) 
algorithm Python code was developed on top of the regular weighted ensemble code by S. Ahn, Stanford. Outlier detection
in spectral clustering has been implemented by J. Birgmeier, Stanford.

Please see the TODO sections in main.py, execute.sh, simulations.sh to get started (e.g. edit main directories, 
etc.). Then edit the initial_values.txt and parameters.py for your simulation.

Note that if more than one initial condition is provided, then write the initial values in order in the same file 
(initial_values.txt) on the next line, i.e., if I have two initial conditions x = 1.0, y = 2.0, and x = 3.0, y = 4.0, 
then my initial_values.txt would look like the following:
1.0 2.0
3.0 4.0

execute.sh is the script to submit to your cluster to run the CAS algorithm.

Please report any bugs and questions/comments to sahn1@stanford.edu.
