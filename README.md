The regular weighted ensemble simulation Python code was originally developed by H. Lee, Stanford. The bash script for 
multiple simultaneous simulations was originally developed by S. Ahn, Stanford. The resampling method implemented in the 
code is from E. Darve, Stanford and J. Izaguirre, University of Notre Dame. The enhanced sampling weighted ensemble 
simulation Python code was developed on top of the regular weighted ensemble code by S. Ahn, Stanford. Improvements and
additions have been made by J. Birgmeier, Stanford.

Please see the TODO sections in we_main.py, we_execute.sh, we_simulations.sh to get started (e.g. edit main directories, 
etc.). Then edit the we_initial_values.txt and we_parameters.py for your simulation.

Note that if more than one initial condition is provided, then write the initial values in order in the same file 
(we_initial_values.txt) on the next line, i.e. if I have two initial conditions x = 1.0, y = 2.0, and x = 3.0, y = 4.0, 
then my we_initial_values.txt would look like the following:
1.0
2.0
3.0
4.0

we_execute.sh is the script to submit to your cluster to run the weighted ensemble simulation.

Please report any bugs and questions/comments to sahn1@stanford.edu.
