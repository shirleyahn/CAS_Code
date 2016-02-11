import numpy as np


def check_state_function(coordinates):
    if -1.0 <= coordinates[0] <= 0.0 and -1.0 <= coordinates[1] <= 1.0 and \
       np.sqrt((coordinates[0]+1.0)**2 + coordinates[1]**2) <= 0.4:
        return 0
    elif 0.0 <= coordinates[0] <= 1.0 and -1.0 <= coordinates[1] <= 1.0 and \
            np.sqrt((coordinates[0]-1.0)**2 + coordinates[1]**2) <= 0.4:
        return 1
    else:
        return -1
