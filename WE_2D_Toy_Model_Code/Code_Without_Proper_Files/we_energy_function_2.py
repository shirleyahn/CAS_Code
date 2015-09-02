import numpy as np


def energy_function(x, y):
    return 4.0*(x**2 + y**2 - 1.0)**2*y**2 - np.exp(-4.0*((x - 1.0)**2 + y**2)) - np.exp(-4.0*((x + 1.0)**2 + y**2)) + \
           np.exp(8.0*(x - 1.5)) + np.exp(-8.0*(x + 1.5)) + np.exp(-4.0*(y + 0.25)) + 0.2*np.exp(-8.0*x**2)
