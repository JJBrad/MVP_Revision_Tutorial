import numpy as np
import p5Lattice as lat
import interactive as interact
import matplotlib.pyplot as pyplot
import sys
from scipy.optimize import curve_fit

# Default values:
params = {"X Dimension":50,
          "Y Dimension":-1,
          "phi0" : 0.5,
          "kappa" : 0.01,
          "D" : 0.1,
          "noise" : 0.1,
          "sigma" : 10.,
          "dx" : 1.,
          "dt" : 0.1,
          "Seed" : None,
          "tMax" : 100000,
          "v0" : 0.1,
          "Interval" : 0, # Pause between frames
          "Rate" : 100 # Number of sweeps per update frame
          }

# Get input from command line
args = sys.argv[1:]
interact.readArgs(args, params)
np.random.seed(params["Seed"]) #None is default, changes each run.

lattice = lat.lattice(params["X Dimension"], 
                      params["Y Dimension"], 
                      params["D"], 
                      params["kappa"], 
                      params["sigma"],
                      params["dx"],
                      params["dt"],
                      params["phi0"], 
                      params["noise"],
                      params["v0"])

lattice.display(tMax=params["tMax"], rate=params["Rate"], interval=params["Interval"])

print("#"*40 + "\nNote: If animation was exited manually then an error may appear above.\nDisregard this error.\n" + "#"*40)
