###### MVP Revision Tutorial ######
Defaults:
Lattice dimensions: 50 x 50
Initial state: 0.5 +/- 0.1
Number of Sweeps: 10000
Random seed: varies
dx = 1
dt = 0.1

Usage: python Ising.py <Tags>
Optional tags:
-x <value>        Lattice x dimension.
-y <value>        Lattice y dimension, defaults to match x.
-D <value>        Value of diffusion coefficient
-k <value>        Value of kappa
-S <value>        Value of sigma, determining source
-N <value>        Amount of noise to add to initial value
-dx <value>       Step size (space)
-dt <value>       Step size (time)
-i <values>       Initial value of the order parameter.
-rs <value>       Set random seed.
-T <values>       Number of sweeps to perform.
-u <value>        Update rate of the animation (sweeps between frames)
-p <value>        Time (ms) to pause between frames
-v <value>        Set a value for v0, part 5 only.
-H                Print this dialogue and exit.

Example with 40 x 25 lattice with 1000 sweeps, and noise 0.3
python Main.py -x 40 -y 25 -N 1000 -p [0.5,0.5,0.5]
