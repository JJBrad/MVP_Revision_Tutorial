"""
VectorCalc.py
A useful class for simple discretised vector calculus with numpy arrays of different dimensions.
Author: J. Bradley
Updated: 01/05/2019
"""

import numpy as np

"""Notes:
The numpy roll function np.roll(arr, N, i) moves all entries in arr N spaces along index i, with periodic
boundary conditions. Moving the whole array is much faster than working componentwise since numpy is
a compiled C++ module.
e.g. lat[x+1] = np.roll(lat, -1, 0) moves all elements one element backwards along 0 axis.

Centered difference rules are used for spatial gradients to avoid biasing specific directions.
"""

def lap3D(lat, dx=1.):
    """
    Returns the (discretised) laplacian at each point in a 3d lattice with periodic BCs
    :param dx:   The size of finite steps in space
    :param lat:  The lattice of values in space for gradient calculation
    :return lap: Numpy array with values of the laplacian at each lattice point.
    """
    lap = np.roll(lat, 1, 0) + np.roll(lat, -1, 0) + \
          np.roll(lat, 1, 1) + np.roll(lat, -1, 1) + \
          np.roll(lat, 1, 2) + np.roll(lat, -1, 2) - \
          6. * lat
    lap = 1./dx**2. * lap
    return(lap)
    

def lap2D(lat, dx=1.):
    """
    Returns the (discretised) laplacian at each point in a 2d lattice with periodic BCs
    :param dx:   The size of finite steps in space
    :param lat:  The lattice of values in space for gradient calculation
    :return lap: Numpy array with values of the laplacian at each lattice point.
    """
    lap = np.roll(lat, 1, 0) + np.roll(lat, -1, 0) + \
          np.roll(lat, 1, 1) + np.roll(lat, -1, 1) - \
          4. * lat
    lap = 1./dx**2. * lap
    return(lap)
    
def lap1D(lat, dx=1.):
    """
    Returns the (discretised) laplacian at each point in a 1d lattice with periodic BCs
    :param dx:   The size of finite steps in space
    :param lat:  The lattice of values in space for gradient calculation
    :return lap: Numpy array with values of the laplacian at each lattice point.
    """
    lap = np.roll(lat, 1, 0) + np.roll(lat, -1, 0) - \
          2. * lat
    lap = 1./dx**2. * lap
    return(lap)

def gradSquared3D(lat, dx):
    """
    Returns the square of the gradient of the lattice at each point in space, using periodic BCs.
    :param dx:   The size of finite steps in space
    :param lat:  The lattice of values in space for gradient calculation
    :return res: Numpy array with values of the gradient squared at each lattice point.
    """
    # Gradient squared of f is sum of (df/dx_i)^2
    res = np.square(np.roll(lat, 1, 0) - np.roll(lat, -1, 0)) + \
          np.square(np.roll(lat, 1, 1) - np.roll(lat, -1, 1)) + \
          np.square(np.roll(lat, 1, 2) - np.roll(lat, -1, 2))
    res = res/(4.*dx**2.)
    return(res)
    
def gradSquared2D(lat, dx):
    """
    Returns the square of the gradient of the lattice at each point in space, using periodic BCs.
    :param dx:   The size of finite steps in space
    :param lat:  The lattice of values in space for gradient calculation
    :return res: Numpy array with values of the gradient squared at each lattice point.
    """
    # Gradient squared of f is sum of (df/dx_i)^2
    res = np.square(np.roll(lat, 1, 0) - np.roll(lat, -1, 0)) + \
          np.square(np.roll(lat, 1, 1) - np.roll(lat, -1, 1))
    res = res/(4.*dx**2.)
    return(res)
    
def gradSquared1D(lat, dx):
    """
    Returns the square of the gradient of the lattice at each point in space, using periodic BCs.
    :param dx:   The size of finite steps in space
    :param lat:  The lattice of values in space for gradient calculation
    :return res: Numpy array with values of the gradient squared at each lattice point.
    """
    # Gradient squared of f is sum of (df/dx_i)^2
    res = np.square(np.roll(lat, 1, 0) - np.roll(lat, -1, 0))
    res = res/(4.*dx**2.)
    return(res)

def grad3D(lat, dx):
    """
    Returns the components of the gradient of the lattice at each point in space, using periodic BCs.
    :param dx:   The size of finite steps in space
    :param lat:  The lattice of values in space for gradient calculation
    :return x:   Numpy array with values of the x component of gradient at each lattice point.
    :return y:   Numpy array with values of the x component of gradient at each lattice point.
    :return z:   Numpy array with values of the x component of gradient at each lattice point.
    """
    x = (np.roll(lat, -1, 0) - np.roll(lat, 1, 0))/(2.*dx)
    y = (np.roll(lat, -1, 1) - np.roll(lat, 1, 1))/(2.*dx)
    z = (np.roll(lat, -1, 2) - np.roll(lat, 1, 2))/(2.*dx)
    return((x, y, z))
    
def grad2D(lat, dx):
    """
    Returns the components of the gradient of the lattice at each point in space, using periodic BCs.
    :param dx:   The size of finite steps in space
    :param lat:  The lattice of values in space for gradient calculation
    :return x:   Numpy array with values of the x component of gradient at each lattice point.
    :return y:   Numpy array with values of the x component of gradient at each lattice point.
    """
    x = (np.roll(lat, -1, 0) - np.roll(lat, 1, 0))/(2.*dx)
    y = (np.roll(lat, -1, 1) - np.roll(lat, 1, 1))/(2.*dx)
    return((x, y))
    
def grad1D(lat, dx):
    """
    Returns the components of the gradient of the lattice at each point in space, using periodic BCs.
    :param dx:   The size of finite steps in space
    :param lat:  The lattice of values in space for gradient calculation
    :return x:   Numpy array with values of the x component of gradient at each lattice point.
    """
    x = (np.roll(lat, -1, 0) - np.roll(lat, 1, 0))/(2.*dx)
    return((x))
