import numpy as np
import matplotlib.pyplot as pyplot
from matplotlib import colors
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit
import os
import VectorCalc as vc
from PIL import Image

""" TODO:
write chi^2 function.
Consider transposing phi in animation
"""

def Exp(x, *p):
    A, r_c = p
    return(A*np.exp(-x/r_c))

def SquareExp(x, *p):
    A, r_c = p
    return(A*np.exp(-np.square(x)/r_c))

def fit(yList, xList, fitFunction, initGuess, minX=-np.inf, maxX=np.inf):
    """
    Function to fit a function (fitFunction) to a set of data, only including data over
    a range of y.
    """
    xVals = np.copy(xList)
    yVals = np.copy(yList)
    # Select only data of interest.
    yVals = yVals[(xVals < maxX) & (xVals > minX)]
    xVals = xVals[(xVals < maxX) & (xVals > minX)]
    # Perform fit
    params, var = curve_fit(fitFunction, xVals, yVals, p0=initGuess)
    # Calculate uncertainty in parameters
    err = np.sqrt(np.diag(var))
    return(params, err)

class lattice(object):
    def __init__(self, xDim=100, yDim=0, D=1., kappa=0.5, sigma=2, dx=1., dt=1., initVal=0.5, noise=0.1, v0=0.01):
        """
        Constructor for the lattice object. Defaults to square lattice
        :param xDim: The x dimension of the lattice. Defaults to 50
        :param yDim: The y dimension of the lattice (optional), defaults to square lattice
        :param D: Diffusion coefficient
        :param kappa: kappa value
        :param initVal: Average value of initial lattice
        :param noise: The noise in the initial order paramter (Uniform)
        :param dx: The step size in space
        :param dt: The step size in time
        """
        self.D = D
        self.kappa = kappa
        self.dx = dx
        self.dt = dt
        self.sigma = sigma
        self.v0 = v0
        # Initialise lattice
        self.n = 0 # Number of sweeps performed.
        self.xDim = xDim
        if yDim > 0:
            self.yDim = yDim
        else:
            self.yDim = xDim
        self.size = self.xDim*self.yDim
                
        self.phi0 = initVal
        self.phi = np.random.uniform(-1.*noise, 1.*noise, size=(self.xDim, self.yDim)) + self.phi0
        
        self.mid = (float(self.xDim)/2., float(self.yDim)/2.)
        
        # Initialise source
        self.rho = np.zeros(shape=(self.xDim, self.yDim))
        for i in range(0, self.xDim):
            for j in range(0, self.yDim):
                rSquare = float(i-self.mid[0])**2. + float(j-self.mid[1])**2.
                self.rho[i, j] = np.exp(-rSquare/self.sigma**2.)

        # Initialise flow
        self.v = np.zeros(shape=(self.xDim, self.yDim))
        for i in range(0, self.xDim):
            for j in range(0, self.yDim):
                self.v[i, j] = - self.v0 * np.sin(2.*np.pi*j/self.yDim)

        # Initialise lists for analysis
        self.nList = []
        self.avPhi = []
        self.nList.append(self.n)
        self.avPhi.append(np.average(self.phi))

        print(self)
        print(self.dt)

            
    def __str__(self):
        """
        Returns string for printing key details of object.
        """
        return("Array has shape {}. Average phi is {:.3e}.".format(self.phi.shape, np.average(self.phi)))
    
    def next(self):
        """
        Perform one sweep.
        """
        gradPhi = vc.grad2D(self.phi, self.dx)[0]
        self.phi = self.phi + self.D*self.dt*vc.lap2D(self.phi, self.dx) + self.dt*self.rho - self.kappa*self.dt*self.phi - self.v * gradPhi
        self.n += 1
        self.nList.append(self.n)
        self.avPhi.append(np.average(self.phi))

    def animate(self, f, tMax, rate):
        for _ in range(0, rate):
            if self.n <= tMax:
                self.next()
        if f % 100 == 0:
            print("{} steps completed".format(f * rate))
            #print("Max: {}, min: {}".format(np.max(self.phi), np.min(self.phi)))
        if self.n == tMax:
            print("Animation complete.")
            self.n += rate
        im = axis.imshow(self.phi, cmap="YlOrRd", interpolation="nearest", vmin=-1, vmax=76) #"seismic"
        pyplot.clim(-1, 76)
        #print(self.n)
        return([im])
        
    def display(self, tMax=20000, rate=10, interval=0):
        global axis, fig, cax
        fig, axis = pyplot.subplots()
        for s in axis.spines: axis.spines[s].set_color("k")
        axis.spines['left'].set_position(('outward', 1))
        axis.spines['top'].set_position(('outward', 0.5))
        # Position colourbar
        im = pyplot.imshow(self.phi, cmap="YlOrRd", interpolation="nearest", vmin=-1, vmax=76) #"bwr"
        pyplot.clim(-1, 76)
        cbar = pyplot.colorbar(im)
        anim = FuncAnimation(fig, self.animate, tMax/rate + 1, interval=interval, blit=True, repeat=False, fargs=(tMax,rate))
        pyplot.show()
        print("Animation finished. Please close the animation window.")

    def run(self, tMax=1000):
        for i in range(0, tMax):
            self.next()
            
    def analyse(self):
        """
        Function to plot all required graphs after a run and do any necessary calculations
        """
        # Print values at final time.
        print("At step {}, maximum phi is {:.2f}, average is {:.2f}".format(self.n, np.max(self.phi), np.average(self.phi)))
        # Plot average phi as a function of step number
        pyplot.plot(self.nList, self.avPhi, "k-")
        pyplot.xlabel("n")
        pyplot.ylabel(r"$\left< \phi \right>$")
        pyplot.show()
        # Plot phi as a function of distance from the centre.
        # Define lists:
        rList = []
        phiList = []
        # Loop through array
        for i in range(0, self.xDim):
            for j in range(0, self.yDim):
                r = np.sqrt((float(i) - self.mid[0])**2. + (float(j) - self.mid[1])**2.) # Calculate distance
                rList.append(r) # append to lists
                phiList.append(self.phi[i, j])
        if False:
            # Plot (linear plot)
            pyplot.plot(rList, phiList, "kx")
            pyplot.xlabel("R")
            pyplot.ylabel(r"$\phi$")
            pyplot.show()
            # Plot (log-linear plot)
            pyplot.plot(rList, phiList, "kx")
            pyplot.xlabel("R")
            pyplot.ylabel(r"$\phi$")
            pyplot.yscale("log")
            pyplot.show()
            # Plot (log-log plot)
            pyplot.plot(rList, phiList, "kx")
            pyplot.xlabel("R")
            pyplot.ylabel(r"$\phi$")
            pyplot.yscale("log")
            pyplot.xscale("log")
            pyplot.show()
            # Plot (log phi vs R^2 plot)
            pyplot.plot(np.square(rList), phiList, "kx")
            pyplot.xlabel(r"$R^2$")
            pyplot.ylabel(r"$\phi$")
            pyplot.yscale("log")
            pyplot.show()
        # Perform fit to exponential decay for R between 10 and 25.
        expVals, expErr = fit(phiList, rList, Exp, [75, 7], maxX = 25, minX = 10.)
        xFit = np.linspace(5, 30)
        yFit = Exp(xFit, *expVals)
        pyplot.plot(rList, phiList, "kx", label="Data")
        pyplot.plot(xFit, yFit, "r--", label=r"$A \exp(-r/r_c)$")
        pyplot.xlabel("R")
        pyplot.ylabel(r"$\phi$")
        pyplot.title("A = {:.3f}, r_c = {:.3f}".format(*expVals))
        print("Error in A: {:.3f}".format(expErr[0]))
        print("Error in r_c: {:.3f}".format(expErr[1]))
        pyplot.show()
        # Perform fit to exponential decay for R^2 up to R^2 = 900 (R = 30)
        exp2Vals, exp2Err = fit(phiList, rList, SquareExp, [75, 10], maxX = 30)
        xFit = np.linspace(0, 30)
        yFit = SquareExp(xFit, *exp2Vals)
        pyplot.plot(rList, phiList, "kx", label="Data")
        pyplot.plot(xFit, yFit, "r--", label=r"$A \exp(-r^2/r_c)$")
        pyplot.xlabel("R")
        pyplot.ylabel(r"$\phi$")
        pyplot.legend()
        pyplot.title("A = {:.3f}, r_c = {:.3f}".format(*exp2Vals))
        print("Error in A: {:.3f}".format(exp2Err[0]))
        print("Error in r_c: {:.3f}".format(exp2Err[1]))
        pyplot.show()
