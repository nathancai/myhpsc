# $MYHPSC/homework3/intersections.py
# Nathan Juswanto


import numpy as np
import matplotlib.pyplot as plt
from newton import solve

def fvals_sqrt(x):
    """
    Return f(x) and f'(x) for applying Newton to find a square root.
    """
    f = x*np.cos(np.pi*x)-1+.6*x**2
    fp = x*np.pi*-np.sin(np.pi*x)+np.cos(np.pi*x)+1.2*x
    return f, fp

def test1(debug_solve=True):
    """
    Test Newton iteration for the square root with different initial
    conditions.
    """
          
    # points to evaluate polynomial
    x = np.linspace(-5., 5, 1000)                           
    g1 = x*np.cos(np.pi*x)
    g2 = 1-.6*x**2

    plt.figure(1)       # open plot figure window
    plt.clf()           # clear figure
    ax = plt.subplot(111)
    ax.plot(x, g1, 'b-', label='g1')  # connect points with a blue line
    ax.plot(x, g2, 'r-', label='g2')
    ax.legend(loc='lower center')

    for x0 in [-2., -1.5, -1., 1.5]:
        print " "  # blank line
        x, iters = solve(fvals_sqrt, x0, debug=debug_solve)
        print "Initial guess: x = %22.15e" % x0
        print "solve returns x = %22.15e after %i iterations " % (x,iters)
        fx = 1-.6*x**2
        print "the value of f(x) is %22.15e" % fx
        # Add data points  (polynomial should go through these points!)
        ax.plot(x, fx, 'ko')   # plot as red circles

    plt.title("Intersections")

    plt.savefig('intersections.png')   # save figure as .png file
    
if __name__=="__main__":
    test1()