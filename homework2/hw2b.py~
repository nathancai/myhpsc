
"""
Demonstration module for quadratic interpolation.
Update this docstring to describe your code.
Modified by: ** Nathan Juswanto **
"""


import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import solve

def quad_interp(xi, yi):
    """
    Quadratic interpolation.  Compute the coefficients of the polynomial
    interpolating the points (xi[i],yi[i]) for i = 0,1,2.
    Returns c, an array containing the coefficients of
      p(x) = c[0] + c[1]*x + c[2]*x**2.

    """

    # check inputs and print error message if not valid:

    error_message = "xi and yi should have type numpy.ndarray"
    assert (type(xi) is np.ndarray) and (type(yi) is np.ndarray), error_message

    error_message = "xi and yi should have length 3"
    assert len(xi)==3 and len(yi)==3, error_message

    # Set up linear system to interpolate through data points:

    ### Fill in this part to compute c ###
    
    # It would be better to define A in terms of the xi points.
    # Doing this is part of the homework assignment.
    A = np.vstack([np.ones(3), xi, xi**2]).T
    b = yi

    # Solve the system:
    c = solve(A,b)

    return c

    
def cubic_interp(xi, yi):
    # check inputs and print error message if not valid:

    error_message = "xi and yi should have type numpy.ndarray"
    assert (type(xi) is np.ndarray) and (type(yi) is np.ndarray), error_message

    error_message = "xi and yi should have length 4"
    assert len(xi)==4 and len(yi)==4, error_message

    # Set up linear system to interpolate through data points:

    ### Fill in this part to compute c ###
    
    # It would be better to define A in terms of the xi points.
    # Doing this is part of the homework assignment.
    A = np.vstack([np.ones(4), xi, xi**2, xi**3]).T
    b = yi

    # Solve the system:
    c = solve(A,b)

    return c
 
 
def poly_interp(xi, yi):
    # check inputs and print error message if not valid:

    error_message = "xi and yi should have type numpy.ndarray"
    assert (type(xi) is np.ndarray) and (type(yi) is np.ndarray), error_message
   
    error_message = "xi and yi should have the same length"
    assert len(xi)==len(yi), error_message

    # Set up linear system to interpolate through data points:

    ### Fill in this part to compute c ###
    
    # It would be better to define A in terms of the xi points.
    # Doing this is part of the homework assignment.
    n = len(xi)
    A = np.ones(n)
    for j in range(1, n, 1):
        A = np.vstack([A, xi**j])
    A = np.transpose(A)
    b = yi

    # Solve the system:
    c = solve(A,b)

    return c

    
def plot_quad(xi, yi):
    # Plot the resulting polynomial:
    
    c = quad_interp(xi, yi)
    
    # points to evaluate polynomial
    x = np.linspace(xi.min() - 1,  xi.max() + 1, 1000)                           
    y = c[0] + c[1]*x + c[2]*x**2

    plt.figure(1)       # open plot figure window
    plt.clf()           # clear figure
    plt.plot(x,y,'b-')  # connect points with a blue line

    # Add data points  (polynomial should go through these points!)
    plt.plot(xi,yi,'ro')   # plot as red circles
    plt.ylim(-2, 8)         # set limits in y for plot

    plt.title("Data points and interpolating polynomial")

    plt.savefig('quadratic.png')   # save figure as .png file

    
def plot_cubic(xi, yi):
    # Plot the resulting polynomial:
    
    c = cubic_interp(xi, yi)
    
    # points to evaluate polynomial
    x = np.linspace(xi.min() - 1,  xi.max() + 1, 1000)                           
    y = c[0] + c[1]*x + c[2]*x**2 + c[3]*x**3

    plt.figure(1)       # open plot figure window
    plt.clf()           # clear figure
    plt.plot(x,y,'b-')  # connect points with a blue line

    # Add data points  (polynomial should go through these points!)
    plt.plot(xi,yi,'ro')   # plot as red circles
    plt.ylim(-3, 8)         # set limits in y for plot

    plt.title("Data points and interpolating polynomial")

    plt.savefig('cubic.png')   # save figure as .png file
    

def plot_poly(xi, yi):
    # Plot the resulting polynomial:
    
    c = poly_interp(xi, yi)
    
    # points to evaluate polynomial
    n = len(xi)    
    x = np.linspace(xi.min() - 1,  xi.max() + 1, 1000)
    y = c[n-1]
    for j in range(n-1, 0, -1):
        y = y*x + c[j-1]
    #y = c[0] + c[1]*x + c[2]*x**2 + c[3]*x**3

    plt.figure(1)       # open plot figure window
    plt.clf()           # clear figure
    plt.plot(x,y,'b-')  # connect points with a blue line

    # Add data points  (polynomial should go through these points!)
    plt.plot(xi,yi,'ro')   # plot as red circles
    plt.ylim(-3, 8)         # set limits in y for plot

    plt.title("Data points and interpolating polynomial")

    plt.savefig('poly.png')   # save figure as .png file
    
    
def test_quad1():
    """
    Test code, no return value or exception if test runs properly.
    """
    xi = np.array([-1.,  0.,  2.])
    yi = np.array([ 1., -1.,  7.])
    c = quad_interp(xi, yi)
    c_true = np.array([-1.,  0.,  2.])
    print "c =      ", c
    print "c_true = ", c_true
    # test that all elements have small error:
    assert np.allclose(c, c_true), \
        "Incorrect result, c = %s, Expected: c = %s" % (c,c_true)
        
    plot_quad(xi, yi)

    
def test_quad2():
    """
    Test code, no return value or exception if test runs properly.
    """
    xi = np.array([-1.,  0.,  2.])
    yi = np.array([ 1.,  5., -1.5])
            
    plot_quad(xi, yi)
        
        
def test_cubic1():
    """
    Test code, no return value or exception if test runs properly.
    """
    xi = np.array([-1.,  0., 1., 2.])
    yi = np.array([ 1., -1., 5., 3.])
  
    plot_cubic(xi, yi)
    
    
def test_poly1():
    """
    Test code, no return value or exception if test runs properly.
    """
    xi = np.array([-1., 0.,   1., 2.])
    yi = np.array([ 1., 4.25, 2., 7.])
  
    plot_poly(xi, yi)
    
    
def test_poly2():
    """
    Test code, no return value or exception if test runs properly.
    """
    xi = np.array([-1.,  0., 1., 2.5, 4.])
    yi = np.array([ 1., -1., 5., 2.,  3.5])
  
    plot_poly(xi, yi)    
    
    
if __name__=="__main__":
    # "main program"
    # the code below is executed only if the module is executed at the command line,
    #    $ python demo2.py
    # or run from within Python, e.g. in IPython with
    #    In[ ]:  run demo2
    # not if the module is imported.
    print "Running test..."
    
    # Homework 2.2
    test_quad1()
    
    # Homework 2.3 & 2.4
    test_quad2()
    
    # Homework 2.6
    test_cubic1()
    
    # Homework 2.7
    test_poly1()
    test_poly2()

