# $MYHPSC/homework3/newton.py
# Nathan Juswanto

def solve(fvals, x0, debug=False):
    k = 0
    maxiter = 20
    
    if debug:
         print "Initial guess: x = %22.15e" % x0
    
    while abs(fvals(x0)[0]) > 1e-14:
	 x0 -= fvals(x0)[0]/fvals(x0)[1]
	 k += 1
         if debug:
              print "After %s iterations, x = %22.15e" % (k, x0)
         if k == maxiter:
	      return x0, k
	 
    return x0, k

# $UWHPSC/codes/homework3/test_code.py 
# To include in newton.py

def fvals_sqrt(x):
    """
    Return f(x) and f'(x) for applying Newton to find a square root.
    """
    f = x**2 - 4.
    fp = 2.*x
    return f, fp

def test1(debug_solve=False):
    """
    Test Newton iteration for the square root with different initial
    conditions.
    """
    from numpy import sqrt
    for x0 in [1., 2., 100.]:
        print " "  # blank line
        x,iters = solve(fvals_sqrt, x0, debug=debug_solve)
        print "solve returns x = %22.15e after %i iterations " % (x,iters)
        fx,fpx = fvals_sqrt(x)
        print "the value of f(x) is %22.15e" % fx
        assert abs(x-2.) < 1e-14, "*** Unexpected result: x = %22.15e"  % x
        
if __name__=="__main__":
    test1()