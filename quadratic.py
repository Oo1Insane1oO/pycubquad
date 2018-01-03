# given two points x0 and x1 with corresponding function values f0 and f1 along
# with derivatives g0 and g1 at these points, find point which minimizes the
# quadratic interpolation through f0, f1 and g0

import sympy as sp

def quadraticSolver():
    """ 
    wrapper function which returns the the solutions to the equation for
    minimization and the symbols used in order x0, x1, f0, f1, g0, a, b, c
    """
    # set symbolic variables for f(x) = ax^2+bx+c, with g=2ax+b being the
    # derivative of f with respect to x, and interpolation points x0 and x1
    a, b, c = sp.symbols('a b c', real=True)
    x0, x1 = sp.symbols('x_0 x_1', real=True, positive=True)
    f0, f1, g0 = sp.symbols('f_0 f_1 g_0', real=True)

    # define linear system as a matrix with the columns corresponding to a, b and c
    # respectively
    A = sp.Matrix((
        [x0**2, x0, 1, f0], 
        [x1**2, x1, 1, f1], 
        [2*x0, x0, 0, g0]))

    # split system in lhs matrix and rhs vector (Ay=r, with A the matrix, y the
    # vector y=(a,b,c) and r the rhs vector r=(f0,f1,g0)) and solve
    solution = tuple(sp.linsolve((A[:, :-1], A[:,-1]), a,b,c))[0]

    # split nominator and denominator in each coefficient (denominator is the same
    # for all coefficients except c)
    nominators = []
    denominators = []
    for s in solution:
        nom, denom = sp.fraction(sp.factor(sp.powsimp(s)))
        denominators.append(denom)
        nominators.append(sp.simplify(nom))

    # update symbolic variables (ignore denominator and c)
    a, b, c = nominators

    # collect common terms in g0 then x0
    a = sp.collect(a, (g0,x0))
    b = sp.collect(b, (g0,x0))
    c = sp.collect(b, (g0,x0))

    # split solution in numerator and denominator and restructure nominator
    snom, sdenom = sp.fraction(-b/(2*a))
    snomArgs = tuple(snom.args)
    snom = sp.factor(snomArgs[0]) + snomArgs[1]

    return (sp.simplify(snom / sdenom), x0, x1, f0, f1, g0, a/denominators[0],
            b/denominators[1], c/denominators[2])
# and function quadraticSolver

if __name__ == "__main__":
    """ print and plot solution """
    sol, x0, x1, f0, f1, g0, a, b, c = quadraticSolver()
    sp.pprint(sol)

    # make a cubic polynomial and substitute symbols with numerical values
    c3, c2, c1, c0 = 1.0, 1.0, 1./2, 1./3
    E = lambda x: c0*x**3 + c1*x**2 + c2*x + c3
    dE = lambda x: 3*c0*x**2 + 2*c1*x + c2
    x0Val = -1.0
    x1Val = 1.0
    f0Val = E(x0Val)
    f1Val = E(x1Val)
    g0Val = dE(x0Val)
    valDict = {x0:x0Val, x1:x1Val, f0:f0Val, f1:f1Val, g0:g0Val}
    f = lambda C: a.subs(valDict) * C**2 + b.subs(valDict) * C + \
            c.subs(valDict)
    sol = sol.subs(valDict)

    import numpy as np
    import matplotlib.pyplot as plt
    x = np.linspace(-10,10,100)
    plt.plot(sol, E(sol), 'bo', label="Extremal of quadratic fit, $%.3f$" %
            sol)
    plt.plot(x0Val, f0Val, 'ro', label="Interpolation point"
            "$(x0,f0)=(%.2f,%.2f)$" % (x0Val, f0Val))
    plt.plot(x1Val, f1Val, 'go', label="Interpolation point"
            "$(x1,f1)=(%.2f,%.2f)$" % (x1Val, f1Val))
    plt.plot(x, f(x), 'c-', label="Quadratic fit")
    plt.plot(x, E(x), 'y-', label="Cubic polynomial")
    plt.legend(loc='best')
    plt.show()
# end ifmain
