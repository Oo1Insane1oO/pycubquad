# given two points x0 and x1 with corresponding function values f0 and f1 along
# with derivatives g0 and g1 at these points, find point which minimizes the
# quadratic interpolation through f0, f1 and g0

import sympy as sp
import sys

class Quad:
    """ class to keep symbols at head """
    def __init__(self):
        """ set symbols in class """
        self.x0, self.x1, self.f0, self.f1, self.g0, self.g1 =  \
                sp.symbols('x_0 x_1 f_0 f_1 g_0 g_1', real=True)
    # end __init__

    def splitter(self, solution):
        """ split numerator and denominator """
        # split nominator and denominator in each coefficient in solution
        # (assumed ordered a, b, c)
        numerators = []
        denominators = []
        for c in solution:
            num, denom = sp.fraction(sp.factor(sp.powsimp(c)))
            numerators.append(num)
            denominators.append(sp.factor(denom))
        # end forc

        return numerators, denominators
    # end function splitter

    def manipulatef0f1g0(self, solution):
        """ manipulate system (f0 = f(x0), f1=f(x1), g0=g(x0) """
        numerators, denominators = self.splitter(solution)

        numerators[0] = sp.collect(numerators[0], self.g0)
        numerators[1] = sp.simplify(sp.collect(numerators[1], (self.g0,
            self.x0)))
        num1Args = tuple(numerators[1].args)
        numerators[1] = sp.factor(num1Args[0]) + num1Args[1]
        numerators[2] = sp.collect(numerators[2], (self.g0, self.f0))
        num2Args = list(numerators[2].args)
        num20Args = tuple(num2Args[0].args)
        num2Args[0] = num20Args[0] * sp.factor(num20Args[1], self.x0)
        num22Args = tuple(num2Args[2].args)
        num2Args[2] = num22Args[0] * sp.factor(num22Args[1], (self.x0, self.x1))
        numerators[2] = sum(num2Args)
 
        return [numerator/denominator for numerator, denominator in
                zip(numerators,denominators)]
    # end function manipulatef0f1g0

    def manipulatef0f1g1(self, solution):
        """ manipulate system (f0 = f(x0), f1=f(x1), g1=g(x1) """
        numerators, denominators = self.splitter(solution)

        numerators[0] = sp.simplify(sp.collect(numerators[0], self.g1))
        numerators[1] = sp.simplify(sp.collect(numerators[1], (self.g1,self.x1)))
        num1Args = tuple(numerators[1].args)
        numerators[1] = sp.factor(num1Args[0]) + num1Args[1]
        numerators[2] = sp.collect(numerators[2], (self.g1, self.f1))
        num2Args = list(numerators[2].args)
        num20Args = tuple(num2Args[0].args)
        num2Args[0] = num20Args[0] * sp.factor(num20Args[1], self.x0)
        num22Args = tuple(num2Args[2].args)
        num2Args[2] = num22Args[0] * sp.factor(num22Args[1], (self.x0, self.x1))
        numerators[2] = sum(num2Args)
        num2Args = tuple(numerators[2].args)
        numerators[2] = num2Args[0] + sp.factor(num2Args[1], self.x0) + \
                num2Args[2]

        return [numerator/denominator for numerator, denominator in
                zip(numerators,denominators)]
    # end function manipulatef0f1g1

    def manipulatef0g0g1(self, solution):
        """ manipulate system (f0 = f(x0), g0=g(x0), g1=g(x1) """
        numerators, denominators = self.splitter(solution)

        numerators[2] = sp.collect(numerators[2], self.f0)
        num2Args = tuple(numerators[2].args)
        numerators[2] = sp.factor(num2Args[0]) + sp.factor(num2Args[1] +
                sp.factor(sum(num2Args[2:]), self.x0), self.x0)

        return [numerator/denominator for numerator, denominator in
                zip(numerators,denominators)]
    # end function manipulatef0g0g1

    def manipulatef1g0g1(self, solution):
        """ manipulate system (f1 = f(x1), g0=g(x0), g1=g(x1) """
        numerators, denominators = self.splitter(solution)

        numerators[2] = sp.simplify(sp.collect(numerators[2], self.f1))
        num2Args = tuple(numerators[2].args)
        numerators[2] = num2Args[2] + sp.factor(sum(num2Args[:2]) +
                num2Args[-1], self.x1)

        return [numerator/denominator for numerator, denominator in
                zip(numerators,denominators)]
    # end function manipulatef1g0g1

    def quadraticSolver(self):
        """ 
        wrapper function which returns the the solutions to the equation for
        minimization and the symbols used in order x0, x1, f0, f1, g0, g1, a,
        b, c
        """
        # set symbolic variables for f(x) = ax^2+bx+c, with g=2ax+b being the
        # derivative of f with respect to x, and interpolation points x0 and x1
        feq = lambda x,f: (x**2, x, 1, f)
        geq = lambda x,g: (2*x, 1, 0, g)

        # define linear system as a matrix with the columns corresponding to a,
        # b and c respectively for all possible system of equations (f0,f1,g0),
        # (f0,f1,g1), (f0,g0,g1), (f1,g0,g1) (in that order)
        M = (
                sp.Matrix((
                    feq(self.x0,self.f0), 
                    feq(self.x1,self.f1), 
                    geq(self.x0,self.g0))),
                sp.Matrix((
                    feq(self.x0,self.f0), 
                    feq(self.x1,self.f1), 
                    geq(self.x1,self.g1))),
                sp.Matrix((
                    feq(self.x0,self.f0), 
                    geq(self.x0,self.g0), 
                    geq(self.x1,self.g1))),
                sp.Matrix((
                    feq(self.x1,self.f1), 
                    geq(self.x0,self.g0), 
                    geq(self.x1,self.g1)))
                )

        # arrange functions in order of systems ((f0,f1,g0), (f0,f1,g1),
        # (f0,g0,g1), (f1,g0,g1)) and apply to each solution (found with sympy
        # linsolve)
        coeffs = [m(s) for m,s in zip((self.manipulatef0f1g0,
            self.manipulatef0f1g1, self.manipulatef0g0g1,
            self.manipulatef1g0g1), [tuple(sp.linsolve((m[:, :-1], m[:,-1]),
                sp.symbols('d1 d2 d3')))[0] for m in M])]

        solutions = [-0.5*cs[1]/cs[0] for cs in coeffs]
        
        return (solutions, self.x0, self.x1, self.f0, self.f1, self.g0,
                self.g1, coeffs, ("(f0,f1,g0)", "(f0,f1,g1)", "(f0,g0,g1)",
                    "(f1,g0,g1)"))
    # and function quadraticSolver
# end class Quad

if __name__ == "__main__":
    """ print and plot solution """
    quad = Quad()
    sols, x0, x1, f0, f1, g0, g1, coeffSols, order = quad.quadraticSolver()
    map(sp.pprint, sols)

    # make a cubic polynomial and substitute symbols with numerical values
    colors = ("r", "g", "y", "c")
    for i,sol in enumerate(sols):
        a = coeffSols[i][0]
        b = coeffSols[i][1]
        c = coeffSols[i][2]
        c3, c2, c1, c0 = 1.0, 1.0, 1./2, 1./3
        E = lambda x: c0*x**3 + c1*x**2 + c2*x + c3
        dE = lambda x: 3*c0*x**2 + 2*c1*x + c2
        x0Val = -1.0
        x1Val = 5.0
        f0Val = E(x0Val)
        f1Val = E(x1Val)
        g0Val = dE(x0Val)
        g1Val = dE(x1Val)
        valDict = {x0:x0Val, x1:x1Val, f0:f0Val, f1:f1Val, g0:g0Val, g1:g1Val}
        f = lambda C: a.subs(valDict) * C**2 + b.subs(valDict) * C + \
                c.subs(valDict)
        sol = sol.subs(valDict)

        import numpy as np
        import matplotlib.pyplot as plt
        x = np.linspace(-10,10,100)
        plt.plot(sol, E(sol), 'o%s' % colors[i], label="Extremal of quadratic "
                "fit, $%.3f$ %s" % (sol, order[i]))
        if i==0:
            """ only plot points once """
            plt.plot(x0Val, f0Val, 'pink', label="Interpolation point"
                    "$(x0,f0)=(%.2f,%.2f)$" % (x0Val, f0Val))
            plt.plot(x1Val, f1Val, 'teal', label="Interpolation point"
                    "$(x1,f1)=(%.2f,%.2f)$" % (x1Val, f1Val))
        # end if
        plt.plot(x, f(x), '%s-' % colors[i], label="Quadratic fit %s" % order[i])
    # end for i,sol
    plt.plot(x, E(x), 'b-', label="Cubic polynomial")
    plt.legend(loc='best')
    plt.show()
# end ifmain
