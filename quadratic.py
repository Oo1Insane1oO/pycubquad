# given two points x0 and x1 with corresponding function values f0 and f1 along
# with derivatives g0 and g1 at these points, find point which minimizes the
# quadratic interpolation through f0, f1 and g0

import sympy as sp

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

a = sp.collect(a, (g0,x0))
b = sp.collect(b, (g0,x0))

# split solution in numerator and denominator and restructure nominator
snom, sdenom = sp.fraction(-b/(2*a))
snomArgs = tuple(snom.args)
snom = sp.factor(snomArgs[0]) + snomArgs[1]

# rearrange numerator and denominator and substitute difference and sums of
# common symbols
sol = sp.simplify(sp.simplify(snom / sdenom).subs({x0-x1:sp.symbols("delta_x"),
    f0-f1:sp.symbols("delta_f"), x0+x1:sp.symbols("sigma_x")}))

# print results in a pretty fashion
sp.pprint(sol)

# test solution
E = lambda a,b,c,d,x: a*x**3 + b*x**2 + c*x + d
dE = lambda a,b,c,x: 3*a*x**2 + 2*b*x + c
c0, c1, c2, c3 = 1.0, 1.0, 0.5, 1./3
x0Val = -1.0
x1Val = 1.0
f0Val = E(c3, c2, c1, c0, x0Val)
f1Val = E(c3, c2, c1, c0, x1Val)
g0Val = dE(c3, c2, c1, x0Val)
valDict = {x0:x0Val, x1:x1Val, f0:f0Val, f1:f1Val, g0:g0Val}
f = lambda C: (a/denominators[0]).subs(valDict) * C**2 + \
        (b/denominators[1]).subs(valDict) * C + \
        (c/denominators[2]).subs(valDict)

# print f(sp.symbols("C"))

import numpy as np
import matplotlib.pyplot as plt
x = np.linspace(-10,10,100)
plt.plot(x0Val, f0Val, 'ro')
plt.plot(x1Val, f1Val, 'go')
plt.plot(x, f(x), 'b-')
plt.plot(x, E(c3, c2, c1, c0, x), 'y-')
plt.show()
