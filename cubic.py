# given two points x0 and x1 with corresponding function values f0 and f1 along
# with derivatives g0 and g1 at these points, find point which minimizes the
# cubic interpolation through f0, f1, g0 and g1

import sympy as sp

# set symbolic variables for f(x) = ax^3+bx^2+cx+d, with g=3ax^2+2bx+c being
# the derivative of f with respect to x, and interpolation points x0 and x1
a, b, c, d = sp.symbols('a b c d', real=True)
x0, x1 = sp.symbols('x_0 x_1', real=True, positive=True)
f0, f1, g0, g1 = sp.symbols('f_0 f_1 g_0 g_1', real=True)

# define linear system as a matrix with the columns corresponding to a, b, c
# and d respectively
A = sp.Matrix((
    [x0**3, x0**2, x0, 1, f0], 
    [x1**3, x1**2, x1, 1, f1], 
    [3*x0**2, 2*x0, 1, 0, g0], 
    [3*x1**2, 2*x1, 1, 0, g1]))

# split system in lhs matrix and rhs vector (Ay=r, with A the matrix, y the
# vector y=(a,b,c,d) and r the rhs vector r=(f0,f1,g0,g1)) and solve
solution = tuple(sp.linsolve((A[:, :-1], A[:,-1]), a,b,c,d))[0]

# split nominator and denominator in each coefficient (denominator is the same
# for all coefficients)
nominators = []
denominators = []
for s in solution:
    nom, denom = sp.fraction(sp.factor(sp.powsimp(s)))
    nominators.append(sp.simplify(nom))
    denominators.append(denom)

# update symbolic variables (ignore denominator)
a, b, c, d = nominators

# manipulate and simplify coefficients (still ignoring denominator)
aArgs = tuple(nominators[0].args)
a = sp.collect(sp.factor(aArgs[0] + aArgs[1]) + sp.simplify(aArgs[2]+aArgs[3])
        + sp.simplify(aArgs[4]+aArgs[5]), g0+g1)

bArgs = tuple(nominators[1].args)
b = sp.simplify(sp.factor(bArgs[2] + bArgs[3] + bArgs[6] + bArgs[7]) +
        sp.factor(bArgs[0] + bArgs[1] + bArgs[4] + bArgs[5] + bArgs[8] +
            bArgs[9]))

cArgs = tuple(nominators[2].args)
c = sp.simplify(sp.factor(cArgs[4] + cArgs[7]) + sp.collect(sp.factor(cArgs[0]
    + cArgs[1] + cArgs[2] + cArgs[3] + cArgs[5] + cArgs[6]), x0*x1))

dArgs = tuple(nominators[3].args)
d = sp.collect(sp.simplify(sp.factor(dArgs[0] + dArgs[1]) + sp.factor(dArgs[6]
    + dArgs[7]) + sp.factor(dArgs[2] + dArgs[3] + dArgs[4] + dArgs[5])), x0*x1)

# substitute difference and sums of common symbols
subser = lambda expr: expr.subs({x0-x1:sp.symbols("delta_x", real=True),
    f0-f1:sp.symbols("delta_f", real=True), g0-g1:sp.symbols("delta_g",
        real=True), x0+x1:sp.symbols("sigma_x", real=True),
    f0+f1:sp.symbols("sigma_f", real=True), g0+g1:sp.symbols("sigma_g",
        real=True)})
# a = subser(a)
# b = subser(b)
# c = subser(c)
# d = subser(d)

# solve quadratic equation to find interpolation points x_t+ and x_t- which
# minimize the function (set derivative g(x_t) to zero and solve for x_t)
sqrtFactor = sp.sqrt(1 - sp.factor(3*a*c/b**2))
xtp = (-1 + sqrtFactor) / (3*a)
xtm = (-1 - sqrtFactor) / (3*a)

# print results in a pretty fashion
sp.pprint(xtp)
sp.pprint(xtm)

# test solution
E = lambda a,b,c,d,e,x: a*x**4 + b*x**3 + c*x + d*x + e
dE = lambda a,b,c,d,x: 4*a*x**3 + 3*b*x**2 + 2*c*x + d
c0, c1, c2, c3, c4 = 1.0, 1.0, 0.5, 1./3, 0.25
x0Val = -10.5
x1Val = 2.0
f0Val = E(c4, c3, c2, c1, c0, x0Val)
f1Val = E(c4, c3, c2, c1, c0, x1Val)
g0Val = dE(c4, c3, c2, c1, x0Val)
g1Val = dE(c4, c3, c2, c1, x1Val)
valDict = {x0:x0Val, x1:x1Val, f0:f0Val, f1:f1Val, g0:g0Val, g1:g1Val}
f = lambda C: (a/denominators[0]).subs(valDict) * C**3 + \
        (b/denominators[1]).subs(valDict) * C**2 + \
        (c/denominators[2]).subs(valDict) * C + \
        (d/denominators[3]).subs(valDict)

import numpy as np
import matplotlib.pyplot as plt
x = np.linspace(-20,10,100)
plt.plot(x0Val, f0Val, 'ro')
plt.plot(x1Val, f1Val, 'go')
plt.plot(x, f(x), 'b-')
plt.plot(x, E(c4, c3, c2, c1, c0, x), 'y-')
plt.show()
