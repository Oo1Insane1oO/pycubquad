# write expression for minimum of interpolation from quadratic or cubic to C++
# code as a C++ template function

from quadratic import *
from cubic import cubicSolver

import sympy as sp
import sys

def makeArgsDict(args):
    """ make a string for function arguments and a dictionary with variable as
    key and expression as value. Return both string and dictionary """
    subsDict = {}
    argString = ""
    for arg in args:
        argName = str(arg)
        argString += "T " + argName + ", "
        subsDict[argName] = arg
    # end forarg
    argString = argString[:-2]

    return subsDict, argString
# end function makeArgsDict

def turnToCppQuad(funcName, expression, sumDiffSubDict, *args):
    """ turn expression into C++ code with sympy and save template function
    expression to string """
    subsDict, argString = makeArgsDict(args)
    subsDict['expr'] = \
            sp.printing.ccode((expression.subs(subsDict)).subs(sumDiffSubDict))

    return '''template<typename T> T ''' + funcName + ('''(''' + argString +
            ''') {\n''' + ''''''.join(["    T "+str(val)+" = "+str(key)+";\n" for
                key,val in sumDiffSubDict.iteritems()]) + '''    return '''
            '''%(expr)s;\n}''') % (subsDict)
# and function turnToCppQuad

def turnToCppCubic(funcName, expressionPlus, expressionMinus, c0, c1,
        sumDiffSubDict, *args):
    """ turn expression into C++ code with sympy and save template function
    expression to string """
    subsDict, argString = makeArgsDict(args)
    subsDict['exprp'] = sp.printing.ccode((expressionPlus.subs(subsDict)) . \
            subs(sumDiffSubDict))
    subsDict['exprm'] = sp.printing.ccode((expressionMinus.subs(subsDict)) . \
            subs(sumDiffSubDict))
    subsDict['secExpr'] = sp.printing.ccode(((6*c0*sp.symbols("xtp",real=True)
        + 2*c1).subs(subsDict)).subs(sumDiffSubDict))

    return '''template<typename T> T ''' + funcName + ('''(''' + argString + ''') {\n''' + ''''''.join(["     T "+str(val)+" = "+str(key)+";\n" for key,val in sumDiffSubDict.iteritems()]) + '''   T xtp = %(exprp)s;\n    if (%(secExpr)s > 0) {\n        return xtp;\n    } else {\n        return %(exprm)s;\n    }\n}''') % (subsDict)
# end function turnToCppCubic

if __name__ == "__main__":
    try:
        filename = sys.argv[1]
        mode = sys.argv[2].lower()
    except IndexError:
        print ("IndexError: USAGE: python writeToCpp.py filename "
                "quadratic/cubic")
        sys.exit(0)
    except ValueError:
        print ("ValueError: USAGE: python writeToCpp.py filename "
                "quadratic/cubic")
        sys.exit(0)
    # end try-except

    # find expression
    if mode == "quadratic" or mode == "quad":
        """ find quadratic interpolation """
        print("finding quadratic fit...")
        quad = Quad()
        returned = quad.quadraticSolver()
        sols = returned[0]
        args = returned[1:7]
        order = [i.replace("(","").replace(")","").replace(",","") for i in
                returned[-1]]

        # substitute difference and sums of common symbols
        x0, x1, f0, f1, g0, g1 = args
        xdiff, fdiff, gdiff = (sp.symbols("diff_x", real=True),
                sp.symbols("diff_f", real=True), sp.symbols("diff_g",
                    real=True))
        funcName = "quadPol"
        print ("turning into C++ code...")
        codes = [
                turnToCppQuad(funcName+order[0], sols[0], {x0-x1:xdiff,
                    f0-f1:fdiff}, x0, x1, f0, f1, g0),
                turnToCppQuad(funcName+order[1], sols[1], {x0-x1:xdiff,
                    f0-f1:fdiff}, x0, x1, f0, f1, g1),
                turnToCppQuad(funcName+order[2], sols[2], {}, x0, x1, f0, g0,
                    g1),
                turnToCppQuad(funcName+order[3], sols[3], {}, x0, x1, f1, g0,
                    g1)
                ]
    elif mode == "cubic":
        """ find quadratic interpolation """
        print("finding cubic fit")
        returned = cubicSolver()
        sol = returned[:2]
        args = returned[2:-4]
        a, b = returned[-4:-2]
        
        # substitute difference and sums of common symbols and second
        x0, x1, f0, f1, g0, g1 = args
        subser = {x0-x1:sp.symbols("diff_x", real=True),
                f0-f1:sp.symbols("diff_f", real=True),
                g0-g1:sp.symbols("diff_g", real=True),
                x0+x1:sp.symbols("sum_x", real=True), g0+g1:sp.symbols("sum_g",
                    real=True)}
        # second derivative for minimum
        funcName = "cubicPol"
        print ("turning into C++ code...")
        codes = [turnToCppCubic(funcName, sol[0], sol[1], a, b, subser, *args)]
    else:
        """ print error message """
        print("Specify interpolation method! (quad/cubic)")
        sys.exit(1);
    # end ifeifelse

    # write code to file filename
    with open(filename, "w") as writeFile:
        """ open file for writing """
        for c in codes:
            writeFile.write(c+"\n")
        # end forc
    # end with open
# and ifmain
