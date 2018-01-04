# write expression for minimum of interpolation from quadratic or cubic to C++
# code as a C++ template function

from quadratic import quadraticSolver
from cubic import cubicSolver

import sympy as sp
import sys

def turnToCpp(funcName, expression, sumDiffSubDict, *args):
    """ turn expression into C++ code with sympy and save template function
    expression to string """
    subsDict = {};
    argString = ""
    for arg in args:
        argName = str(arg)
        argString += "T " + argName + ", "
        subsDict[argName] = arg
    # end forarg
    argString = argString[:-2]
    subsDict['expr'] = \
            sp.printing.ccode((expression.subs(subsDict)).subs(sumDiffSubDict))
    return '''template<typename T> T ''' + funcName + ('''(''' + argString +
            ''') {\n''' + ''''''.join(["     T "+str(val)+" = "+str(key)+";\n" for
                key,val in sumDiffSubDict.iteritems()]) + '''    return '''
            '''%(expr)s;\n}''') % (subsDict)
# and function turnToCpp

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
        returned = quadraticSolver()
        sol = returned[0]
        args = returned[1:-3]

        # substitute difference and sums of common symbols
        x0, x1, f0, f1, g0 = args
        subser = {x0-x1:sp.symbols("diff_x", real=True),
                f0-f1:sp.symbols("diff_f", real=True),
                x0+x1:sp.symbols("sum_x", real=True)}
        funcName = "quadPol"
    elif mode == "cubic":
        """ find quadratic interpolation """
        print("finding cubic fit")
        returned = cubicSolver()
        sol = returned[:2]
        args = returned[2:-4]
        
        # substitute difference and sums of common symbols
        x0, x1, f0, f1, g0, g1 = args
        subser = {x0-x1:sp.symbols("diff_x", real=True),
                f0-f1:sp.symbols("diff_f", real=True),
                g0-g1:sp.symbols("diff_g", real=True),
                x0+x1:sp.symbols("sum_x", real=True), g0+g1:sp.symbols("sum_g",
                    real=True)}
        funcName = "cubicPol"
    else:
        """ print error message """
        print("Specify interpolation method! (quad/cubic)")
        sys.exit(1);
    # end ifeifelse
       
    # turn expression into string representing C++ template and write to file
    print ("turning into C++ code...")
    if mode == "cubic":
        """ two solutions to write with cubic interpolation """
        codes = [turnToCpp(funcName+"Plus", sol[0], subser, *args),
                turnToCpp(funcName+"Minus", sol[1], subser, *args)]
    else:
        """ only one solution to write with quadratic interpolation """
        codes = [turnToCpp(funcName, sol, subser, *args)]
    # end ifelse

    # write code to file filename
    with open(filename, "w") as writeFile:
        """ open file for writing """
        for c in codes:
            writeFile.write(c+"\n")
        # end forc
    # end with open
# and ifmain
