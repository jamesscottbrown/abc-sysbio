from math import *
import re

def mathMLConditionParser(mathMLstring):
    """
    Replaces and and or with and_ and or_ in a MathML string.
    Returns the string with and and or replaced by and_ and or_

    ***** args *****

    mathMLstring:

            A mathMLstring

    """
    
    andString = re.compile("and")
    orString = re.compile("or")
    mathMLstring = andString.sub("and_", mathMLstring)
    mathMLstring = orString.sub("or_", mathMLstring)
    return mathMLstring

def eq(a, b):
    """Check for equality between a and b.
    Return boolean.

    **** args ****

    a, b:

        python objects on which a test for equality will work

    """
    if a == b:
        return True
    else:
        return False

def neq(a, b):
    """Check for inequality between a and b.
    Return boolean

    **** args ****

    a, b:

        python objects on which a test for inequality will work

    """
    if a != b:
        return True
    else:
        return False

def gt(a, b):
    """Check for inequality a > b
    Return boolean

    **** args ****

    a, b:

        python objects on which a test for inequality will work

    """
    if a > b:
        return True
    else:
        return False

def lt(a, b):
    """Check for inequality a < b
    Return boolean

    **** args ****

    a, b:

        python objects on which a test for inequality will work

    """
    if a < b:
        return True
    else:
        return False

def geq(a, b):
    """Check for inequality a>=b
    Return boolean

    **** args ****

    a, b:

        python objects on which a test for inequality will work

    """

    if a >= b:
        return True
    else:
        return False

def leq(a, b):
    """Check for inequality a<=b
    Return boolean

    **** args ****

    a, b:

        python objects on which a test for inequality will work

    """
    if a <= b:
        return True
    else:
        return False

def and_(a, b):
    """Logical AND a and b
    Return boolean

    **** args ****

    a, b:

        python objects on which a python logical test will work.

    """
    if a & b:
        return True
    else:
        return False

def or_(a, b):
    """Logical OR a or b
    Return boolean

    **** args ****

    a, b:

        python objects on which a python logical test will work.

    """
    if a or b:
        return True
    else:
        return False

#NB piecewise functions have to be passed through MathMLConditionParser because
#they may contain and and or statements (which need to be made into and_ and or_
#statements)

def piecewise(*args):
    """implements MathML piecewise function passed as args
    NB piecewise functions have to be passed through MathMLConditionParser because
    they may contain and and or statements (which need to be made into and_ and or_
    statements)
    """
    
    result = None

    for i in range(1, len(args), 2):
        if args[i]:
            result = args[i-1]
            break
            
        else:
            result = args[-1]

    return result

