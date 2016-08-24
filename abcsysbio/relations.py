from math import *
import re


def mathMLConditionParser(string):
    """
    Replaces and and or with and_ and or_ in a string and returns the result.
    """

    string = re.compile("and").sub("and_", string)
    string = re.compile("or").sub("or_", string)
    return string


def eq(a, b):
    return a == b


def neq(a, b):
    return a != b


def gt(a, b):
    return a > b


def lt(a, b):
    return a < b


def geq(a, b):
    return a >= b


def leq(a, b):
    return a <= b


def and_(a, b):
    return bool(a) & bool(b)


def or_(a, b):
    return bool(a) or bool(b)


def piecewise(*args):
    """
    Implements MathML piecewise function passed as args.

    N.B. piecewise functions have to be passed through MathMLConditionParser because they may contain and and or
    statements (which need to be made into and_ and or_ statements).
    """

    result = None

    for i in range(1, len(args), 2):
        if args[i]:
            result = args[i - 1]
            break

        else:
            result = args[-1]

    return result
