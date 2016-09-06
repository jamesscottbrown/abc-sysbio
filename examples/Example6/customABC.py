from numpy import *


def get_es(ts):
    t_o1 = [95, 96, 97, 98, 99]
    t_o2 = [196, 197, 198, 199, 200]
    t_sig = 100

    o1 = mean(ts[t_o1])
    o2 = mean(ts[t_o2])
    o_max = max(ts[t_sig:])
    o_min = min(ts[t_sig:])

    skip = False

    # check that we reached steady state before and after
    if sqrt(var(ts[t_o1])) > 0.01 * abs(o1):
        skip = True
    if sqrt(var(ts[t_o2])) > 0.01 * abs(o2):
        skip = True

    if abs(o_max - o1) > abs(o_min - o1):
        op = o_max
        imax = argmax(ts[t_sig:])
        o_next_min = min(ts[imax:])
        if abs(o_next_min - o1) > 0.5 * abs(op - o1):
            skip = True

    else:
        op = o_min
        imin = argmin(ts[t_sig:])
        o_next_max = max(ts[imin:])
        if abs(o_next_max - o1) > 0.5 * abs(op - o1):
            skip = True

    if not skip:
        e = abs(o2 - o1) / (0.2 * o1)
        s = abs(op - o1) / (0.2 * o1)
    else:
        e = None
        s = None

    return [e, s]


def distance(data1, data2, parameters, model):
    # data1 is simulated, and has shape npoints x beta
    # data2 is real

    d1, s = get_es(data1[:, 0])

    d2 = None
    if s > 0:
        d2 = 1 / s

    return [d1, d2]
