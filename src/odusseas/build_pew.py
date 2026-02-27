import numpy as np


def cut_data(w, f, w1=None, w2=None):
    if w1 is None:
        w1 = w[0]
    if w2 is None:
        w2 = w[-1]
    idx = (w1 <= w) & (w <= w2)
    return w[idx], f[idx]


def area_between(f, g, dx):
    h = abs(g - f) / g
    return np.trapezoid(h, dx=dx)
