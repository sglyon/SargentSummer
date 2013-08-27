from __future__ import division
import sympy as sym

k = sym.symbols('k', integer=True, positive=True)
z = sym.symbols('z', complex=True)


def sym_ztrans(f, k=k, z=z):
    """
    f needs to be a function of k
    """

    return sym.summation(f * z ** (-k), (k, 0, sym.oo))
