import numpy as np
import sys
from rmt_utils import olrp, ricatti
from matrix2latex import matrix2latex as to_tex

bmat_tex = lambda i: to_tex(i, None, 'bmatrix', formatColumn=['%.3f'] * 15)

args = sys.argv
probs = args[1:] if len(args) > 0 else []


def p5_3(beta, r, b, gamma, rho_1, rho_2):
    """
    Using rmt_utils.olrp to solve the optimal linear regulator problem.
    from RMT4 problem 5.3b
    """
    A = np.array([[1., 0., 0., 0.],
                  [0., rho_1, rho_2, 0.],
                  [0., 1., 0., 0.],
                  [0., 0., 0., 1.]])

    B = np.array([1., 0., 0., 0.]).reshape(4, 1)

    R = np.array([[r**2, r, 0., -b*r],
                 [r, 1., 0, -b],
                 [0., 0., 0., 0.],
                 [-b*r, -b, 0, b**2]])

    Q = 1 + gamma

    H = np.array([[-r, -1., 0., b]]).reshape(1, 4)

    return ricatti(beta, A, B, R, Q, H)

if '53' in probs:
    in_53 = [0.95, 1 / 0.95 - 1., 30, 1., 1.2, -0.3]
    f, p = p5_3(*in_53)


def p5_4(beta, r, b, gamma, rho_1, rho_2, lamb, pi, delta, theta):
    """
    Using rmt_utils.olrp to solve the optimal linear regulator problem.
    from RMT4 problem 5.4(b,c)
    """
    A = np.array([[1, 0, 0, 0, 0],
                 [0, rho_1, rho_2, 0, 0],
                 [0,  1,  0, 0,  0],
                 [theta * r, theta, 0, delta, 0],
                 [0,  0,  0,  0,  1]]
                 )

    B = np.array([1, 0, 0, -theta, 0.])

    Q = pi**2 + gamma
    R = np.array([[(pi*r)**2, pi**2*r, 0, lamb*pi*r, -b*pi*r],
                  [pi**2*r, pi**2, 0, lamb*pi, -b*pi],
                  [0, 0, 0, 0, 0],
                  [lamb*pi*r, lamb*pi, 0, lamb**2, -b*lamb],
                  [-b*pi*r, -b*pi, 0, -b*lamb, b**2]])

    H = np.array([-pi**2*r, -pi**2, 0, -pi*lamb, b*pi])

    return ricatti(beta, A, B, R, Q, H)

if '54' in probs:
    in_54b = [0.95, 1 / 0.95 - 1, 30, 1, 1.2, -0.3, 1, 0.5, 0.95, 1.0]
    in_54c = [0.95, 1 / 0.95 - 1, 30, 1, 1.2, -0.3, -1, 1, 0.95, 1.0]
    f4b, p4b = p5_4(*in_54b)
    f4c, p4c = p5_4(*in_54c)
