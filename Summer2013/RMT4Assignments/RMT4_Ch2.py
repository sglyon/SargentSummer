from __future__ import division
import sys
import numpy as np
from scipy.linalg import inv, eig
import pandas as pd
from rmt_utils import doublej
from matrix2latex import matrix2latex as to_tex

import sympy as sym
if len(sys.argv) == 2:
    verb = sys.argv[1]
else:
    verb = False

bmat_tex = lambda i: to_tex(i, None, 'bmatrix', formatColumn=['%.3f'] * 15)


class StochasticLinearDiff(object):
    """
    Represents and computes various things for a model in the form
    of the canonical stochastic linear difference equation:

    .. math::

        x_{t+1} = A x_t + C w_{t+1}
    """

    def __init__(self, A, C):
        self.A = A
        self.C = C

        # Evaluate eigenvalues and vectors for use later on. Check boundedness
        evals, evecs = eig(self.A, left=False, right=True)
        self.evals, self.evecs = evals, evecs
        self.unbounded = np.abs(evals).max() > 1

    def Cx(self, j=0):
        "Covariance stationary covariance matrix"
        if not self.unbounded:
            c_x = doublej(self.A, self.C.dot(self.C.T))

            # Return if we want C_x(0)
            if j == 0:
                return c_x
            else:
                # Or evaluate C_x(abs(j))
                c_xj = np.linalg.matrix_power(self.A, abs(j)).dot(c_x)
            if j < 0:
                return c_xj.T  # transpose if j < 0
            else:
                return c_xj

        else:
            msg = 'This computation will not work because the eigenvalues'
            msg += '\nof A are not all below 1 in modulus.'
            raise ValueError(msg)

    @property
    def mu(self):
        "Covariance stationary mean"
        if self.unbounded:
            msg = 'This computation will not work because the eigenvalues {0}'
            msg += '\nof A are not all below 1 in modulus.'
            raise ValueError(msg.format(self.evals))

        # Try to get index of unit eigenvalue
        try:
            ind = np.where(self.evals == 1)[0][0]
        except IndexError:
            raise ValueError("The A matrix doesn't have any unit eigenvalues")

        # compute Stationary mean using the eigenvector for unit eigenvalue
        return self.evecs[:, ind] / self.evecs[-1, ind]


class Markov(object):
    """
    Do basic things with Markov matrices.
    """

    def __init__(self, P, verbose=False):
        self.P = P
        self.verbose = verbose

    def __repr__(self):
        msg = "Markov process with transition matrix P = \n{0}"
        return msg.format(self.P)

    def stationary_distributions(self):
        evals, l_evecs, r_evecs = eig(self.P, left=True, right=True)
        self.evals, self.l_evecs, self.r_evecs = evals, l_evecs, r_evecs
        units = np.where(evals == 1)[0]
        stationary = []
        for i, ind in enumerate(units):
            sd_name = 'sd{0}'.format(i + 1)
            sd_vec = l_evecs[:, ind]

            # Normalize to be probability vector
            sd_vec = sd_vec * (-1) if all(sd_vec <= 0) else sd_vec
            sd_vec /= sd_vec.sum()
            self.__setattr__(sd_name, sd_vec)
            stationary.append(sd_vec)
            if self.verbose:
                msg = 'Set instance variable %s for stationary distribution'
                print(msg % sd_name)
        return stationary

    def invariant_distributions(self):
        units = np.where(self.evals == 1)[0]
        invariant = []
        for i, ind in enumerate(units):
            id_name = 'id{0}'.format(i + 1)
            id_vec = self.r_evecs[:, ind]

            # Normalize to be probability vector
            id_vec = id_vec * (-1) if all(id_vec <= 0) else id_vec
            id_vec /= id_vec.sum()
            self.__setattr__(id_name, id_vec)
            invariant.append(id_vec)
            if self.verbose:
                msg = 'Set instance variable %s for invariant distribution'
                print(msg % id_name)
        return invariant


class SymMarkov(object):
    """
    Do basic things with Markov matrices. The matrix P that is passed
    to the constructor for this class is assumed to be a sympy matrix.
    If it isn't, then it is cast as such.
    """

    def __init__(self, P, verbose=False):
        self.P = P if isinstance(P, sym.Matrix) else sym.Matrix(P)
        self.verbose = verbose

    def stationary_distributions(self, subs, normalize=True):
        """
        Find the stationary distributions associated with the Markov
        process, by substituting parameters into the transition matrix

        Parameters
        ==========
        subs : dist
            A dictionary of substitutions to be passed to self.P before
            doing the computation

        normalize : bool, optional(default=True)
            Whether or not the stationary distributions should be
            normalized so they sum to 1 before returning.

        Returns
        =======
        pi0s : list
            A list of stationary distributions.

        """
        # Make the substitutions
        PN = self.P.subs(subs)

        # Transpose gives left eigenvectors
        l_vecs = PN.T.eigenvects()

        # keep only unit eigenvalues around, grab the vectors
        units = filter(lambda x: x[0] == 1, l_vecs)
        pi0s = units[0][2] if len(units) != 0 else []

        # Normalize so they sum to 1
        if normalize:
            pi0s = [i / sum(i) for i in pi0s]

        return pi0s


def isoelastic_util(c, gamma):
    """
    Implementation of the isoelastic utility function
    :math:`u = \\frac{c^{1 - \gamma}}{1 - \gamma}`
    """
    return c ** (1 - gamma) / (1 - gamma)


def p2_3b(P, beta=0.95, pi_0=[.5, .5], cbar=[1, 5], gamma=[2.5, 4.],
          u=isoelastic_util):
    """
    Solution to problem 2.3b for RMT4

    Parameters
    ==========
    P : iterable of array_like, dtype=float
        A list of transition matrices to be used as different processes
        consumption could follow.

    beta : float, optional(default=0.95)
        The value of the discount factor :math:`\\beta`

    pi_0: array_like, dtype=float, optional(default=[.5, 5])
        The value of the initial distribution :math:`\\pi_0`

    cbar : array_like, dtype=float, optional(default=[.5, 5])
        An array representing the consumption state space
        :math:`\\bar{c}`

    gamma : array_like, dtype=float, optional(default=[2.5, 4])
        Different values for the coeeficient of relative risk aversion
        :math:`gamma` in the utility function.

    u : function, optional(default=isoelastic_util)
        A python function representing the utility function. It is
        assumed that this function takes two arguments: 1. the value of
        consumption and 2. the value of gamma

    Returns
    =======
    V : pandas.DataFrame
        A pandas DataFrame representing the discounted expected utility
        :math:`V` given the parameters.
    """
    ## Define parameters
    beta = 0.95
    pi_0 = np.asarray(pi_0)
    cbar = np.asarray(cbar)

    # Prepare DataFrame to hold results
    c_names = ['Process%i' % (i + 1) for i in range(len(P))]
    V = pd.DataFrame(index=gamma, columns=c_names)

    for gam in gamma:
        u = cbar ** (1 - gam) / (1 - gam)
        for i, p in enumerate(P):
            v = inv(np.eye(2) - beta * p).dot(u)
            V.ix[gam, i] = v.dot(pi_0)

    return V

V = p2_3b(P=[np.eye(2), np.ones((2, 2)) * 0.5])

if verb:
    print(V.to_latex())

# Problem 2.4
case4_1 = ([1.2, -0.3, 0, 0], 10, 1)
case4_2 = ([1.2, -0.3, 0, 0], 10, 2)
case4_3 = ([0.9, 0, 0, 0], 5, 1)
case4_4 = ([0.2, 0, 0, 0.5], 5, 1)
case4_5 = ([0.8, 0.3, 0, 0], 5, 1)

cases4 = [case4_1, case4_2, case4_3, case4_4, case4_5]


def _beta2tex(beta):
    bet = pd.DataFrame(beta)
    return bet.T.applymap(lambda x: '%.3f' % x).to_latex(index=False,
                                                         header=False)


def p2_4(rho, mu, c, print_results=False):
    """
    Solution to problem 2.4 parts b-e for RMT4

    Problem Text
    ============

    Consider the univariate stochastic process

    .. math::

        y_{t+1} = \\alpha + \sum_{j=1}^4 \\rho_j y_{t+1 - j} + c w_{t+1}

    where :math:`w_{t+1}` is a scalar martingale difference sequence
    adapted to :math:`J_t = \left[ w_t, \dots w_1, y_0, y_{-1}, y_{-2},
    y_{-3} \\right], \\alpha = \mu \left( a - \sum_j \\rho_j \\right)` and
    the :math:`\\rho_j` 's are such that the matrix

    ..math ::

        \\begin{bmatrix}
        \\rho_1 & \\rho_2 & \\rho_3 & \\rho_4 & \\alpha \\
        1 &0 & 0 & 0 & 0 \\
        0 & 1 & 0 & 0 & 0 \\
        0 & 0 & 1 & 0 & 0 \\
        0 & 0 & 0 & 0 & 1
        \end{bmatrix}

    has all of its eigenvalues in modulus bounded below unity.

    # TODO: This is incomplete...

    Parameters
    ==========
    rho : array_like, dtype=float
        The array-like object that corresponds to
        :math:`[\\rho_1 \\rho_2 \\rho_3 \\rho_4]` in the problem

    mu : scalar, dtype=float
        The sclar :math:`\mu` from the problem

    c : scalar, dtype=float
        The value of the :math:`c` in the problem.

    print_results : bool
        Whether or not the function should print its results.

    Returns
    =======
    mu : array, dtype=float
        The array of average values. Represents the vector :math:`\mu`

    Cx : array, dtype=float
        The stationary variance/co-variance matrix for the problem
    """
    # Solve for alpha
    alpha = mu - mu * sum(rho)

    # Create a1 (A) matrix
    r1, r2, r3, r4 = rho

    a1 = np.array([[r1, r2, r3, r4, alpha],
                   [1, 0, 0, 0, 0],
                   [0, 1, 0, 0, 0],
                   [0, 0, 1, 0, 0],
                   [0, 0, 0, 0, 1],
                   ], dtype=float)

    # Check to make sure all eigenvalues are bounded below q
    w, vr = eig(a1, left=False, right=True)
    unbounded = np.abs(w).max() > 1
    if unbounded:
        msg = 'This computation will not work because the eigenvalues {0}\n'
        msg += 'of A are not all below 1 in modulus.'
        raise ValueError(msg.format(w))

    # Try to get index of unit eigenvalue
    try:
        ind = np.where(w == 1)[0][0]
    except IndexError:
        raise ValueError('The A matrix does not have any unit eigenvalues')

    # compute Stationary mean using the eigenvector for unit eigenvalue
    mu = vr[:, ind] / vr[-1, ind]

    # Solve for the stationary covariance Cx
    # Create b1 (CC') matrix
    _c = np.array([c, 0, 0, 0, 0])
    b1 = np.outer(_c, _c)

    Cx = doublej(a1, b1)

    # Part c
    mumu = np.outer(mu, mu)
    G = np.array([1, 0, 0, 0, 0])
    term1 = inv(Cx + mumu)
    term2 = (np.dot(Cx.T, np.linalg.matrix_power(a1, 5).T) + mumu).dot(G)
    beta = np.dot(term1, term2)

    # Part d
    beta_d = inv(np.eye(5) - 0.95 * a1).T.dot(G.T)

    # Part e
    auto_covar = [np.linalg.matrix_power(a1, 1).dot(Cx)[0, 0],
                  np.linalg.matrix_power(a1, 5).dot(Cx)[0, 0],
                  np.linalg.matrix_power(a1, 10).dot(Cx)[0, 0]]

    if print_results:
        msg = 'Results for when rho = {0}, mu = {1} and c = {2}:'
        print(msg.format(rho, mu, c))
        evals = np.abs(w) if all(w == np.abs(w)) else w
        print('The eigenvalues are: %s' % evals)
        print('The stationary mean for this problem is:\n%s' % mu)
        print('The stationary covariance for this problem is:\n%s\n' % (Cx))
        print('The value of beta in part c is: \n%s\n' % (beta))
        print('The value of beta in part d is: \n%s\n' % (beta_d))
        print('The autocovariances are k=1: %.3f, k=5: %.3f, k=10: %.3f'
              % tuple(auto_covar))

    return mu, Cx, beta, beta_d, auto_covar

roman_numerals = ['i', 'ii', 'iii', 'iv', 'v']
auto_covars = pd.DataFrame(index=[1, 5, 10], columns=roman_numerals)
auto_covars.index.name = 'k'
auto_covars.columns.name = 'Parameterization'

for i, case in enumerate(cases4):
    try:
        # Get solution
        ux, Cx, beta_c, beta_d, a_cvar = p2_4(*case, print_results=verb)
        auto_covars[roman_numerals[i]] = a_cvar

        if verb:
            # Print latex form for copy/paste
            print('------------ CASE {0} -------------'.format(i+1))
            msg = '\n\nmu:\n{0}\n\nCx:\n{1}\n\nbeta c:\n{2}\n\nbeta d:\n{3}'
            print(msg.format(to_tex(np.abs(ux)),
                             to_tex(Cx),
                             _beta2tex(beta_c),
                             _beta2tex(beta_d)))
    except ValueError as e:  # don't let the eigenvalue error stop computation
        if verb:
            print('\n\n' + '*' * 72)
            err_msg = 'When rho = {0}, mu = {1} and c = {2} we had this error:'
            print(err_msg.format(*case))
            print(e)

if verb:
    print('\n\nThe autocovariances for part e are:\n ')
    print(auto_covars.to_latex())


# Problem 2.5

# Define parameterizations. give (rho, alpha, delta, gamma, phi, psi1, psi2)
case5_1 = ([0.8, -0.3], 1., [0.2, 0], [0., 0.], [0.7, -0.2], 1., 1.)
case5_2 = ([0.8, -0.3], 1., [0.2, 0], [0., 0.], [0.7, -0.2], 2., 1.)

cases5 = [case5_1, case5_2]


def p2_5(rho, alpha, delta, gamma, phi, psi1, psi2):
    """
    Solution to problem 2.5 part b for RMT4

    Problem Text
    ============

    TODO: Fill this in

    Parameters
    ==========
    rho : array_like, dtype=float
        The value of rho from the problem description

    alpha : float
        The value of alpha from the problem description

    delta : array_like, dtype=float
        The value of delta from the problem description

    gamma : array_like, dtype=float
        The value of gamma from the problem description

    phi : array_like, dtype=float
        The value of phi from the problem description

    psi1 : float
        The value of psi1 from the problem description

    psi2 : float
        The value of psi2 from the problem description

    Returns
    =======
    R : array-like, dtype=float
        The value of the matrix R that determines the solution for V_0

    xi : float
        The value of the scalar xi that determines the solution for V_0

    """
    A = np.array([
                 [rho[0], rho[1], delta[0], delta[1], alpha],
                 [1, 0, 0, 0, 0],
                 [gamma[0], gamma[1], phi[0], phi[1], 0],
                 [0, 0, 1, 0, 0],
                 [0, 0, 0, 0, 1]
                 ], dtype=float)

    C = np.array([[psi1, 0],
                  [0, 0],
                  [0, psi2],
                  [0, 0],
                  [0, 0]
                  ], dtype=float)

    G = np.array([
                 [-0.5, 0, 0, 0, 30],
                 [0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0],
                 [30, 0, 0, 0, -1800]
                 ], dtype=float)

    R = doublej(0.95 ** (1. / 2) * A, G)
    xi = (0.95) / (1 - 0.95) * np.dot(C.T.dot(R), C).trace()

    return R, xi


for i, case in enumerate(cases5):
    R_i, xi_i = p2_5(*case)
    if verb:
        msg = '\n\nR:\n{0:}\n\nxi:\n{1:.3e}'
        print('------------ CASE {0} -------------'.format(i+1))
        print(msg.format(
              bmat_tex(R_i, formatColumn=['$%.3e$'] * 5), xi_i))


def p2_6(rho, alpha, delta, gamma, phi, psi1, psi2, verbose=False):
    """
    Solution to problem 2.5 part b for RMT4

    Problem Text
    ============

    TODO: Fill this in

    Parameters
    ==========
    rho : array_like, dtype=float
        The value of rho from the problem description

    alpha : float
        The value of alpha from the problem description

    delta : array_like, dtype=float
        The value of delta from the problem description

    gamma : array_like, dtype=float
        The value of gamma from the problem description

    phi : array_like, dtype=float
        The value of phi from the problem description

    psi1 : float
        The value of psi1 from the problem description

    psi2 : float
        The value of psi2 from the problem description

    """
    A = np.array([
                 [rho[0], rho[1], delta[0], delta[1], alpha],
                 [1, 0, 0, 0, 0],
                 [gamma[0], gamma[1], phi[0], phi[1], 0],
                 [0, 0, 1, 0, 0],
                 [0, 0, 0, 0, 1]
                 ], dtype=float)

    C = np.array([[psi1, 0],
                  [0, 0],
                  [0, psi2],
                  [0, 0],
                  [0, 0]
                  ], dtype=float)

    # part A
    diff_eq = StochasticLinearDiff(A, C)

    cx = diff_eq.Cx()
    mu = diff_eq.mu

    # part b
    e_xx = np.array([[1., 0., 0.],
                     [0., cx[2, 2], diff_eq.Cx(4)[2, 2]],
                     [0., diff_eq.Cx(-4)[2, 2], cx[2, 2]]])

    e_xy = np.array([2, diff_eq.Cx(2)[0, 2], diff_eq.Cx(5)[0, 3]])

    beta = e_xy.dot(inv(e_xx))

    if verbose:
        print('\n\nPROBLEM 2.6 %s\n\n' % ('#' * 65))
        vals = "C_X(0) = {0}\nmu = {1}\nE(XX') = {2}\nE(YX') = {3}\nbeta= {4}"
        print(vals.format(bmat_tex(cx),
                          bmat_tex(np.abs(mu)),
                          bmat_tex(e_xx),
                          bmat_tex(e_xy),
                          bmat_tex(beta)))

    return cx, mu, beta

cx, mu, beta = p2_6(*case5_2, verbose=False)

P_14 = np.array([[.5, .5, .0, .0],
                 [.1, .9, .0, .0],
                 [.0, .0, .9, .1],
                 [.0, .0, .0, 1]])


def p2_14(P=P_14, verbose=False):
    """
    Solves for stationary distributions and long-run averages for
    :math:`y_t = \\bar{y} x_t = [1 2 3 4] x_t` given a transition matrix
    P

    Parameters
    ==========
    P : array_like, dtype=float
        The 2d numpy array representing the Markov transition matrix

    Returns
    =======
    markov : Markov
        The object of type Markov representing the process
    stationary : array_like, dtype=float
        The list of arrays representing stationary distributions of the
        process P
    invariant : array_like, dtype=float
        The list of arrays representing the invariant distributions of
        the process P
    lr_means : array_like, dtype=float
        The long-run mean of the process under each stationary
        distribution

    """
    markov = Markov(P)
    stationary = markov.stationary_distributions()
    invariant = markov.invariant_distributions()

    ybar = np.array([1, 2, 3, 4])

    # Long run mean values
    lr_means = [markov.sd1.dot(ybar), markov.sd2.dot(ybar)]

    if verbose:
        print('\n\nPROBLEM 2.14 %s\n\n' % ('#' * 65))
        vals = "\\pi_0 = {0}\n\\pi_1 = {1}"
        print(vals.format(bmat_tex(stationary[0]), bmat_tex(stationary[1])))

    return markov, stationary, invariant, lr_means


markov_14, _, _, means = p2_14(verbose=False)


def p2_17(ll=0.05, dd=0.25):
    lamb, delta = sym.symbols('lambda delta')
    P = sym.Matrix([[1 - lamb, lamb], [delta, 1 - delta]])
    markov = SymMarkov(P)

    # Normalize so they sum to 1
    pi0s = markov.stationary_distributions({lamb: ll, delta: dd})

    # If there is only one, return it, not a list
    pi0s = pi0s if len(pi0s) > 1 else pi0s[0] if len(pi0s) > 0 else None

    # part c

    return markov, pi0s


def p2_20(sigma_0=1.36602540378444, t=5000, plot=True, save_fig=True):
    sigma = np.zeros(t)
    sigma[0] = sigma_0
    for i in range(1, t):
        sigma[i] = (3 * sigma[i - 1] + 1) / (2 * sigma[i - 1] + 1)

    a_gk = 1 / (2 * sigma + 1)
    if plot:
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(range(t), a_gk, 'k--')
        ax.set_title(r'Stability of $A - K_t G$')
        ax.set_ylabel(r'$A - K_t G$')
        ax.set_xlabel('Time (t)')
        if save_fig:
            plt.savefig('./Ch2/p2_20b.eps', format='eps', dpi=1000)
        plt.show()

    return a_gk.max()



def trans(chain, a, b):
    n = 0
    for i in range(1, chain.size):
        if chain[i-1] == a and chain[i] == b:
            n += 1
    return n / float(chain.size - 1)

markov_17, pi_0 = p2_17(sym.Rational(1, 20), sym.Rational(1, 4))
