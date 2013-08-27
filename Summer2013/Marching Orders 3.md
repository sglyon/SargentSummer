# RMT4 Chapter 5 - Linear Quadratic Dynamic Programming {#ch5}

This chapter will talk about dynamic programming problems where the return functions are quadratic and the transition functions are linear.

## The optimal linear regulator problem {#optimlinreg}

The undiscounted optimal linear regulator problem is to maximize over choice of $m\{u_t\}_{t=0}^{\infty}$ the criterion

$$ - \sum_{t=0}^{\infty} \{x_t' R x_t + u_t' Q u_t \} $$


We guess that the value function is quadratic, $V(x) = - x'Px$. Some algebra covered in the book yields the algebraic matrix Riccati equation and the solution for $P$:

$$P = R + A'PA - A'RB(Q + B'PB)^{-1}B'PA$$

### Value function iteration {#vfi}

In specific cases there is a unique positive definite matrix $P$ that defines the solution to the matrix Riccati equation above. This is found by iterating on that equation by saying

$$P_{j+1} = R + A'P_jA - A'RB(Q + B'P_jB)^{-1}B'P_jA$$

Starting at $P_0 = 0$.

### Discounted linear regulator problem {#dislinreg}

The discounted optimal linear regulator problem and its matrix Riccati difference equation appear below

$$ - \sum_{t=0}^{\infty} \beta_t \{x_t' R x_t + u_t' Q u_t \}, \quad 0 < \beta <  1 $$

$$P_{j+1} = R + \beta A'P_jA - \beta^2 A'RB(Q + B'P_jB)^{-1}B'P_jA$$

The Matlab program `olrp.m` solves the discounted optimal linear regulator problem.

### Policy improvement algorithm {#policyalgo}

The policy improvement algorithm (Howards' improvement algorithm) can be applied by iterating on the following two equations, assuming we have a starting $F_0$ and know that all the eigenvalues of $A - BF_0$ are bounded below $\frac{1}{\sqrt{\beta}}$ in modulus:

$$P_j = R + F_j'Qf_j + \beta (A - BF_j)'P(A-BF_j)$$

$$F_{j+1} = \beta(Q + \beta B'P_jB)^{-1}B'P_jA$$

The equation for $P_j$ above is a discrete Lyapunov or Sylvester equation and can be solved for $P_j$ using the function `doublej`. This $P_j$ gives the associated value ($V(x) = -x'P_jx$) of applying the policy $F_j$ forever. The solution to this system can be represented as

$$P_j =\sum_{k=0}^{\infty} \beta^k (A - B F_j)^{'k}(R + F_j'QF_j) (A - BF_j)^k$$

As long as the eigenvalues of $A - BF_j$ are bounded below $\frac{1}{\sqrt{\beta}}$ in modulus, then a solution exists. The function `policyi` solves the un-discounted optimal linear regulator problem using policy iteration.


## The stochastic optimal linear regulator problem {#stochlinreg}

The stochastic discounted optimal linear regulator problem is to choose a decision rule for $u_t$ to maximize

$$ - E_0\sum_{t=0}^{\infty} \beta_t \{x_t' R x_t + u_t' Q u_t \}, \quad 0 < \beta <  1$$

subject to a given $x_0$ and the following law of motion for x:

$$x_{t+1} = Ax_t + Bu_t  + C \epsilon_{t+1}$$

where $\epsilon_{t+1} \sim N(0, I)$ is an $(n \times 1)$ vector of iid random variables.

The new value function for this problem is

$$ v(x) = -x'Px - d $$

Where $P$ is the equation from the [discounted section](#dislinreg) and the scalar $d$ is given by

$$d = \beta (1 - \beta)^{-1} \text{trace}(PCC')$$

The value of $u_t = - F x_t$ is still defined by the matrix $F$ from the [policy algorithm section](#policyalgo).

The **Certainty Equivalence Principle** holds here because the solutions $F$ for the stochastic and non-stochastic discounted optimal linear regulator problems are the same (stochastic part doesn't change the policy).

### Discussion of certainty equivalence {#certiantyequiv}

It is quite remarkable that although the objective function depends on $CC'$ (through $d$), the solution is independent of $CC'$. Another way to think of this is that the optimal decision rule $u_t = h(x_t)$ is *independent* of the problem's noise statistics.

Certainty equivalence is a special property we get here because we designed the problem as such; certainty equivalence is *not* generally true of stochastic control problems!

## Shadow prices in the linear regulator {#shadows}

The gradient of the value function (either $v(x) = -x'Px -d$ or $v(x) = - x'Px$) is $-2Px$. This gradient can be interpreted as a shadow price or Lagrangian multiplier. To show this, we will express the Bellman function as a Lagrangian:

$$ - x_t'Px_t = V(x_t) -= \min_{\mu_{t + 1}} \max_{u_t, x_{t+1}} - \left\{x_t' R x_t  + u_t'Q u_t + x_{t+1}' P x_{t+1} + 2 \mu_{t+1}' \left[A x_t + B u_t  -x_{t+1} \right] \right\}$$

where $2 \mu_{t+1}$ is a vector of Lagrange multipliers. The first order necessary conditions for an optimum with respect to $u_t$ and $x_{t+1}$ are, respectively

$$
\begin{align}
2 Q u_t + 2 B' \mu_{t+1} &= 0\\
2 P x_{t+1} - 2 \mu_{t+1}  &= 0
\end{align}
$$

This implies the following value for the shadow price vector:

$$\mu_{t+1} = P x_{t+1}$$

### Stability {#stable}

When the optimal control selection $u_t = - F x_t$ is substituted into the law of motion $x_{t+1} = A x_t + B u_t$, we can obtain the optimal *closed loop system* represented as $x_{t+1} = (A - BF) x_t$. This governs the evolution of $x_t$.

As explained in this section of the book, the matrix $(A - BF)$ is considered stable if all its eigenvalues are strictly less than unity in absolute value. When that matrix is stable, $\lim_{x \rightarrow \infty} x_t = 0$.

A case optimal control literature is dedicated to pinning down conditions for matrices $A, B,Q$ and $R$ so that $F$ leads to an optimal control closed-loop system (See @sargent1996 for more info).


## A Lagrangian formulation

Formulating the optimal linear regulator in this way has computational benefits as well as the ability to make connections between stability and optimality.

For the undiscounted problem, the Lagrangian is

$$ \mathscr{L} = - \sum_{t = 0}^{\infty} \left\{x_t' R x_t + u_t' Q u_t + 2 \mu_{t+1} \left[A x_t  + B u_t - x_{t+1} \right] \right\}$$

The first order conditions with respect to $u_t$ and $x_{t+1}$ are

$$
\begin{align}
0 &= 2 Q u_t + 2 B' \mu_{t+1} \\
\mu_t &= R x_t + A'\mu_{t+1}, \quad t \ge 0
\end{align}
$$

The Lagrange multiplier vector $\mu_{t+1}$ is often called the costate vector.

We can use the expression for $\mu_{t+1}$ that we obtained in the [shadow pricing section](#shadows): $\mu_{t+1} = P x_{t+1}$. We can then solve the equation above for $\mu_t$ in terms of $\mu_{t+1}$, substitute in the optimal condition and the law of motion for $x$ to get the system above into the following matrix form:

$$ L \begin{pmatrix}x_{t+1} \\ \mu_{t+1} \end{pmatrix} = N \begin{pmatrix}x_t \\ \mu_t \end{pmatrix}, \quad t \ge 0$$

where

$$L= \begin{pmatrix}I & B Q^{-1}B' \\ 0 & A' \end{pmatrix}, \quad N = \begin{pmatrix}A & 0 \\ -R & I \end{pmatrix}$$

When $L$ is invertible, we can write this system as

$$ \begin{pmatrix}x_{t+1} \\ \mu_{t+1} \end{pmatrix} = M \begin{pmatrix}x_t \\ \mu_t \end{pmatrix} $$

where $M$ is the $(2 \times 2)$ matrix:

$$ M = L^{-1} N = \begin{pmatrix} A + BQ^{-1}B'A^{'-1}R & -BQ^{-1}B'A^{'-1} \\ -A^{'-1}R & A^{'-1}\end{pmatrix}$$

I now introduce a new matrix $J$, which is of rank $2n$:

$$J = \begin{pmatrix}0 & -I_n \\ I_n & 0 \end{pmatrix}$$

The matrix $M$ is called *symplectic* if

$$MJM' = J$$

Our matrix $M$ is symplectic, which means that its eigenvalues come in reciprocal pairs (this comes from the above equation -- think about it). We can define a new variable $y_t = \left(\begin{smallmatrix} x_t \\ \mu_t\end{smallmatrix} \right)$ and write the following triangularization of $M$:

$$V^{-1}MV = \begin{pmatrix}W_{11} & W_{12} \\ 0 & W_{22} \end{pmatrix}$$

where each block on the RHS is and $(n \times n)$ matrix, $V$ is non-singular, all the eigenvalues of $W_{22}$ are greater than 1 in modulus and all eigenvalues of $W_{11}$ are less than 1 in modulus. The Schur and eigenvalue decomposition both follow these properties. Note also that the function `scipy.linalg.schur` provides the Shcur decomposition.

We can then write our first order conditions in terms of $y, V$, and $W$ as follows.

$$
\begin{align}
\begin{pmatrix}x_{t+1} \\ \mu_{t+1} \end{pmatrix} &= M \begin{pmatrix}x_t \\ \mu_t \end{pmatrix}\\
y_{t+1} &= M y_t \\
y_{t+1} &= V W V^{-1} y_t
\end{align}
$$

There is some more algebra in the book, but the result is that we can find the optimal value of $P$ via the expression

$$P + V{21}V_{11}^{-1}$$

This method of finding $P$ is generally very efficient computationally.


## The Kalman filter again {#kalman2}

As described in RMT4 chapter 2, the Kalman filter is a recursive algorithm for computing the mathematical expectation $E\left[x_t|y_{t+1}, \dots, y_0 \right]$ of a hidden state vector $x_t$ conditional on observing a history of a vector of noisy signals on the state $y_t$.

It turns out that the same recursion we used in the [value function iteration section](#vfi) to solve the optimal linear regulator also determines the Kalman filter. The required algebra to prove this is quite tedious, but it can be done.

The Matlab program `kfilter.m` computes the Kalman filter.

## Matrix equations {#matrixeq}

Let (z,x,a) each be $n \times 1$ vectors, A,C,D, and V each be ($n \times n$) matrices, B an ($m \times n$) matrix,and y an ($m \times 1$) vector. The the following derivatives are valid:

- $\frac{\partial a'x}{\partial x} = a$
- $\frac{\partial x'A x}{\partial x} = (A + A')x$
- $\frac{\partial ^2 x'A x}{\partial x \partial x'} = (A + A')$
- $\frac{\partial x'Ax}{\partial A} = xx'$
- $\frac{\partial y'Bz}{\partial y} = Bz$
- $\frac{\partial y'Bz}{\partial z} B'y$
- $\frac{\partial y'Bz}{\partial B} = yz'$

The **disrete Lyapuov equation** is defined as follows when trying to solve for $V$:

$$ A'V A + C = V$$

A generalization of this is the **discrete Sylvester equation**:

$$A'V D + C = V$$

The discrete Sylvester equation has a unique solution iff the eigenvalues $\{\lambda_i\}$ of $A$ and $\{\delta_i\}$ of $D$ satisfy the condition that $\lambda_i \delta_i \ne 1 \forall i, j$.

### Solving with `slycot`

`slycot` is a Python wrapper around the Fortran routines in `SLICOT` (Subroutine Library In COntrol Theory). `slycot` exposes highly optimized, stable, and numerically accurate Fortran routines for solving Riccati, Sylvester, and Lyapuov equations. I describe the pertinent routines here.

- `slycot.sb04qd`: This routine solve *discrete Sylvester equations*. The `SLICOT` docs say it solves for $X$ in the form $X + AXB = C$, where the `slycot` docstrings say the form is$AXB + X + C = 0$. I don't think it matters (TODO: verify this). Also, it is presented in a slightly different form than above. The equations above can be transformed into this form, but I am not sure whether or not that is necessary (TODO: verify this also).
- TODO: Document other routines for Lyapuov and Riccati equations

# References
