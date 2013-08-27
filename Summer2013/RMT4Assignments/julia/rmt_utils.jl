function doubleo(A, C, Q, R)

    a0=A'
    b0=C'*(R\C)
    g0=Q
    tol=1e-15
    dd=1
    ss=max(size(A))
    v=eye(ss)

    while dd > tol

        a1=a0 *((v+b0*g0)\a0)
        b1=b0+a0*((v+b0*g0)\(b0*a0'))
        g1=g0+a0'*g0*((v+b0*g0)\a0)
        k1=A*g1*C'/(C*g1*C'+R)
        k0=A*g0*C'/(C*g0*C'+R)
        dd=max(abs(k1-k0))
        a0=a1
        b0=b1
        g0=g1

    end

    return k1, g1

end


function olrp(beta, A, B, Q, R, W=nothing)
    #function [f,p] = olrp(beta,A,B,Q,R,W)
    #%OLRP can have arguments: (beta,A,B,Q,R) if there are no cross products
    #     (i.e. W=0).  Set beta=1, if there is no discounting.
    #
    #     OLRP calculates f of the feedback law:
    #
    #       u = -fx
    #
    #  that maximizes the function:
    #
    #          sum {beta^t [x'Qx + u'Ru +2x'Wu] }
    #
    #  subject to
    #       x[t+1] = Ax[t] + Bu[t]
    #
    #  where x is the nx1 vector of states, u is the kx1 vector of controls,
    #  A is nxn, B is nxk, Q is nxn, R is kxk, W is nxk.
    #
    #    Also returned is p, the steady-state solution to the associated
    #  discrete matrix Riccati equation.
    #

    m=max(size(A))

    if ndims(B) == 1
        rb = length(B)
        cb = 1
        cb1 = true
    else
        rb, cb = size(B)
        cb1 = false
    end

    if W == nothing
        W = zeros(m, cb)
    end

    if min(abs(eigvals(R))) > 1e-5

        A=sqrt(beta)*(A-B*(R\W'))
        B=sqrt(beta)*B
        Q=Q-W*(R\W')

        k, s = doubleo(A',B',Q,R)

        f=k'+(R\W')

        p=s

    else

        p0=-.01*eye(m)
        dd=1
        it=1
        maxit=1000

        # check tolerance; for greater accuracy set it to 1e-10

        while (dd>1e-6 && it<=maxit)

            f0_1 = cb1 ? (B'*p0*B)[1] : (B'*p0*B)

            f0=   (R+beta*f0_1)\(beta*B'*p0*A+W')
            p1=beta*A'*p0*A + Q -(beta*A'*p0*B+W)*f0'
            f1=   (R+beta*B'*p1*B)\(beta*B'*p1*A+W')
            dd=max(max(abs(f1-f0)))
            it=it+1
            p0=p1

        end
        f=f1
        p=p0
        if it>=maxit; warn("Iteration limit of 1000 reached in OLRP"); end

    end

    return f, p

end
