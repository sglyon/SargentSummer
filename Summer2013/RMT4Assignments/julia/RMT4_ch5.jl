require("rmt_utils.jl")

function p5_3(beta, r, b, gamma, rho_1, rho_2)
    #
    A = [rho_1 rho_2 0 0
         0 1 0 0
         0 0 1 0
         0 0 0 1]

    B = [0.
         0.
         1.
         0.]

    R = [1 0 r -b
         0 0 0 0
         r 0 r^2 -b * r
         -b 0 -b*r b^2]

    Q = 1 + gamma

    H = [-1
         0
         -r
         b]

    return A, B, R, Q, H
end

in_53 = [0.95, 1 / 0.95 - 1., 30, 1., 1.2, -0.3]
A, B, R, Q, H = in_dblo = p5_3(in_53...)
