# Mock the exponential return mapping algorithm from CVEN7511
using LinearAlgebra

function coaxial(A, b)
    eig_vec = eigvecs(A)
    C = zeros(3, 3)
    for i = 1:3
        C[:, :] = C[:, :] + b[i] * ((eig_vec[:, i] * eig_vec[:, i]'))
    end
    return C
end

function grad_basis_linear(xi)
    # derivatives of shape functions with respect to xi
    dN1_dxi = -0.125 * (1 - xi[2]) * (1 - xi[3])
    dN2_dxi = -dN1_dxi
    dN3_dxi = 0.125 * (1 + xi[2]) * (1 - xi[3])
    dN4_dxi = -dN3_dxi
    dN5_dxi = -0.125 * (1 - xi[2]) * (1 + xi[3])
    dN6_dxi = -dN5_dxi
    dN7_dxi = 0.125 * (1 + xi[2]) * (1 + xi[3])
    dN8_dxi = -dN7_dxi

    # derivatives of shape functions with respect to eta
    dN1_deta = -0.125 * (1 - xi[1]) * (1 - xi[3])
    dN2_deta = -0.125 * (1 + xi[1]) * (1 - xi[3])
    dN3_deta = -dN2_deta
    dN4_deta = -dN1_deta
    dN5_deta = -0.125 * (1 - xi[1]) * (1 + xi[3])
    dN6_deta = -0.125 * (1 + xi[1]) * (1 + xi[3])
    dN7_deta = -dN6_deta
    dN8_deta = -dN5_deta

    # derivatives of shape functions with respect to zeta
    dN1_dzeta = -0.125 * (1 - xi[1]) * (1 - xi[2])
    dN2_dzeta = -0.125 * (1 + xi[1]) * (1 - xi[2])
    dN3_dzeta = -0.125 * (1 + xi[1]) * (1 + xi[2])
    dN4_dzeta = -0.125 * (1 - xi[1]) * (1 + xi[2])
    dN5_dzeta = -dN1_dzeta
    dN6_dzeta = -dN2_dzeta
    dN7_dzeta = -dN3_dzeta
    dN8_dzeta = -dN4_dzeta

    # Populate dNdX
    dNdX = [dN1_dxi dN2_dxi dN3_dxi dN4_dxi dN5_dxi dN6_dxi dN7_dxi dN8_dxi
        dN1_deta dN2_deta dN3_deta dN4_deta dN5_deta dN6_deta dN7_deta dN8_deta
        dN1_dzeta dN2_dzeta dN3_dzeta dN4_dzeta dN5_dzeta dN6_dzeta dN7_dzeta dN8_dzeta]
    return dNdX
end


# -------------------------------------------------------------------------------------------------
#                                  Return mapping function
# -------------------------------------------------------------------------------------------------
function return_map_exp(coordsx, d, params, c_tau_n, Fp_n)
    # N-R specifications
    fyield_atol = 1e-6
    fyield_rtol = 1e-8
    iter_break_local = 10

    # params = [ lambda mu Aphi Bphi Apsi Bpsi Hc ];
    lambda = params[1]
    mu = params[2]
    Aphi = params[3]
    Bphi = params[4]
    Apsi = params[5]
    Bpsi = params[6]
    Hc = params[7]

    # define principal Kirchhoff stress
    tau_princ_coef(J, eigval) = lambda * log(J) - mu + mu * eigval * eigval

    # Mapping from natural coordinates to physical coordinates (X -> x)
    X = [-1 -1 1] / sqrt(3)
    dNdX = grad_basis_linear(X)
    Je = dNdX * coordsx
    dN_dx = transpose(dNdX) / Je

    # strain-displacement matrix
    Bu = zeros(9, 24)
    for k = 1:8
        Bu[:, (k-1)*3+1:k*3] = kron(I(3), dN_dx[k, :])
    end

    # total deformation tensors
    dudX = Bu * d
    Fdef = I(3) + transpose(reshape(dudX, (3, 3)))

    # calculate trial elastic deformation gradient and left elastic Cauchy Green tensor
    Fe_tr = Fdef / Fp_n
    be_tr = Fe_tr * Fe_tr'

    # find eigenvalues and eigenvectors of be_tr
    D = eigvals(be_tr)
    # D = D[[1, 3, 2]]
    lambda_e_tr = sqrt.(D)
    Jdef_e_tr = prod(lambda_e_tr)

    @show lambda_e_tr

    # calculate trial Kirchhoff stress
    tau_princ_tr = zeros(3)
    for i = 1:3
        tau_princ_tr[i] = tau_princ_coef(Jdef_e_tr, lambda_e_tr[i])
    end
    tau_tr = coaxial(be_tr, tau_princ_tr)
    p_tau_tr = sum(tau_princ_tr) / 3
    dev_tau_princ_tr = tau_princ_tr .- p_tau_tr

    # calculate trial elastic left stretch and rotation
    ve_tr = sqrt(be_tr)
    Re_tr = inv(ve_tr) * Fe_tr

    # calculate trial yield function
    c_tau = c_tau_n
    dev_tau_tr_norm = norm(dev_tau_princ_tr)
    f_tr = dev_tau_tr_norm - Aphi * c_tau_n + Bphi * p_tau_tr

    # dgdtau_tr
    dgdtau_princ_tr = dev_tau_princ_tr ./ dev_tau_tr_norm .+ Bpsi / 3
    trace_dgdtau_tr = sum(dgdtau_princ_tr)
    dgdtau_tr = coaxial(be_tr, dgdtau_princ_tr)

    # check if elastic, or elasto-plastic
    if (f_tr < 0)  # elastic
        Dg = 0 # increment of plastic multiplier
        tau = tau_tr
        Fp = Fp_n
        c_tau = c_tau_n
    else # plastic
        lambda_e = lambda_e_tr
        Jdef_e = prod(lambda_e)
        @show lambda_e

        # calculate Kirchhoff stress and deviatoric part
        tau_princ = zeros(3)
        for i = 1:3
            tau_princ[i] = tau_princ_coef(Jdef_e, lambda_e[i])
        end
        p_tau = sum(tau_princ) / 3
        dev_tau_princ = tau_princ .- p_tau
        dev_tau_norm = norm(dev_tau_princ)

        # solve for Dg (Delta_gamma)
        k = 0
        Dg_iter = 0
        fyield = f_tr
        while (abs(fyield) > fyield_atol && abs(fyield / f_tr) > fyield_rtol)
            # update iterator
            println(" ")
            k += 1
            @show k

            # derivative of devtau wrt Dg
            dtau_princ_dDg = zeros(3)
            for i = 1:3
                dtau_princ_dDg[i] = -lambda * trace_dgdtau_tr - 2 * mu * lambda_e[i] * lambda_e[i] * dgdtau_princ_tr[i]
            end
            trace_dtaudDg = sum(dtau_princ_dDg)
            dev_dtau_princ_dDg = dtau_princ_dDg .- trace_dtaudDg / 3

            # derivative of devtau_norm wrt Dg
            ddevtau_norm_dDg = dot(dev_tau_princ, dev_dtau_princ_dDg) / dev_tau_norm

            # yield tangent
            dfdDg = ddevtau_norm_dDg - Hc * Aphi * Apsi + Bphi * trace_dtaudDg / 3

            # increment of Dg
            del_Dg = -fyield / dfdDg
            Dg_iter = Dg_iter + del_Dg
            @show del_Dg
            @show Dg_iter

            # update elastic stretches
            for i = 1:3
                lambda_e[i] = lambda_e_tr[i] * exp(-Dg_iter * dgdtau_princ_tr[i])
            end
            Jdef_e = prod(lambda_e)

            @show lambda_e

            # calculate Kirchhoff stress and deviatoric part
            for i = 1:3
                tau_princ[i] = tau_princ_coef(Jdef_e, lambda_e[i])
            end
            p_tau = sum(tau_princ) / 3
            dev_tau = tau_princ .- p_tau
            dev_tau_norm = norm(dev_tau)

            @show tau_princ

            # update c_tau
            c_tau = c_tau_n + Dg_iter * Hc * Apsi

            # update fyield
            fyield = dev_tau_norm - (Aphi * c_tau + Bphi * p_tau)

            if k == iter_break_local
                println(" ")
                println("Reached max number of iterations (", k, ")")
                break
            end
        end

        # update Fp
        Fp = transpose(Re_tr) * exp(Dg_iter * dgdtau_tr) * Re_tr * Fp_n

        # update tau
        tau = coaxial(be_tr, tau_princ)

        # update Dg
        Dg = Dg_iter
    end

    return Fp, tau, c_tau, Dg
end

# -------------------------------------------------------------------------------------------------
#  Test the implementation
# -------------------------------------------------------------------------------------------------
# elastic
#d = [0, 0, 0,
#    -4.0123e-04, 1.5702e-03, 0,
#    -1.9717e-03, 1.3689e-03, 0,
#    -1.5705e-03, -2.0123e-04, 0,
#    0, 0, 1.9148e-04,
#    -4.0123e-04, 1.5702e-03, 1.9148e-04,
#    -1.9717e-03, 1.3689e-03, 1.9148e-04,
#    -1.5705e-03, -2.0123e-04, 1.9148e-04];
# plastic
d = [0, 0, 0,
    -0.00121108996046273, 0.004706716693526170, 0,
    -0.00592063407691252, 0.004095620071092810, 0,
    -0.00470954411644980, -0.000611096622433368, 0,
    0, 0, 0.000574665721779025,
    -0.00121108996046273, 0.004706716693526170, 0.000574665721779095,
    -0.00592063407691252, 0.004095620071092810, 0.000574665721779159,
    -0.00470954411644980, -0.000611096622433368, 0.000574665721779107];

coordsx = [
    0.0 0.0 0.0
    1.0 0.0 0.0
    1.0 1.0 0.0
    0.0 1.0 0.0
    0.0 0.0 1.0
    1.0 0.0 1.0
    1.0 1.0 1.0
    0.0 1.0 1.0
];

params = [1.5e+09 1.6e+09 1.633 0 1.633 0 306.19];

c_tau_n = 2.5e6;

# Fp
Fp_el_n = I(3);

Fp, tau, c_tau, Dg = return_map_exp(coordsx, d, params, c_tau_n, Fp_el_n)
