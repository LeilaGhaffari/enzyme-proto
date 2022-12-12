% Return mapping algorithm - exponential
% Test with:
%   meld q_stress_src.txt q_stress.log; meld q_Fp_src.txt q_Fp.log ; meld q_isv_src.txt q_isv.log

clear all

% plastic
d = [   0	                  0	                     0                    ...
       -0.00121108996046273   0.00470671669352617	 0                    ...
       -0.00592063407691252	  0.00409562007109281	 0                    ...
       -0.00470954411644980  -0.000611096622433368	 0                    ...
        0                 	  0	                     0.000574665721779025 ...
       -0.00121108996046273	  0.00470671669352617	 0.000574665721779095 ...
       -0.00592063407691252	  0.00409562007109281	 0.000574665721779159 ...
       -0.00470954411644980  -0.000611096622433368	 0.000574665721779107 ];

% elastic
%d = [   0	                  0	                     0 ...
%       -0.000401233206816348  0.00157016736255855	 0 ...
%       -0.00197171472851106	  0.00136893390900214	 0 ...
%       -0.00157048152169471	 -0.000201233453556407	 0 ...
%        0	                  0	                     0 ...
%       -0.000401233206816348  0.00157016736255855	 0 ...
%       -0.00197171472851106	  0.00136893390900214	 0 ...
%       -0.00157048152169471  -0.000201233453556407	 0 ];

coordsx = [ 0.0  0.0  0.0;
            1.0  0.0  0.0;
            1.0  1.0  0.0;
            0.0  1.0  0.0;
            0.0  0.0  1.0;
            1.0  0.0  1.0;
            1.0  1.0  1.0;
            0.0  1.0  1.0 ];

params = [ 1.5e+09 1.6e+09 1.633 0 1.633 0 306.19 8 2 2 3 ];

c_tau_n = 2.5e6;

% Fp
Fp_el_n = eye(3);

nstress = round(params(9));
ndim    = round(params(11));

[stresses, isv, Fp_el] = return_map_exp(coordsx, d, params, c_tau_n, Fp_el_n);

% Test the algorithm
stress_q = zeros(ndim*nstress, ndim);
for i=1:nstress
    stress_q((i-1)*ndim+1 : i*ndim, :) = stresses(:, :, i);
end
isv_q = isv';
Fp_q = Fp_el;

save q_stress.log stress_q -ASCII
save    q_isv.log    isv_q -ASCII
save     q_Fp.log     Fp_q -ASCII

% -------------------------------------------------------------------------------------------------
%                                  Return mapping function
% -------------------------------------------------------------------------------------------------
function [stresses, isv, Fp_el] = return_map_exp(coordsx, d, params, c_tau_n, Fp_n)

    % Helper function: 3x3 identity matrix
    eye_mat = eye(3);

    % params = [ lambda mu Aphi Bphi Apsi Bpsi Hc numips nstress nisv ndim ];
    lambda  = params(1);
    mu      = params(2);
    Aphi    = params(3);
    Bphi    = params(4);
    Apsi    = params(5);
    Bpsi    = params(6);
    Hc      = params(7);
    nstress = round(params(9));
    nisv    = round(params(10));
    ndim    = round(params(11));

    % initialize
    isv      = zeros(nisv);
    Fp_el    = zeros(ndim, ndim);
    stresses = zeros(ndim, ndim, nstress);

    % Mapping from natural coordinates to physical coordinates (X -> x)
    X = [-1 -1  1] / sqrt(3);
    dNdX = grad_basis_linear(X);
    Je = dNdX * coordsx;
    dN_dx = dNdX' / Je;

    % strain-displacement matrix
    Bu = zeros(9, 24);
    for k=1:8
        Bu(:, (k-1)*3+1:k*3) = kron(eye(3), dN_dx(k,:)');
    end

    % total deformation tensors
    dudX = Bu * d';
    Fdef = eye(3) + reshape(dudX, [3, 3])';
    bdef = Fdef  * Fdef';

    % total strain
    estrain  = (eye_mat - inv(bdef)) / 2;

    % calculate trial elastic deformation gradient and left elastic Cauchy Green tensor
    Fe_tr = Fdef / Fp_n;
    be_tr = Fe_tr * Fe_tr';

    % find eigenvalues and eigenvectors of be_tr
    [V, D]      = eig(be_tr);
    lambda_e_tr = sqrt(diag(D))
    Jdef_e_tr   = prod(lambda_e_tr);

    % calculate trial Kirchhoff stress
    tau_princ_tr   = zeros(3, 1);
    tau_princ_coef = @(J, eigval) lambda * log(J) - mu + mu * eigval * eigval;
    for i=1:3
        tau_princ_tr(i) = tau_princ_coef(Jdef_e_tr, lambda_e_tr(i));
    end
    tau_tr           = coaxial(be_tr, tau_princ_tr);
    p_tau_tr         = sum(tau_princ_tr) / 3;
    dev_tau_princ_tr = tau_princ_tr - p_tau_tr;

    % calculate trial elastic left stretch and rotation
    ve_tr = sqrtm(be_tr);
    Re_tr = inv(ve_tr) * Fe_tr;

    % calculate trial yield function
    c_tau           = c_tau_n;
    dev_tau_tr_norm = norm(dev_tau_princ_tr);
    f_tr            = dev_tau_tr_norm - Aphi * c_tau_n + Bphi * p_tau_tr;

    % dgdtau_tr
    dgdtau_princ_tr = dev_tau_princ_tr / dev_tau_tr_norm + Bpsi / 3;
    trace_dgdtau_tr = sum(dgdtau_princ_tr);
    dgdtau_tr       = coaxial(be_tr, dgdtau_princ_tr);

    % check if elastic, or elasto-plastic
    if (f_tr < 0)  % elastic
        Dg    = 0; % increment of plastic multiplier
        tau   = tau_tr;
        Fp    = Fp_n;
        c_tau = c_tau_n;
        be    = be_tr;
    else % plastic
        lambda_e = lambda_e_tr
        Jdef_e   = prod(lambda_e);

        % calculate Kirchhoff stress and deviatoric part
        tau_princ = zeros(3, 1);
        for i=1:3
            tau_princ(i) = tau_princ_coef(Jdef_e, lambda_e(i));
        end
        p_tau         = sum(tau_princ) / 3;
        dev_tau_princ = tau_princ - p_tau;
        dev_tau_norm  = norm(dev_tau_princ);

        % solve for Dg (Delta_gamma)
        fyield_atol      = 1e-6;
        fyield_rtol      = 1e-8;
        k                = 0;
        iter_break_local = 10;
        fyield           = f_tr;
        Dg_iter          = 0;
        while ( abs(fyield) > fyield_atol && abs(fyield/f_tr) > fyield_rtol )
            k = k+1

            % derivative of devtau wrt Dg
            dtau_princ_dDg = zeros(3, 1);
            for i=1:3
                dtau_princ_dDg(i) = -lambda * trace_dgdtau_tr - 2 * mu * lambda_e(i) * lambda_e(i) * dgdtau_princ_tr(i);
            end
            trace_dtaudDg      = sum(dtau_princ_dDg);
            dev_dtau_princ_dDg = dtau_princ_dDg - trace_dtaudDg / 3;

            % derivative of devtau_norm wrt Dg
            ddevtau_norm_dDg = dot(dev_tau_princ, dev_dtau_princ_dDg) / dev_tau_norm;

            % yield tangent
            dfdDg = ddevtau_norm_dDg - Hc * Aphi * Apsi + Bphi * trace_dtaudDg / 3;

            % increment of Dg
            del_Dg  = -fyield / dfdDg
            Dg_iter = Dg_iter + del_Dg

            % update elastic stretches
            for i=1:3
                lambda_e(i) = lambda_e_tr(i) * exp(-Dg_iter * dgdtau_princ_tr(i));
            end
            lambda_e
            Jdef_e = prod(lambda_e);

            % calculate Kirchhoff stress and deviatoric part
            for i=1:3
                tau_princ(i) = tau_princ_coef(Jdef_e, lambda_e(i));
            end
            tau_princ
            p_tau        = sum(tau_princ) / 3;
            dev_tau      = tau_princ - p_tau;
            dev_tau_norm = norm(dev_tau);

            % update c_tau
            c_tau = c_tau_n + Dg_iter * Hc * Apsi;

            % update fyield
            dev_tau_norm;
            Aphi_mult_ctau = Aphi * c_tau;
            Bphi_mult_p_tau = Bphi * p_tau;
            fyield = dev_tau_norm - Aphi_mult_ctau + Bphi * p_tau;


            if k == iter_break_local
                k
                error('reached max number of local iterations for Dg')
            end
        end

        % update Fp
        Fp = Re_tr' * expm(Dg_iter * dgdtau_tr) * Re_tr * Fp_n;

        % update tau
        tau = coaxial(be_tr, tau_princ);

        % update Dg
        Dg = Dg_iter;
    end

    stresses(:, :, 1) = estrain;
    stresses(:, :, 2) = tau;

    % Store isv
    isv = [c_tau Dg];

    % return Fp
    Fp_el = Fp;
end


function C = coaxial(A, b)
    [eig_vec, ~] = eig(A);
    C            = zeros(3, 3);
    for i=1:3
        C = C + b(i) * eig_vec(:, i) * eig_vec(:, i)';
    end
end


function dNdX = grad_basis_linear(xi)
    % derivatives of shape functions with respect to xi
    dN1_dxi = -0.125*(1-xi(2))*(1-xi(3));
    dN2_dxi = -dN1_dxi;
    dN3_dxi = 0.125*(1+xi(2))*(1-xi(3));
    dN4_dxi = -dN3_dxi;
    dN5_dxi = -0.125*(1-xi(2))*(1+xi(3));
    dN6_dxi = -dN5_dxi;
    dN7_dxi = 0.125*(1+xi(2))*(1+xi(3));
    dN8_dxi = -dN7_dxi;

    % derivatives of shape functions with respect to eta
    dN1_deta = -0.125*(1-xi(1))*(1-xi(3));
    dN2_deta = -0.125*(1+xi(1))*(1-xi(3));
    dN3_deta = -dN2_deta;
    dN4_deta = -dN1_deta;
    dN5_deta = -0.125*(1-xi(1))*(1+xi(3));
    dN6_deta = -0.125*(1+xi(1))*(1+xi(3));
    dN7_deta = -dN6_deta;
    dN8_deta = -dN5_deta;

    % derivatives of shape functions with respect to zeta
    dN1_dzeta = -0.125*(1-xi(1))*(1-xi(2));
    dN2_dzeta = -0.125*(1+xi(1))*(1-xi(2));
    dN3_dzeta = -0.125*(1+xi(1))*(1+xi(2));
    dN4_dzeta = -0.125*(1-xi(1))*(1+xi(2));
    dN5_dzeta = -dN1_dzeta;
    dN6_dzeta = -dN2_dzeta;
    dN7_dzeta = -dN3_dzeta;
    dN8_dzeta = -dN4_dzeta;

    % Populate dNdX
    dNdX = [dN1_dxi dN2_dxi dN3_dxi dN4_dxi dN5_dxi dN6_dxi dN7_dxi dN8_dxi;
            dN1_deta dN2_deta dN3_deta dN4_deta dN5_deta dN6_deta dN7_deta dN8_deta;
            dN1_dzeta dN2_dzeta dN3_dzeta dN4_dzeta dN5_dzeta dN6_dzeta dN7_dzeta dN8_dzeta];
end
