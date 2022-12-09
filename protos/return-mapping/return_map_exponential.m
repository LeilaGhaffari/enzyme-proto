% Return mapping algorithm - exponential
% Test with:
%   meld q_stress_src.txt q_stress.log; meld q_Fp_src.txt q_Fp.log ; meld q_isv_src.txt q_isv.log

clear all

coordsx = [ 0.0  0.0  0.0;
            1.0  0.0  0.0;
            1.0  1.0  0.0;
            0.0  1.0  0.0;
            0.0  0.0  1.0;
            1.0  0.0  1.0;
            1.0  1.0  1.0;
            0.0  1.0  1.0 ];

d = [ 0             0            0          ...
     -4.0123e-04    1.5702e-03   0          ...
     -1.9717e-03    1.3689e-03   0          ...
     -1.5705e-03   -2.0123e-04   0          ...
      0             0            1.9148e-04 ...
     -4.0123e-04    1.5702e-03   1.9148e-04 ...
     -1.9717e-03    1.3689e-03   1.9148e-04 ...
     -1.5705e-03   -2.0123e-04   1.9148e-04 ];

params = [ 1.5e+09 1.6e+09 1.633 0 1.633 0 306.19 8 13 5 3 ];

c_tau_n = 2.5e6;

% Fp
Fp_el_n = eye(3);

% rotation angle
theta = pi / 2;

nstress = round(params(9));
ndim    = round(params(11));

[stress_el, isv_el, Fp_el] = element_stress_isv(coordsx, d, params, c_tau_n, Fp_el_n, theta);

% Test the algorithm
stress_q = zeros(ndim*nstress, ndim);
for i=1:nstress
    stress_q((i-1)*ndim+1 : i*ndim, :) = stress_el(:, :, i);
end
isv_q = isv_el';
Fp_q = Fp_el;

save q_stress.log stress_q -ASCII
save    q_isv.log    isv_q -ASCII
save     q_Fp.log     Fp_q -ASCII

% -------------------------------------------------------------------------------------------------
%                                  Return mapping function
% -------------------------------------------------------------------------------------------------
function [stress_el,isv_el,Fp_el] = element_stress_isv(coordsx,d,params,c_tau_n,Fp_el_n,theta)

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
    isv_el    = zeros(nisv);
    Fp_el     = zeros(ndim, ndim);
    stress_el = zeros(ndim, ndim, nstress);

    % integration points in natural coordinates (this is for the fifth quadrature point)
    X = [-1 -1  1] / sqrt(3);

    % gradient of shape functions with respect to the natural coordinates
    dNdX = grad_basis_linear(X);

    % calculate jacobian
    Je = dNdX * coordsx;

    % gradient of shape functions with respect to the physical coordinates
    dN_dx = dNdX' / Je;

    % strain-displacement matrix
    Bu = zeros(9, 24);
    for k=1:8
        Bu(:, (k-1)*3+1:k*3) = kron(eye(3), dN_dx(k,:)');
    end

    % total deformation tensors
    dudX = Bu * d';
    Fdef = eye(3) + reshape(dudX, [3, 3])';
    Jdef = det(Fdef);
    Cdef = Fdef' * Fdef;
    bdef = Fdef  * Fdef';
    vdef = sqrtm(bdef);

    % total strains
    Estrain  = (Cdef - eye_mat) / 2;
    estrain  = (eye_mat - inv(bdef)) / 2;
    eHstrain = logm(vdef);

    % total principal Hencky strain
    lambda_eig     = eig(eHstrain);
    eHstrain_princ = [ max(lambda_eig)  0              0 ;
                       0                lambda_eig(2)  0 ;
                       0                0              min(lambda_eig) ];

    % calculate trial elastic deformation gradient and left elastic Cauchy Green tensor
    Fp_n  = Fp_el_n(:, :);
    Fe_tr = Fdef / Fp_n;
    be_tr = Fe_tr * Fe_tr';

    % find eigenvalues and eigenvectors of be_tr
    [V, D]      = eig(be_tr);
    lambda_e_tr = sqrt(diag(D));
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

    % dfdtau_tr
    dfdtau_princ_tr = dev_tau_princ_tr / dev_tau_tr_norm + Bphi / 3;
    trace_dfdtau_tr = sum(dfdtau_princ_tr);
    dfdtau_tr       = coaxial(be_tr, dfdtau_princ_tr);;

    % check if elastic, or elasto-plastic
    if (f_tr < 0)  % elastic
        Dg    = 0; % increment of plastic multiplier
        tau   = tau_tr;
        Fp    = Fp_n;
        c_tau = c_tau_n;
        be    = be_tr;
    else % plastic
        lambda_e = lambda_e_tr;
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
            k = k+1;

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
            del_Dg  = -fyield / dfdDg;
            Dg_iter = Dg_iter + del_Dg;
            if Dg_iter < 0
                Dg_iter
            end

            % update elastic stretches
            for i=1:3
                lambda_e(i) = lambda_e_tr(i) * exp(-Dg_iter * dgdtau_princ_tr(i));
            end
            Jdef_e = prod(lambda_e);

            % calculate Kirchhoff stress and deviatoric part
            for i=1:3
                tau_princ(i) = tau_princ_coef(Jdef_e, lambda_e(i));
            end
            p_tau        = sum(tau_princ) / 3;
            dev_tau      = tau_princ - p_tau;
            dev_tau_norm = norm(dev_tau);

            % update c_tau
            c_tau = c_tau_n + Dg_iter * Hc * Apsi;

            % update fyield
            fyield = dev_tau_norm - Aphi * c_tau + Bphi * p_tau;

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

        % update be
        Fe = Fdef * inv(Fp);
        be = Fe * Fe';
    end

    % stresses
    sigma   = tau / Jdef;            % Cauchy stress => \sigma = J^{-1} \tau
    Pstress = tau / Fdef';           % First Piola   => P = \tau F^{-T}
    Sstress = inv(Fdef) * Pstress;   % Second Piola  => S = F^{-1} P

    % rotate stress for non-rotating case
    Qrot = [  cos(theta) sin(theta) 0;
             -sin(theta) cos(theta) 0;
              0          0          1 ];
    sigma_rot = Qrot * sigma * Qrot';

    % principal Cauchy stress
    lambda_eig  = eig(sigma);
    sigma_princ = [ max(lambda_eig)  0              0 ;
                    0                lambda_eig(2)  0 ;
                    0                0              min(lambda_eig) ];

    % order of storage: S,P,sig,E,e,le,tau,be,dfdtau,dgdtau,sig_princ,sig_rot,le_princ
    stress_el(:, :,  1) = Sstress;
    stress_el(:, :,  2) = Pstress;
    stress_el(:, :,  3) = sigma;
    stress_el(:, :,  4) = Estrain;
    stress_el(:, :,  5) = estrain;
    stress_el(:, :,  6) = eHstrain;
    stress_el(:, :,  7) = tau;
    stress_el(:, :,  8) = be;
    stress_el(:, :,  9) = dfdtau_tr;
    stress_el(:, :, 10) = dgdtau_tr;
    stress_el(:, :, 11) = sigma_princ;
    stress_el(:, :, 12) = sigma_rot;
    stress_el(:, :, 13) = eHstrain_princ;

    % calculate stress invariants for output
    sig_mean = trace(sigma) / 3;
    sig_dev  = sigma - sig_mean * eye_mat;
    sig_VM   = sqrt(3/2) * norm(sig_dev);

    % order of storage: p, q, c_tau, Dg, Jdef
    isv_el = [ sig_mean sig_VM c_tau Dg Jdef ];

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
