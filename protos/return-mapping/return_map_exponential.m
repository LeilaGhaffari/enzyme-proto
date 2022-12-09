% Return mapping algorithm - exponential
% Test with:
%   meld q_stress_src.txt q_stress.log; meld q_Fp_src.txt q_Fp.log ; meld q_isv_src.txt q_isv.log

clc
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

    %params = [ lambda mu Aphi Bphi Apsi Bpsi Hc numips nstress nisv ndim ];
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

    % tolerances for local Newton-Raphson algorithm to solve for Delta gamma
    fyield_atol = 1e-6;
    fyield_rtol = 1e-8;


    % initialize size of stress, and ISVs at Gauss pts, to be returned
    isv_el    = zeros(nisv);
    Fp_el     = zeros(ndim, ndim);
    stress_el = zeros(ndim, ndim, nstress);

    % integration points in natural coordinates (this is for the fifth quadrature point)
    xi   = -1 / sqrt(3);
    eta  = -1 / sqrt(3);
    zeta =  1 / sqrt(3);

    % gradient of shape functions
    dNdX = grad_basis_linear(xi, eta, zeta);

    % calculate jacobian, its determinant, and its inverse
    Je    = dNdX * coordsx;
    Jeinv = inv(Je);

    %shape function derivatives with respect to X,Y,Z
    dN_dx = dNdX' * Jeinv;

    % strain-displacement matrix
    Bu = zeros(9, 24);
    for k=1:8
        Bu(:, (k-1)*3+1:k*3) = kron(eye(3), dN_dx(k,:)');
    end

    % total deformation tensors
    dudX     = Bu*d';
    Fdef     = eye(3) + [ dudX(1) dudX(2) dudX(3); dudX(4) dudX(5) dudX(6); dudX(7) dudX(8) dudX(9) ];
    Jdef     = det(Fdef);
    Fdef_inv = inv(Fdef);
    Cdef     = Fdef' * Fdef;
    Cdef_inv = inv(Cdef);
    bdef     = Fdef * Fdef';
    bdef_inv = inv(bdef);
    vdef     = sqrtm(bdef);

    % 3x3 identity matrix
    eye_mat = eye(3);

    % total strains
    Estrain  = (Cdef-eye_mat)/2;
    estrain  = (eye_mat-inv(bdef))/2;
    eHstrain = logm(vdef);

    % total principal Hencky strain
    lambda_eig     = eig(eHstrain);
    eHstrain_princ = [ max(lambda_eig)  0              0 ;
                           0                lambda_eig(2)  0 ;
                           0                0              min(lambda_eig) ];

    % calculate trial elastic deformation gradient and left elastic Cauchy Green tensor
    Fp_n      = Fp_el_n(:, :);
    Fp_n_inv  = inv(Fp_n);
    Fe_tr     = Fdef * Fp_n_inv;
    be_tr     = Fe_tr * Fe_tr';
    Jdef_e_tr = det(Fe_tr);

    % find eigenvalues and eigenvectors of be_tr
    [V,D]       = eig(be_tr);
    lambda1e_tr = sqrt(D(1, 1));
    lambda2e_tr = sqrt(D(2, 2));
    lambda3e_tr = sqrt(D(3, 3));
    Jdef_e_tr   = lambda1e_tr * lambda2e_tr * lambda3e_tr;
    n1_tr       = V(:, 1);
    n2_tr       = V(:, 2);
    n3_tr       = V(:, 3);

    % calculate trial Kirchhoff stress
    tau1_tr    = lambda * log(Jdef_e_tr) - mu + mu * lambda1e_tr * lambda1e_tr;
    tau2_tr    = lambda * log(Jdef_e_tr) - mu + mu * lambda2e_tr * lambda2e_tr;
    tau3_tr    = lambda * log(Jdef_e_tr) - mu + mu * lambda3e_tr * lambda3e_tr;
    tau_tr     = tau1_tr*(n1_tr*n1_tr') + tau2_tr*(n2_tr*n2_tr') + tau3_tr*(n3_tr*n3_tr');
    p_tau_tr   = (tau1_tr + tau2_tr + tau3_tr) / 3;
    devtau1_tr = tau1_tr - p_tau_tr;
    devtau2_tr = tau2_tr - p_tau_tr;
    devtau3_tr = tau3_tr - p_tau_tr;

    % calculate trial elastic left stretch and rotation
    ve_tr = sqrtm(be_tr);
    Re_tr = inv(ve_tr) * Fe_tr;

    % calculate trial yield function
    c_tau          = c_tau_n;
    devtau_tr_norm = sqrt(devtau1_tr*devtau1_tr + devtau2_tr*devtau2_tr + devtau3_tr*devtau3_tr);
    f_tr           = devtau_tr_norm - Aphi*c_tau_n + Bphi*p_tau_tr;

    %dgdtau_tr
    dgdtau1_tr      = (devtau1_tr/devtau_tr_norm) + Bpsi/3;
    dgdtau2_tr      = (devtau2_tr/devtau_tr_norm) + Bpsi/3;
    dgdtau3_tr      = (devtau3_tr/devtau_tr_norm) + Bpsi/3;
    trace_dgdtau_tr = dgdtau1_tr + dgdtau2_tr + dgdtau3_tr;
    dgdtau_tr       = dgdtau1_tr*(n1_tr*n1_tr') + dgdtau2_tr*(n2_tr*n2_tr') + dgdtau3_tr*(n3_tr*n3_tr');

    %dfdtau_tr
    dfdtau1_tr      = (devtau1_tr/devtau_tr_norm) + Bphi/3;
    dfdtau2_tr      = (devtau2_tr/devtau_tr_norm) + Bphi/3;
    dfdtau3_tr      = (devtau3_tr/devtau_tr_norm) + Bphi/3;
    trace_dfdtau_tr = dfdtau1_tr + dfdtau2_tr + dfdtau3_tr;
    dfdtau_tr       = dfdtau1_tr*(n1_tr*n1_tr') + dfdtau2_tr*(n2_tr*n2_tr') + dfdtau3_tr*(n3_tr*n3_tr');

    % check if elastic, or elasto-plastic
    if (f_tr < 0)  % elastic
        Dg    = 0; % increment of plastic multiplier
        tau   = tau_tr;
        Fp    = Fp_n;
        c_tau = c_tau_n;
        be    = be_tr;
    else % plastic
        lambda1e = lambda1e_tr;
        lambda2e = lambda2e_tr;
        lambda3e = lambda3e_tr;
        Jdef_e   = lambda1e*lambda2e*lambda3e;

        % calculate Kirchhoff stress and deviatoric part
        tau1        = lambda*log(Jdef_e) - mu + mu*lambda1e*lambda1e;
        tau2        = lambda*log(Jdef_e) - mu + mu*lambda2e*lambda2e;
        tau3        = lambda*log(Jdef_e) - mu + mu*lambda3e*lambda3e;
        p_tau       = (tau1 + tau2 + tau3)/3;
        devtau1     = tau1 - p_tau;
        devtau2     = tau2 - p_tau;
        devtau3     = tau3 - p_tau;
        devtau_norm = sqrt(devtau1*devtau1 + devtau2*devtau2 + devtau3*devtau3);

        %solve for Dg (Delta_gamma)
        k                = 0;
        iter_break_local = 10;
        fyield           = f_tr;
        Dg_iter          = 0;
        while ( abs(fyield) > fyield_atol && abs(fyield/f_tr) > fyield_rtol )
            k = k+1;

            % derivative of devtau wrt Dg
            dtau1dDg      = -lambda*trace_dgdtau_tr - 2*mu*lambda1e*lambda1e*dgdtau1_tr;
            dtau2dDg      = -lambda*trace_dgdtau_tr - 2*mu*lambda2e*lambda2e*dgdtau2_tr;
            dtau3dDg      = -lambda*trace_dgdtau_tr - 2*mu*lambda3e*lambda3e*dgdtau3_tr;
            trace_dtaudDg = dtau1dDg + dtau2dDg + dtau3dDg;
            devdtau1dDg   = dtau1dDg - trace_dtaudDg/3;
            devdtau2dDg   = dtau2dDg - trace_dtaudDg/3;
            devdtau3dDg   = dtau3dDg - trace_dtaudDg/3;

            % derivative of devtau_norm wrt Dg
            ddevtau_norm_dDg = (devtau1*devdtau1dDg + devtau2*devdtau2dDg + devtau3*devdtau3dDg) / devtau_norm;

            % yield tangent
            dfdDg = ddevtau_norm_dDg - Hc*Aphi*Apsi + (Bphi/3)*trace_dtaudDg;

            % increment of Dg
            del_Dg  = -fyield / dfdDg;
            Dg_iter = Dg_iter + del_Dg;
            if Dg_iter < 0
                Dg_iter
            end

            %update elastic stretches
            lambda1e = lambda1e_tr * exp(-Dg_iter*dgdtau1_tr);
            lambda2e = lambda2e_tr * exp(-Dg_iter*dgdtau2_tr);
            lambda3e = lambda3e_tr * exp(-Dg_iter*dgdtau3_tr);
            Jdef_e   = lambda1e*lambda2e*lambda3e;

            % calculate Kirchhoff stress and deviatoric part
            tau1        = lambda*log(Jdef_e) - mu + mu*lambda1e*lambda1e;
            tau2        = lambda*log(Jdef_e) - mu + mu*lambda2e*lambda2e;
            tau3        = lambda*log(Jdef_e) - mu + mu*lambda3e*lambda3e;
            p_tau       = (tau1 + tau2 + tau3)/3;
            devtau1     = tau1 - p_tau;
            devtau2     = tau2 - p_tau;
            devtau3     = tau3 - p_tau;
            devtau_norm = sqrt(devtau1*devtau1 + devtau2*devtau2 + devtau3*devtau3);

            % update c_tau
            c_tau = c_tau_n + Dg_iter*Hc*Apsi;

            % update fyield
            fyield = devtau_norm - Aphi*c_tau + Bphi*p_tau;

            if k==iter_break_local
                k
                error('reached max number of local iterations for Dg')
            end
        end
        %update Fp
        Fp = Re_tr' * expm(Dg_iter * dgdtau_tr) * Re_tr * Fp_n;

        %update tau
        tau = tau1*(n1_tr*n1_tr') + tau2*(n2_tr*n2_tr') + tau3*(n3_tr*n3_tr');

        %update Dg
        Dg = Dg_iter;

        %update be
        Fe_ = Fdef*inv(Fp);
        be  = Fe_ * Fe_';
        %be = (lambda1e^2)*(n1_tr*n1_tr') + (lambda2e^2)*(n2_tr*n2_tr') + (lambda3e^2)*(n3_tr*n3_tr'); % this is correct too
    end

    % stresses
    Pstress = tau*Fdef_inv';
    Sstress = Fdef_inv*Pstress;
    sigma   = (1/Jdef)*tau;

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
    stress_el(:, :,  1)  = Sstress;
    stress_el(:, :,  2)  = Pstress;
    stress_el(:, :,  3)  = sigma;
    stress_el(:, :,  4)  = Estrain;
    stress_el(:, :,  5)  = estrain;
    stress_el(:, :,  6)  = eHstrain;
    stress_el(:, :,  7)  = tau;
    stress_el(:, :,  8)  = be;
    stress_el(:, :,  9)  = dfdtau_tr;
    stress_el(:, :, 10)  = dgdtau_tr;
    stress_el(:, :, 11)  = sigma_princ;
    stress_el(:, :, 12)  = sigma_rot;
    stress_el(:, :, 13)  = eHstrain_princ;

    % calculate stress invariants for output
    sig_mean = trace(sigma)/3;
    sig_dev  = sigma - sig_mean*eye_mat;
    sig_VM   = sqrt(3/2)*norm(sig_dev);

    % order of storage: p, q, c_tau, Dg, Jdef
    isv_el = [ sig_mean sig_VM c_tau Dg Jdef ];

    % return Fp
    Fp_el = Fp;
end
