# Mock the exponential return mapping algorithm from CVEN7511

function  element_stress_isv(coordsx, d, params, isv_el_n, Fp_el_n, theta)
    # nodal coordinates of elements in local node numbering
    x1=coordsx[1];  y1=coordsx[2];  z1=coordsx[3]; 
    x2=coordsx[4];  y2=coordsx[5];  z2=coordsx[6];
    x3=coordsx[7];  y3=coordsx[8];  z3=coordsx[9];
    x4=coordsx[10]; y4=coordsx[11]; z4=coordsx[12];
    x5=coordsx[13]; y5=coordsx[14]; z5=coordsx[15];
    x6=coordsx[16]; y6=coordsx[17]; z6=coordsx[18];
    x7=coordsx[19]; y7=coordsx[20]; z7=coordsx[21];
    x8=coordsx[22]; y8=coordsx[23]; z8=coordsx[24];

    #params = [ lambda mu Aphi Bphi Apsi Bpsi Hc numips nstress nisv ndim ];
    # recall parameters
    lambda=params[1];
    mu=params[2];
    Aphi=params[3];
    Bphi=params[4];
    Apsi=params[5];
    Bpsi=params[6];
    Hc=params[7];

    # tolerances for local Newton-Raphson algorithm to solve for Delta gamma
    fyield_atol = 1e-6;
    fyield_rtol = 1e-8;

    numips=round(params[8]);
    nstress=round(params[9]);
    nisv=round(params[10]);
    ndim=round(params[11]);

    # initialize size of stress, and ISVs at Gauss pts, to be returned
    stress_el=zeros(ndim,ndim,nstress,numips);
    isv_el=zeros(nisv,numips);
    Fp_el=zeros(ndim,ndim,numips);

    # set Gauss point coordinates in xi,eta,zeta space
    constant = 1/sqrt(3);
    xi_vect[1,:] = [-constant -constant -constant];
    xi_vect[2,:] = [constant -constant -constant];
    xi_vect[3,:] = [constant constant -constant];
    xi_vect[4,:] = [-constant constant -constant];
    xi_vect[5,:] = [-constant -constant constant];
    xi_vect[6,:] = [constant -constant constant];
    xi_vect[7,:] = [constant constant constant];
    xi_vect[8,:] = [-constant constant constant];

    # loop through the Gauss points
    for ip=1:numips
        # integration points in natural coordinates
        xi=xi_vect(ip,1);
        eta=xi_vect(ip,2);
        zeta=xi_vect(ip,3);
    
        # derivatives of shape functions with respect to xi
        dN1_dxi = -0.125*(1-eta)*(1-zeta);
        dN2_dxi = -dN1_dxi;
        dN3_dxi = 0.125*(1+eta)*(1-zeta);
        dN4_dxi = -dN3_dxi;
        dN5_dxi = -0.125*(1-eta)*(1+zeta);
        dN6_dxi = -dN5_dxi;
        dN7_dxi = 0.125*(1+eta)*(1+zeta);
        dN8_dxi = -dN7_dxi;

        # derivatives of shape functions with respect to eta
        dN1_deta = -0.125*(1-xi)*(1-zeta);
        dN2_deta = -0.125*(1+xi)*(1-zeta);
        dN3_deta = -dN2_deta;
        dN4_deta = -dN1_deta;
        dN5_deta = -0.125*(1-xi)*(1+zeta);
        dN6_deta = -0.125*(1+xi)*(1+zeta);
        dN7_deta = -dN6_deta;
        dN8_deta = -dN5_deta;

        # derivatives of shape functions with respect to zeta
        dN1_dzeta = -0.125*(1-xi)*(1-eta);
        dN2_dzeta = -0.125*(1+xi)*(1-eta);
        dN3_dzeta = -0.125*(1+xi)*(1+eta);
        dN4_dzeta = -0.125*(1-xi)*(1+eta);
        dN5_dzeta = -dN1_dzeta;
        dN6_dzeta = -dN2_dzeta;
        dN7_dzeta = -dN3_dzeta;
        dN8_dzeta = -dN4_dzeta;

        # calculate jacobian, its determinant, and its inverse
        dx_dxi = dN1_dxi*x1 + dN2_dxi*x2 + dN3_dxi*x3 + dN4_dxi*x4 + dN5_dxi*x5 + dN6_dxi*x6 + dN7_dxi*x7 + dN8_dxi*x8;
        dx_deta = dN1_deta*x1 + dN2_deta*x2 + dN3_deta*x3 + dN4_deta*x4 + dN5_deta*x5 + dN6_deta*x6 + dN7_deta*x7 + dN8_deta*x8;
        dx_dzeta = dN1_dzeta*x1 + dN2_dzeta*x2 + dN3_dzeta*x3 + dN4_dzeta*x4 + dN5_dzeta*x5 + dN6_dzeta*x6 + dN7_dzeta*x7 + dN8_dzeta*x8;
        dy_dxi = dN1_dxi*y1 + dN2_dxi*y2 + dN3_dxi*y3 + dN4_dxi*y4 + dN5_dxi*y5 + dN6_dxi*y6 + dN7_dxi*y7 + dN8_dxi*y8;
        dy_deta = dN1_deta*y1 + dN2_deta*y2 + dN3_deta*y3 + dN4_deta*y4 + dN5_deta*y5 + dN6_deta*y6 + dN7_deta*y7 + dN8_deta*y8;
        dy_dzeta = dN1_dzeta*y1 + dN2_dzeta*y2 + dN3_dzeta*y3 + dN4_dzeta*y4 + dN5_dzeta*y5 + dN6_dzeta*y6 + dN7_dzeta*y7 + dN8_dzeta*y8;
        dz_dxi = dN1_dxi*z1 + dN2_dxi*z2 + dN3_dxi*z3 + dN4_dxi*z4 + dN5_dxi*z5 + dN6_dxi*z6 + dN7_dxi*z7 + dN8_dxi*z8;
        dz_deta = dN1_deta*z1 + dN2_deta*z2 + dN3_deta*z3 + dN4_deta*z4 + dN5_deta*z5 + dN6_deta*z6 + dN7_deta*z7 + dN8_deta*z8;
        dz_dzeta = dN1_dzeta*z1 + dN2_dzeta*z2 + dN3_dzeta*z3 + dN4_dzeta*z4 + dN5_dzeta*z5 + dN6_dzeta*z6 + dN7_dzeta*z7 + dN8_dzeta*z8;
        Je = [dx_dxi dx_deta dx_dzeta ; dy_dxi dy_deta dy_dzeta ; dz_dxi dz_deta dz_dzeta];
        jdet=det(Je);
        Jeinv=inv(Je);

        #shape function derivatives with respect to X,Y,Z
        dN1_dx_vect = [dN1_dxi dN1_deta dN1_dzeta]*Jeinv;
        dN2_dx_vect = [dN2_dxi dN2_deta dN2_dzeta]*Jeinv;
        dN3_dx_vect = [dN3_dxi dN3_deta dN3_dzeta]*Jeinv;
        dN4_dx_vect = [dN4_dxi dN4_deta dN4_dzeta]*Jeinv;    
        dN5_dx_vect = [dN5_dxi dN5_deta dN5_dzeta]*Jeinv;
        dN6_dx_vect = [dN6_dxi dN6_deta dN6_dzeta]*Jeinv;
        dN7_dx_vect = [dN7_dxi dN7_deta dN7_dzeta]*Jeinv;
        dN8_dx_vect = [dN8_dxi dN8_deta dN8_dzeta]*Jeinv;

        # strain-displacement matrix
        Bu=[ dN1_dx_vect[1] 0 0 dN2_dx_vect[1] 0 0 dN3_dx_vect[1] 0 0 dN4_dx_vect[1] 0 0 dN5_dx_vect[1] 0 0 dN6_dx_vect[1] 0 0 dN7_dx_vect[1] 0 0 dN8_dx_vect[1] 0 0;
             dN1_dx_vect[2] 0 0 dN2_dx_vect[2] 0 0 dN3_dx_vect[2] 0 0 dN4_dx_vect[2] 0 0 dN5_dx_vect[2] 0 0 dN6_dx_vect[2] 0 0 dN7_dx_vect[2] 0 0 dN8_dx_vect[2] 0 0;
             dN1_dx_vect[3] 0 0 dN2_dx_vect[3] 0 0 dN3_dx_vect[3] 0 0 dN4_dx_vect[3] 0 0 dN5_dx_vect[3] 0 0 dN6_dx_vect[3] 0 0 dN7_dx_vect[3] 0 0 dN8_dx_vect[3] 0 0;
             0 dN1_dx_vect[1] 0 0 dN2_dx_vect[1] 0 0 dN3_dx_vect[1] 0 0 dN4_dx_vect[1] 0 0 dN5_dx_vect[1] 0 0 dN6_dx_vect[1] 0 0 dN7_dx_vect[1] 0 0 dN8_dx_vect[1] 0;
             0 dN1_dx_vect[2] 0 0 dN2_dx_vect[2] 0 0 dN3_dx_vect[2] 0 0 dN4_dx_vect[2] 0 0 dN5_dx_vect[2] 0 0 dN6_dx_vect[2] 0 0 dN7_dx_vect[2] 0 0 dN8_dx_vect[2] 0;
             0 dN1_dx_vect[3] 0 0 dN2_dx_vect[3] 0 0 dN3_dx_vect[3] 0 0 dN4_dx_vect[3] 0 0 dN5_dx_vect[3] 0 0 dN6_dx_vect[3] 0 0 dN7_dx_vect[3] 0 0 dN8_dx_vect[3] 0;
             0 0 dN1_dx_vect[1] 0 0 dN2_dx_vect[1] 0 0 dN3_dx_vect[1] 0 0 dN4_dx_vect[1] 0 0 dN5_dx_vect[1] 0 0 dN6_dx_vect[1] 0 0 dN7_dx_vect[1] 0 0 dN8_dx_vect[1];
             0 0 dN1_dx_vect[2] 0 0 dN2_dx_vect[2] 0 0 dN3_dx_vect[2] 0 0 dN4_dx_vect[2] 0 0 dN5_dx_vect[2] 0 0 dN6_dx_vect[2] 0 0 dN7_dx_vect[2] 0 0 dN8_dx_vect[2];
             0 0 dN1_dx_vect[3] 0 0 dN2_dx_vect[3] 0 0 dN3_dx_vect[3] 0 0 dN4_dx_vect[3] 0 0 dN5_dx_vect[3] 0 0 dN6_dx_vect[3] 0 0 dN7_dx_vect[3] 0 0 dN8_dx_vect[3]
             ];

        # total deformation tensors
        dudX = Bu*d';
        Fdef = eye(3) + [ dudX[1] dudX[2] dudX[3]; dudX[4] dudX[5] dudX[6]; dudX[7] dudX[8] dudX[9] ];
        Jdef = det(Fdef);
        Fdef_inv = inv(Fdef);
        Cdef = Fdef' * Fdef;
        Cdef_inv = inv(Cdef);
        bdef = Fdef * Fdef';
        bdef_inv = inv(bdef);
        vdef = sqrtm(bdef);

        # 3x3 identity matrix
        eye_mat = eye(3);

        # total strains
        Estrain = (Cdef-eye_mat)/2;
        estrain = (eye_mat-inv(bdef))/2;
        eHstrain = logm(vdef);

        # total principal Hencky strain
        lambda_eig = eig(eHstrain);
        eHstrain_princ = [max(lambda_eig) 0         0 ;
                          0           lambda_eig[2] 0 ;
                          0           0         min(lambda_eig) ];

        # calculate trial elastic deformation gradient and left elastic Cauchy Green tensor
        Fp_n = Fp_el_n(:,:,ip);
        Fp_n_inv = inv(Fp_n);
        Fe_tr = Fdef * Fp_n_inv;
        be_tr = Fe_tr * Fe_tr';
        Jdef_e_tr = det(Fe_tr);

        # find eigenvalues and eigenvectors of be_tr
        V,D = eig(be_tr);
        lambda1e_tr = sqrt(D(1,1));
        lambda2e_tr = sqrt(D(2,2));
        lambda3e_tr = sqrt(D(3,3));
        Jdef_e_tr = lambda1e_tr*lambda2e_tr*lambda3e_tr;
        n1_tr = V(:,1);
        n2_tr = V(:,2);
        n3_tr = V(:,3);

        # calculate trial Kirchhoff stress
        tau1_tr = lambda*log(Jdef_e_tr) - mu + mu*lambda1e_tr*lambda1e_tr;
        tau2_tr = lambda*log(Jdef_e_tr) - mu + mu*lambda2e_tr*lambda2e_tr;
        tau3_tr = lambda*log(Jdef_e_tr) - mu + mu*lambda3e_tr*lambda3e_tr;
        tau_tr = tau1_tr*(n1_tr*n1_tr') + tau2_tr*(n2_tr*n2_tr') + tau3_tr*(n3_tr*n3_tr');
        p_tau_tr = (tau1_tr + tau2_tr + tau3_tr)/3;
        devtau1_tr = tau1_tr - p_tau_tr;
        devtau2_tr = tau2_tr - p_tau_tr;
        devtau3_tr = tau3_tr - p_tau_tr;
        # or (comment out, for debugging)
        #tau_tr_alt = (lambda*log(Je_tr)-mu)*eye_mat + mu*be_tr;

        # calculate trial elastic left stretch and rotation
        ve_tr = sqrtm(be_tr);
        Re_tr = inv(ve_tr) * Fe_tr;

        # calculate trial yield function
        c_tau_n = isv_el_n(3,ip);
        c_tau = c_tau_n;
        devtau_tr_norm = sqrt(devtau1_tr*devtau1_tr+devtau2_tr*devtau2_tr+devtau3_tr*devtau3_tr);
        f_tr = devtau_tr_norm - Aphi*c_tau_n + Bphi*p_tau_tr;
        #f_tr = -1; #to debug

        #dgdtau_tr
        dgdtau1_tr = (devtau1_tr/devtau_tr_norm) + Bpsi/3;
        dgdtau2_tr = (devtau2_tr/devtau_tr_norm) + Bpsi/3;
        dgdtau3_tr = (devtau3_tr/devtau_tr_norm) + Bpsi/3;
        trace_dgdtau_tr = dgdtau1_tr+dgdtau2_tr+dgdtau3_tr;
        dgdtau_tr = dgdtau1_tr*(n1_tr*n1_tr') + dgdtau2_tr*(n2_tr*n2_tr') + dgdtau3_tr*(n3_tr*n3_tr');

        #dfdtau_tr
        dfdtau1_tr = (devtau1_tr/devtau_tr_norm) + Bphi/3;
        dfdtau2_tr = (devtau2_tr/devtau_tr_norm) + Bphi/3;
        dfdtau3_tr = (devtau3_tr/devtau_tr_norm) + Bphi/3;
        trace_dfdtau_tr = dfdtau1_tr+dfdtau2_tr+dfdtau3_tr;
        dfdtau_tr = dfdtau1_tr*(n1_tr*n1_tr') + dfdtau2_tr*(n2_tr*n2_tr') + dfdtau3_tr*(n3_tr*n3_tr');

        # check if elastic, or elasto-plastic
        if (f_tr < 0) #elastic
            Dg = 0; #increment of plastic multiplier
            tau = tau_tr;
            Fp = Fp_n;
            c_tau = c_tau_n;
            be = be_tr;
        else #plastic
            #Dg = 0; #to debug
            #
            lambda1e = lambda1e_tr;
            lambda2e = lambda2e_tr;
            lambda3e = lambda3e_tr;
            Jdef_e = lambda1e*lambda2e*lambda3e;
            #
            # calculate Kirchhoff stress and deviatoric part
            tau1 = lambda*log(Jdef_e) - mu + mu*lambda1e*lambda1e;
            tau2 = lambda*log(Jdef_e) - mu + mu*lambda2e*lambda2e;
            tau3 = lambda*log(Jdef_e) - mu + mu*lambda3e*lambda3e;
            p_tau = (tau1 + tau2 + tau3)/3;
            devtau1 = tau1 - p_tau;
            devtau2 = tau2 - p_tau;
            devtau3 = tau3 - p_tau;
            devtau_norm = sqrt(devtau1*devtau1+devtau2*devtau2+devtau3*devtau3);
            #
            #solve for Dg
            k=0;
            iter_break_local=10;
            fyield = f_tr;
            Dg_iter = 0;
            while ( abs(fyield) > fyield_atol && abs(fyield/f_tr) > fyield_rtol )
                k=k+1;
                #
                # derivative of devtau wrt Dg
                dtau1dDg = -lambda*trace_dgdtau_tr - 2*mu*lambda1e*lambda1e*dgdtau1_tr;
                dtau2dDg = -lambda*trace_dgdtau_tr - 2*mu*lambda2e*lambda2e*dgdtau2_tr;
                dtau3dDg = -lambda*trace_dgdtau_tr - 2*mu*lambda3e*lambda3e*dgdtau3_tr;
                trace_dtaudDg = dtau1dDg+dtau2dDg+dtau3dDg;
                devdtau1dDg = dtau1dDg - trace_dtaudDg/3;
                devdtau2dDg = dtau2dDg - trace_dtaudDg/3;
                devdtau3dDg = dtau3dDg - trace_dtaudDg/3;
                #
                # derivative of devtau_norm wrt Dg
                ddevtau_norm_dDg = (devtau1*devdtau1dDg+devtau2*devdtau2dDg+devtau3*devdtau3dDg)/devtau_norm;
                #
                # yield tangent
                dfdDg = ddevtau_norm_dDg - Hc*Aphi*Apsi + (Bphi/3)*trace_dtaudDg;
                #
                # increment of Dg
                del_Dg = -fyield/dfdDg;
                Dg_iter = Dg_iter + del_Dg;
                if Dg_iter < 0
                    Dg_iter
                end
                #
                #update elastic stretches
                lambda1e = lambda1e_tr*exp(-Dg_iter*dgdtau1_tr);
                lambda2e = lambda2e_tr*exp(-Dg_iter*dgdtau2_tr);
                lambda3e = lambda3e_tr*exp(-Dg_iter*dgdtau3_tr);
                Jdef_e = lambda1e*lambda2e*lambda3e;
                #
                # calculate Kirchhoff stress and deviatoric part
                tau1 = lambda*log(Jdef_e) - mu + mu*lambda1e*lambda1e;
                tau2 = lambda*log(Jdef_e) - mu + mu*lambda2e*lambda2e;
                tau3 = lambda*log(Jdef_e) - mu + mu*lambda3e*lambda3e;
                p_tau = (tau1 + tau2 + tau3)/3;
                devtau1 = tau1 - p_tau;
                devtau2 = tau2 - p_tau;
                devtau3 = tau3 - p_tau;
                devtau_norm = sqrt(devtau1*devtau1+devtau2*devtau2+devtau3*devtau3);
                #
                # update c_tau
                c_tau = c_tau_n + Dg_iter*Hc*Apsi;
                #
                # update fyield
                fyield = devtau_norm - Aphi*c_tau + Bphi*p_tau;
                #
                if k==iter_break_local
                    k
                    println("reached max number of local iterations for Dg")
                    exit()
                end
            end
            #update Fp
            Fp = Re_tr' * expm(Dg_iter * dgdtau_tr) * Re_tr * Fp_n;
            #update tau
            tau = tau1*(n1_tr*n1_tr') + tau2*(n2_tr*n2_tr') + tau3*(n3_tr*n3_tr');
            #update Dg
            Dg = Dg_iter;
            #update be
            Fe_ = Fdef*inv(Fp);
            be = Fe_ * Fe_';
            #be = (lambda1e^2)*(n1_tr*n1_tr') + (lambda2e^2)*(n2_tr*n2_tr') + (lambda3e^2)*(n3_tr*n3_tr'); # this is correct too
        end

        # stresses
        Pstress = tau*Fdef_inv';
        Sstress = Fdef_inv*Pstress;
        sigma = (1/Jdef)*tau;

        # rotate stress for non-rotating case
        Qrot = [cos(theta) sin(theta) 0;
               -sin(theta) cos(theta) 0;
                0          0          1];
        sigma_rot = Qrot * sigma * Qrot';    

        # principal Cauchy stress
        lambda_eig = eig(sigma);
        sigma_princ = [max(lambda_eig) 0         0 ;
                       0           lambda_eig[2] 0 ;
                       0           0         min(lambda_eig) ];

        # order of storage: S,P,sig,E,e,le,tau,be,dfdtau,dgdtau,sig_princ,sig_rot,le_princ
        stress_el[:,:,1,ip] = Sstress;
        stress_el[:,:,2,ip] = Pstress;
        stress_el[:,:,3,ip] = sigma;
        stress_el[:,:,4,ip] = Estrain;
        stress_el[:,:,5,ip] = estrain;
        stress_el[:,:,6,ip] = eHstrain;
        stress_el[:,:,7,ip] = tau;
        stress_el[:,:,8,ip] = be;
        stress_el[:,:,9,ip] = dfdtau_tr;
        stress_el[:,:,10,ip] = dgdtau_tr;
        stress_el[:,:,11,ip] = sigma_princ;
        stress_el[:,:,12,ip] = sigma_rot;
        stress_el[:,:,13,ip] = eHstrain_princ;

        # calculate stress invariants for output
        sig_mean = trace(sigma)/3;
        sig_dev = sigma - sig_mean*eye_mat;
        sig_VM = sqrt(3/2)*norm(sig_dev);

        # order of storage: p, q, c_tau, Dg, Jdef
        isv_el[1,ip] = sig_mean;
        isv_el[2,ip] = sig_VM;
        isv_el[3,ip] = c_tau;
        isv_el[4,ip] = Dg;    # Δγ
        isv_el[5,ip] = Jdef;

        # return Fp
        Fp_el[:,:,ip] = Fp;
    end

    return stress_el, isv_el, Fp_el
end 
