% This function is the main switching function for Extended Eight-Chain model 
function [P11_ut, P11_et ,P11_ps,w_par, x_par, P11_be1, P22_be1, err_ut, err_et, err_ps, err_be] = ...
abm(n_ut,n_et,n_ps,lam_ut,lam_et,lam_ps,P_ut,P_et,P_ps,wt,w,experiment,lam1_be,lam2_be,P11_be,P22_be,n_be)

c = sym('c',[1 6]);

% Deformation gradient, Left Cauchy Green tensor and Invariant calculation
[F_ut, F_et, F_ps] = F_calc(lam_ut,lam_et,lam_ps,n_ut,n_et,n_ps);
F_be = F_calc_be(lam1_be,lam2_be,n_be);
[b_ut, b_et, b_ps] = b_calc(F_ut,F_et,F_ps,n_ut,n_et,n_ps);
b_be = b_calc_be(F_be,n_be);
[I1_ut, I1_et, I1_ps] = I1_calc(b_ut,b_et,b_ps,n_ut,n_et,n_ps);
[I2_ut, I2_et, I2_ps] = I2_calc(b_ut,b_et,b_ps,n_ut,n_et,n_ps);
I1_be = I1_calc_be(b_be, n_be);
I2_be = I2_calc_be(b_be, n_be);
% This part calculates errors for uniaxial, equaibiaxial and pure shear
% deformations SYMBOLICALLY using function named arruda_boyce_function
if experiment(1) == 1
    [err_ut,~] = abm_function(F_ut,b_ut,n_ut,P_ut,I1_ut,I2_ut,c);
else
    err_ut = 0;
end
if experiment(2) == 1
    [err_et,~] = abm_function(F_et,b_et,n_et,P_et,I1_et,I2_et,c);
else 
    err_et = 0;
end
if experiment(3) == 1
    [err_ps,~] = abm_function(F_ps,b_ps,n_ps,P_ps,I1_ps,I2_ps,c);
else
    err_ps = 0;
end
if experiment(4) == 1
    [err_be1,~] = abm_function(F_be,b_be,n_be,P11_be,I1_be,I2_be,c);
    [~,err_be2] = abm_function(F_be,b_be,n_be,P22_be,I1_be,I2_be,c);
    err_be = err_be1 + err_be2;
else
    err_be = 0;
end
% Switching btw automatic and manual optimization of parameters

switch wt
    case 1
     % Automatic optimization function which takes weight factors as
     % parameters
        [x] = abm_auto(err_ut,err_et,err_ps,err_be,experiment);
        
        mu   = x(1);
        N    = x(2);
        mu2  = x(3);
        w_ut = x(4);
        w_et = x(5);
        w_ps = x(6);

        x_par = [mu N mu2];
        w_par = [w_ut w_et w_ps];

    case 2
        w_ut = w(1);
        w_et = w(2);
        w_ps = w(3);
        w_be = w(4);
        [x]  = abm_manual(err_ut,err_et,err_ps,err_be,w_ut,w_et,w_ps,w_be,experiment);
        mu   = x(1);
        N    = x(2);
        mu2  = x(3);
        
        x_par = [mu N mu2];
        w_par = [w_ut w_et w_ps];
end

% This part calculates errors for uniaxial, equaibiaxial and pure shear
% deformations NUMERICALLY using function named arruda_boyce_function_
[P11_ut, taubar_ut, tau_ut, err_ut,~,~] = abm_function_(F_ut,b_ut,n_ut,P_ut,I1_ut,I2_ut,x_par);
[P11_et, taubar_et, tau_et, err_et,~,~] = abm_function_(F_et,b_et,n_et,P_et,I1_et,I2_et,x_par);
[P11_ps, taubar_ps, tau_ps, err_ps,~,~] = abm_function_(F_ps,b_ps,n_ps,P_ps,I1_ps,I2_ps,x_par);
[~,taubar_be,tau_be,~,err_be2,P22_be1]   = abm_function_(F_be,b_be,n_be,P22_be,I1_be,I2_be,x_par);
[P11_be1,taubar_be,tau_be,err_be1,~,~]   = abm_function_(F_be,b_be,n_be,P11_be,I1_be,I2_be,x_par);
