% This function is the main switching function for ogden model 

function [P11_ut, P11_et ,P11_ps,w_par, x_par, P11_be1, P22_be1, err_ut, err_et, err_ps, err_be] = ...
ogden(n_ut,n_et,n_ps,lam_ut,lam_et,lam_ps,P_ut,P_et,P_ps,wt,w,experiment,lam1_be,lam2_be,P11_be,P22_be,n_be)

% Parameters including weight factors 
c = sym('c',[1 9]);

% Deformation gradient tensor calculation
[F_ut, F_et, F_ps] = F_calc(lam_ut,lam_et,lam_ps,n_ut,n_et,n_ps);
F_be = F_calc_be(lam1_be,lam2_be,n_be);

% This part calculates errors for uniaxial, equaibiaxial and pure shear
% deformations SYMBOLICALLY using function named ogden_function
if experiment(1) ~= 0
    [err_ut,~] = ogden_function(F_ut,n_ut,P_ut,c);
else
    err_ut = 0;
end
if experiment(2) ~= 0
    [err_et,~] = ogden_function(F_et,n_et,P_et,c);
else
    err_et = 0;
end
if experiment(3) ~= 0
    [err_ps,~] = ogden_function(F_ps,n_ps,P_ps,c);
else
    err_ps = 0;
end
if experiment(4) ~= 0
    [err_be1,~] = ogden_function(F_be,n_be,P11_be,c);
    [~,err_be2] = ogden_function(F_be,n_be,P22_be,c);
    err_be = err_be1 + err_be2;
else
    err_be = 0;
end

% Switching btw automatic and manual optimization of parameters
switch wt
    case 1
     % Automatic optimization function which takes weight factors as
     % parameters
     [x] = ogden_auto(err_ut,err_et,err_ps,err_be,experiment);

     a1  = x(1);
     a2  = x(2);
     a3  = x(3);
     mu1 = x(4);
     mu2 = x(5);
     mu3 = x(6);        
     w_ut = x(7);
     w_et = x(8);
     w_ps = x(9);

     x_par = [a1 a2 a3 mu1 mu2 mu3];  % Parameters of ogden model
     w_par = [w_ut w_et w_ps];        % Weight factors

    case 2
     % Manual optimization function which takes weight factors as
     % inputs     
     w_ut = w(1);
     w_et = w(2);
     w_ps = w(3);
     w_be = w(4);
     [x] = ogden_manual(err_ut,err_et,err_ps,err_be,w_ut,w_et,w_ps,w_be,experiment);

     a1  = x(1);
     a2  = x(2);
     a3  = x(3);
     mu1 = x(4);
     mu2 = x(5);
     mu3 = x(6);        

     x_par = [a1 a2 a3 mu1 mu2 mu3];  % Parameters of ogden model
     w_par = [w_ut w_et w_ps];        % Weight factors
end

% This part calculates errors for uniaxial, equaibiaxial and pure shear
% deformations NUMERICALLY using function named ogden_function_
[P11_ut, taubar_ut, tau_ut, err_ut,~,~] = ogden_function_(F_ut,n_ut,P_ut,x_par);
[P11_et, taubar_et, tau_et, err_et,~,~] = ogden_function_(F_et,n_et,P_et,x_par);
[P11_ps, taubar_ps, tau_ps, err_ps,~,~] = ogden_function_(F_ps,n_ps,P_ps,x_par);
[P11_be1, taubar_be, tau_be, err_be1, ~, ~] = ogden_function_(F_be,n_be,P11_be,x_par);
[~, taubar_be, tau_be, ~, err_be2, P22_be1] = ogden_function_(F_be,n_be,P22_be,x_par);
err_be = err_be1+err_be2;
