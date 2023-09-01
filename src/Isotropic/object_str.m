%======================================================================%
%           Computational Micromechanics Laboratory (CMML)             %
%              Middle East Technical University (METU)                 %
%----------------------------------------------------------------------%
%                        Hyper-Data Toolbox - v1.0                     %
%----------------------------------------------------------------------%
%   Authors: Recep Durna, Alp Kağan Açan, Oğuz Ziya Tikenoğulları      %
%   Software License: GNU GPL v3.0   Copyright (C) 2023 Hüsnü Dal      %
%======================================================================%

% This function will calculate the error function SYMBOLICALLY and the 
% Symbolic expressions will be used in datadriven function to optimize the
% parameters
function [err_ut, err_et,err_ps,err_be, Pdata_ut, Pdata_et,Pdata_ps,Pdata1_be,Pdata2_be] = ...
   object_str(control_pts,degree,lam_ut,lam_et,lam_ps,lam1_be,lam2_be,P_ut,P_et,P_ps,P1_be,P2_be)

n_ut = length(lam_ut);
n_et = length(lam_et);
n_ps = length(lam_ps);
n_be = length(lam1_be);

% Deformation gradient, Left Cauchy Green tensor and Invariant calculation
[F_ut, F_et, F_ps] = F_calc(lam_ut,lam_et,lam_ps,n_ut,n_et,n_ps);
F_be = F_calc_be(lam1_be,lam2_be,n_be);
[b_ut, b_et, b_ps] = b_calc(F_ut,F_et,F_ps,n_ut,n_et,n_ps);
b_be = b_calc_be(F_be,n_be);

% Stretch calculations
lambda1_ut = zeros(n_ut,1);
lambda2_ut = zeros(n_ut,1);
lambda3_ut = zeros(n_ut,1);
 for m = 1:n_ut
      lambda1_ut(m,1) = lambda1_ut(m,1) + F_ut(1,1,m);
      lambda2_ut(m,1) = lambda2_ut(m,1) + F_ut(2,2,m);
      lambda3_ut(m,1) = lambda3_ut(m,1) + F_ut(3,3,m);
 end

lambda1_et = zeros(n_et,1);
lambda2_et = zeros(n_et,1);
lambda3_et = zeros(n_et,1);
 for m = 1:n_et
      lambda1_et(m,1) = lambda1_et(m,1) + F_et(1,1,m);
      lambda2_et(m,1) = lambda2_et(m,1) + F_et(2,2,m);
      lambda3_et(m,1) = lambda3_et(m,1) + F_et(3,3,m);
 end

lambda1_ps = zeros(n_ps,1);
lambda2_ps = zeros(n_ps,1);
lambda3_ps = zeros(n_ps,1);
 for m = 1:n_ps
      lambda1_ps(m,1) = lambda1_ps(m,1) + F_ps(1,1,m);
      lambda2_ps(m,1) = lambda2_ps(m,1) + F_ps(2,2,m);
      lambda3_ps(m,1) = lambda3_ps(m,1) + F_ps(3,3,m);
 end

lambda1_be = zeros(n_be,1);
lambda2_be = zeros(n_be,1);
lambda3_be = zeros(n_be,1);
 for m = 1:n_be
      lambda1_be(m,1) = lambda1_be(m,1) + F_be(1,1,m);
      lambda2_be(m,1) = lambda2_be(m,1) + F_be(2,2,m);
      lambda3_be(m,1) = lambda3_be(m,1) + F_be(3,3,m);
 end

lambda_ut = [lambda1_ut lambda2_ut lambda3_ut];
lambda_et = [lambda1_et lambda2_et lambda3_et];
lambda_ps = [lambda1_ps lambda2_ps lambda3_ps];
lambda_be = [lambda1_be lambda2_be lambda3_be];

lambda_max = max([lambda1_ut; lambda2_ut; lambda3_ut; ...
                  lambda1_et; lambda2_et; lambda3_et; ...
                  lambda1_ps; lambda2_ps; lambda3_ps; ...
                  lambda1_be; lambda2_be; lambda3_be]);
lambda_min = min([lambda1_ut; lambda2_ut; lambda3_ut; ...
                  lambda1_et; lambda2_et; lambda3_et; ...
                  lambda1_ps; lambda2_ps; lambda3_ps; ...
                  lambda1_be; lambda2_be; lambda3_be]);
% Inverse stretch calculations
lambda1_inv_ut = zeros(n_ut,1);                  
lambda2_inv_ut = zeros(n_ut,1);
lambda3_inv_ut = zeros(n_ut,1);
 for m = 1:n_ut
      lambda1_inv_ut(m,1) = lambda1_inv_ut(m,1) + 1/F_ut(1,1,m);
      lambda2_inv_ut(m,1) = lambda2_inv_ut(m,1) + 1/F_ut(2,2,m);
      lambda3_inv_ut(m,1) = lambda3_inv_ut(m,1) + 1/F_ut(3,3,m);
 end

lambda1_inv_et = zeros(n_et,1);                  
lambda2_inv_et = zeros(n_et,1);
lambda3_inv_et = zeros(n_et,1);
 for m = 1:n_et
      lambda1_inv_et(m,1) = lambda1_inv_et(m,1) + 1/F_et(1,1,m);
      lambda2_inv_et(m,1) = lambda2_inv_et(m,1) + 1/F_et(2,2,m);
      lambda3_inv_et(m,1) = lambda3_inv_et(m,1) + 1/F_et(3,3,m);
 end

lambda1_inv_ps = zeros(n_ps,1);                  
lambda2_inv_ps = zeros(n_ps,1);
lambda3_inv_ps = zeros(n_ps,1);
 for m = 1:n_ps
      lambda1_inv_ps(m,1) = lambda1_inv_ps(m,1) + 1/F_ps(1,1,m);
      lambda2_inv_ps(m,1) = lambda2_inv_ps(m,1) + 1/F_ps(2,2,m);
      lambda3_inv_ps(m,1) = lambda3_inv_ps(m,1) + 1/F_ps(3,3,m);
 end

lambda1_inv_be = zeros(n_be,1);                  
lambda2_inv_be = zeros(n_be,1);
lambda3_inv_be = zeros(n_be,1);
 for m = 1:n_be
      lambda1_inv_be(m,1) = lambda1_inv_be(m,1) + 1/F_be(1,1,m);
      lambda2_inv_be(m,1) = lambda2_inv_be(m,1) + 1/F_be(2,2,m);
      lambda3_inv_be(m,1) = lambda3_inv_be(m,1) + 1/F_be(3,3,m);
 end

lambda_inv_ut = [lambda1_inv_ut lambda2_inv_ut lambda3_inv_ut];
lambda_inv_et = [lambda1_inv_et lambda2_inv_et lambda3_inv_et];
lambda_inv_ps = [lambda1_inv_ps lambda2_inv_ps lambda3_inv_ps];
lambda_inv_be = [lambda1_inv_be lambda2_inv_be lambda3_inv_be];

invlambda_max = max([lambda1_inv_ut; lambda2_inv_ut; lambda3_inv_ut; ...
                     lambda1_inv_et; lambda2_inv_et; lambda3_inv_et; ...
                     lambda1_inv_ps; lambda2_inv_ps; lambda3_inv_ps; ...
                     lambda1_inv_be; lambda2_inv_be; lambda3_inv_be]);
invlambda_min = min([lambda1_inv_ut; lambda2_inv_ut; lambda3_inv_ut; ...
                     lambda1_inv_et; lambda2_inv_et; lambda3_inv_et; ...
                     lambda1_inv_ps; lambda2_inv_ps; lambda3_inv_ps; ...
                     lambda1_inv_be; lambda2_inv_be; lambda3_inv_be]);
% ---------------------------------------------%

% Generate bsplines
[p1_ut_1]=Bspline_generator(control_pts(:,1),[lambda_min,lambda_max],degree,lambda1_ut(:,1));
[p1_ut_2]=Bspline_generator(control_pts(:,1),[lambda_min,lambda_max],degree,lambda2_ut(:,1));
[p1_ut_3]=Bspline_generator(control_pts(:,1),[lambda_min,lambda_max],degree,lambda3_ut(:,1));
[p2_ut_1]=Bspline_generator(control_pts(:,2),[invlambda_min,invlambda_max],degree,lambda1_inv_ut(:,1));
[p2_ut_2]=Bspline_generator(control_pts(:,2),[invlambda_min,invlambda_max],degree,lambda2_inv_ut(:,1));
[p2_ut_3]=Bspline_generator(control_pts(:,2),[invlambda_min,invlambda_max],degree,lambda3_inv_ut(:,1));
p1_ut = [p1_ut_1 p1_ut_2 p1_ut_3];
p2_ut = [p2_ut_1 p2_ut_2 p2_ut_3];

[p1_et_1]=Bspline_generator(control_pts(:,1),[lambda_min,lambda_max],degree,lambda1_et(:,1));
[p1_et_2]=Bspline_generator(control_pts(:,1),[lambda_min,lambda_max],degree,lambda2_et(:,1));
[p1_et_3]=Bspline_generator(control_pts(:,1),[lambda_min,lambda_max],degree,lambda3_et(:,1));
[p2_et_1]=Bspline_generator(control_pts(:,2),[invlambda_min,invlambda_max],degree,lambda1_inv_et(:,1));
[p2_et_2]=Bspline_generator(control_pts(:,2),[invlambda_min,invlambda_max],degree,lambda2_inv_et(:,1));
[p2_et_3]=Bspline_generator(control_pts(:,2),[invlambda_min,invlambda_max],degree,lambda3_inv_et(:,1));
p1_et = [p1_et_1 p1_et_2 p1_et_3];
p2_et = [p2_et_1 p2_et_2 p2_et_3];

[p1_ps_1]=Bspline_generator(control_pts(:,1),[lambda_min,lambda_max],degree,lambda1_ps(:,1));
[p1_ps_2]=Bspline_generator(control_pts(:,1),[lambda_min,lambda_max],degree,lambda2_ps(:,1));
[p1_ps_3]=Bspline_generator(control_pts(:,1),[lambda_min,lambda_max],degree,lambda3_ps(:,1));
[p2_ps_1]=Bspline_generator(control_pts(:,2),[invlambda_min,invlambda_max],degree,lambda1_inv_ps(:,1));
[p2_ps_2]=Bspline_generator(control_pts(:,2),[invlambda_min,invlambda_max],degree,lambda2_inv_ps(:,1));
[p2_ps_3]=Bspline_generator(control_pts(:,2),[invlambda_min,invlambda_max],degree,lambda3_inv_ps(:,1));
p1_ps = [p1_ps_1 p1_ps_2 p1_ps_3];
p2_ps = [p2_ps_1 p2_ps_2 p2_ps_3];

[p1_be_1]=Bspline_generator(control_pts(:,1),[lambda_min,lambda_max],degree,lambda1_be(:,1));
[p1_be_2]=Bspline_generator(control_pts(:,1),[lambda_min,lambda_max],degree,lambda2_be(:,1));
[p1_be_3]=Bspline_generator(control_pts(:,1),[lambda_min,lambda_max],degree,lambda3_be(:,1));
[p2_be_1]=Bspline_generator(control_pts(:,2),[invlambda_min,invlambda_max],degree,lambda1_inv_be(:,1));
[p2_be_2]=Bspline_generator(control_pts(:,2),[invlambda_min,invlambda_max],degree,lambda2_inv_be(:,1));
[p2_be_3]=Bspline_generator(control_pts(:,2),[invlambda_min,invlambda_max],degree,lambda3_inv_be(:,1));
p1_be = [p1_be_1 p1_be_2 p1_be_3];
p2_be = [p2_be_1 p2_be_2 p2_be_3];

% Run the constitutive model
if n_ut > 0
    [Pdata_ut,~] = data_driven3_function(F_ut,lambda_ut,lambda_inv_ut,b_ut,n_ut,p1_ut,p2_ut);
else
    Pdata_ut = [];
end
if n_et > 0
    [Pdata_et,~] = data_driven3_function(F_et,lambda_et,lambda_inv_et,b_et,n_et,p1_et,p2_et);
else
    Pdata_et = [];
end
if n_ps > 0
    [Pdata_ps,~] = data_driven3_function(F_ps,lambda_ps,lambda_inv_ps,b_ps,n_ps,p1_ps,p2_ps);
else
    Pdata_ps = [];
end
if n_be>0
    [Pdata1_be,Pdata2_be] = data_driven2_function(F_be,lambda_be,lambda_inv_be,b_be,n_be,p1_be,p2_be);
else
    Pdata1_be = [];
    Pdata2_be = [];
end

% Calculate mean square errors
if n_ut > 0
    err_ut = 1/n_ut * sum((Pdata_ut-P_ut).^2);
else
    err_ut = 0;
end
if n_et > 0
    err_et = 1/n_et * sum((Pdata_et-P_et).^2);
else
    err_et = 0;
end
if n_ps > 0
    err_ps = 1/n_ps * sum((Pdata_ps-P_ps).^2);
else
    err_ps = 0;
end
if n_be>0
    err_be = 1/n_be * sum((Pdata2_be-P2_be).^2) + 1/n_be * sum((Pdata1_be-P1_be).^2);
else
    err_be = 0;
end

end
