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
   object(control_pts,degree,lam_ut,lam_et,lam_ps,lam1_be,lam2_be,P_ut,P_et,P_ps,P1_be,P2_be)

n_ut = length(lam_ut);
n_et = length(lam_et);
n_ps = length(lam_ps);
n_be = length(lam1_be);

% Deformation gradient, Left Cauchy Green tensor and Invariant calculation
[F_ut, F_et, F_ps] = F_calc(lam_ut,lam_et,lam_ps,n_ut,n_et,n_ps);
F_be = F_calc_be(lam1_be,lam2_be,n_be);
[b_ut, b_et, b_ps] = b_calc(F_ut,F_et,F_ps,n_ut,n_et,n_ps);
b_be = b_calc_be(F_be,n_be);
[I1_ut, I1_et, I1_ps] = I1_calc(b_ut,b_et,b_ps,n_ut,n_et,n_ps);
[I2_ut, I2_et, I2_ps] = I2_calc(b_ut,b_et,b_ps,n_ut,n_et,n_ps);
I1_be = I1_calc_be(b_be,n_be);
I2_be = I2_calc_be(b_be,n_be);

I1_max= max([I1_ut;I1_et;I1_ps;I1_be]);
I1_min= min([I1_ut;I1_et;I1_ps;I1_be]);

I2_max= max([I2_ut;I2_et;I2_ps;I2_be]);
I2_min= min([I2_ut;I2_et;I2_ps;I2_be]);

% Generate bsplines
[p1_ut]=Bspline_generator(control_pts(:,1),[I1_min, I1_max],degree,I1_ut);
[p2_ut]=Bspline_generator(control_pts(:,2),[I2_min, I2_max],degree,I2_ut);

[p1_et]=Bspline_generator(control_pts(:,1),[I1_min, I1_max],degree,I1_et);
[p2_et]=Bspline_generator(control_pts(:,2),[I2_min, I2_max],degree,I2_et);

[p1_ps]=Bspline_generator(control_pts(:,1),[I1_min, I1_max],degree,I1_ps);
[p2_ps]=Bspline_generator(control_pts(:,2),[I2_min, I2_max],degree,I2_ps);

[p1_be]=Bspline_generator(control_pts(:,1),[I1_min, I1_max],degree,I1_be);
[p2_be]=Bspline_generator(control_pts(:,2),[I2_min, I2_max],degree,I2_be);


% Run the constitutive model
if n_ut > 0
    [Pdata_ut,~,wdI1_ut,wdI2_ut] = data_driven_function(F_ut,I1_ut,I2_ut,b_ut,n_ut,p1_ut,p2_ut);
else
    Pdata_ut = [];
end
if n_et > 0
    [Pdata_et,~,wdI1_et,wdI2_et] = data_driven_function(F_et,I1_et,I2_et,b_et,n_et,p1_et,p2_et);
else
    Pdata_et = [];
end
if n_ps > 0
    [Pdata_ps,~,wdI1_ps,wdI2_ps] = data_driven_function(F_ps,I1_ps,I2_ps,b_ps,n_ps,p1_ps,p2_ps);
else
    Pdata_ps = [];
end
if n_be>0
    [Pdata1_be,Pdata2_be,wdI1_be,wdI2_be] = data_driven_function(F_be,I1_be,I2_be,b_be,n_be,p1_be,p2_be);
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
