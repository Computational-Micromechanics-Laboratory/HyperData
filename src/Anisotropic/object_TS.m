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
% Objective function for triaxial shear data
function [err_fs, err_sf,err_fn, err_nf,err_sn, err_ns,P_fs,P_sf,P_fn,P_nf,P_sn,P_ns, p1, pf, ps, pfs] = ...
    object_TS(control_pts,degree,len,F_fs,F_sf,F_fn,F_nf,F_sn,F_ns,Pexp_fs,Pexp_sf,Pexp_fn,Pexp_nf,Pexp_sn,Pexp_ns,I1,I4s,I4f,I8fs)

% Calculations according to:
% niestrawska et al 2016 microstructure 
% holzapfel et al 2015 modelling 

% Boundaries
I1max= max(I1); I1min= min(I1);
I4fmax= max(I4f); I4fmin= min(I4f);
I4smax= max(I4s); I4smin= min(I4s);
I8fsmax= max(I8fs); I8fsmin= min(I8fs);

% Generate B-splines
p1 =Bspline_generator(control_pts(:,1),[I1min,I1max],    degree,I1);
pf =Bspline_generator(control_pts(:,2),[I4fmin,I4fmax],  degree,I4f);
ps =Bspline_generator(control_pts(:,3),[I4smin,I4smax],  degree,I4s);
pfs=Bspline_generator(control_pts(:,4),[I8fsmin,I8fsmax],degree,I8fs);

% Run the constitutive model
[P_fs,~] = triaxial_stress(p1,pf,ps,pfs,F_fs,len,12) ;
[P_sf,~] = triaxial_stress(p1,pf,ps,pfs,F_sf,len,21) ;
[P_fn,~] = triaxial_stress(p1,pf,ps,pfs,F_fn,len,13) ;
[P_nf,~] = triaxial_stress(p1,pf,ps,pfs,F_nf,len,31) ;
[P_sn,~] = triaxial_stress(p1,pf,ps,pfs,F_sn,len,23) ;
[P_ns,~] = triaxial_stress(p1,pf,ps,pfs,F_ns,len,32) ;

% Calculate mean square error
err_fs = 1/len * sum((P_fs-Pexp_fs).^2);
err_sf = 1/len * sum((P_sf-Pexp_sf).^2);
err_fn = 1/len * sum((P_fn-Pexp_fn).^2);
err_nf = 1/len * sum((P_nf-Pexp_nf).^2);
err_sn = 1/len * sum((P_sn-Pexp_sn).^2);
err_ns = 1/len * sum((P_ns-Pexp_ns).^2);

end
