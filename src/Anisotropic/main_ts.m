%======================================================================%
%           Computational Micromechanics Laboratory (CMML)             %
%              Middle East Technical University (METU)                 %
%----------------------------------------------------------------------%
%                        Hyper-Data Toolbox - v1.0                     %
%----------------------------------------------------------------------%
%   Authors: Recep Durna, Alp Kağan Açan, Oğuz Ziya Tikenoğulları      %
%   Software License: GNU GPL v3.0   Copyright (C) 2023 Hüsnü Dal      %
%======================================================================%

% Calculate deformation gradient tensor for 6 different shear mode.;
[F_fs,I1,~,I4f,I8fs]= shear_defo(shear,12);
[F_sf,~,I4s,~,~]    = shear_defo(shear,21);
[F_fn,~,~,~,~]      = shear_defo(shear,13);
[F_nf,~,~,~,~]      = shear_defo(shear,31);
[F_sn,~,~,~,~]      = shear_defo(shear,23);
[F_ns,~,~,~,~]      = shear_defo(shear,32);

%============ Define the objective function ============%
control_pts = sym('c',[n_ctrl_pts n_free_ener]);  % Control point vertices
% Error functions
[err_fs, err_sf,err_fn, err_nf,err_sn, err_ns,~,~,~,~,~,~,~,~ ] = ...
    object_TS(control_pts,degree,len,F_fs,F_sf,F_fn,F_nf,F_sn,F_ns,Pexp_fs,Pexp_sf,Pexp_fn,Pexp_nf,Pexp_sn,Pexp_ns,I1,I4s,I4f,I8fs);
%============ Setup the optimization problem ============%
% Weight of data in each experiment
w_1=1/6;
w_2=1/6;
w_3=1/6;
w_4=1/6;
w_5=1/6;
w_6=1/6;
% Total error
err_tot = w_1*err_fs + w_2*err_sf + w_3*err_fn + w_4*err_nf + + w_5*err_sn + w_6*err_ns ;
err_tot = matlabFunction(err_tot,'vars',{control_pts}) ; 

% Initial guess of the control pts C 
init_control_pts=ones(n_ctrl_pts,n_free_ener);
% Lower and upper bounds
LB = zeros(n_ctrl_pts,n_free_ener);
UB = inf*ones(n_ctrl_pts,n_free_ener);
% Polyconvexity and normalization?? constraints
[A,b]     = convexity_constraint(n_free_ener,n_ctrl_pts);
[Aeq,beq] = zero_stress_constraint(n_free_ener,n_ctrl_pts);
% Optimization
opt_control_pts = fmincon(err_tot,init_control_pts,A,b,Aeq,beq,LB,UB);

%============ Results at the end of optimization ============%
% Stresses and first derivatives or B-splines
[~,~,~,~,~,~,P_fs,P_sf,P_fn,P_nf,P_sn,P_ns, p1, pf, ps, pfs] = ...
    object_TS(opt_control_pts,degree,len,F_fs,F_sf,F_fn,F_nf,F_sn,F_ns,Pexp_fs,Pexp_sf,Pexp_fn,Pexp_nf,Pexp_sn,Pexp_ns,I1,I4s,I4f,I8fs);

% Quality of fit metrics
qof_fs = qof_function(shear, P_fs, Pexp_fs);
qof_sf = qof_function(shear, P_sf, Pexp_sf);
qof_fn = qof_function(shear, P_fn, Pexp_fn);
qof_nf = qof_function(shear, P_nf, Pexp_nf);
qof_sn = qof_function(shear, P_sn, Pexp_sn);
qof_ns = qof_function(shear, P_ns, Pexp_ns);
qof_tot = qof_fs+qof_sf+qof_fn+qof_nf+qof_sn+qof_ns;

% Stresses and first derivatives or B-splines / extended
[F_fs,I1,~,I4f,I8fs]= shear_defo(shear_ext,12);
[F_sf,~,I4s,~,~]    = shear_defo(shear_ext,21);
[F_fn,~,~,~,~]      = shear_defo(shear_ext,13);
[F_nf,~,~,~,~]      = shear_defo(shear_ext,31);
[F_sn,~,~,~,~]      = shear_defo(shear_ext,23);
[F_ns,~,~,~,~]      = shear_defo(shear_ext,32);
[~,~,~,~,~,~,P_fs,P_sf,P_fn,P_nf,P_sn,P_ns, p1, pf, ps, pfs] = ...
    object_TS(opt_control_pts,degree,len_ext,F_fs,F_sf,F_fn,F_nf,F_sn,F_ns,0,0,0,0,0,0,I1,I4s,I4f,I8fs);

