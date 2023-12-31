%======================================================================%
%           Computational Micromechanics Laboratory (CMML)             %
%              Middle East Technical University (METU)                 %
%----------------------------------------------------------------------%
%                        Hyper-Data Toolbox - v1.0                     %
%----------------------------------------------------------------------%
%   Authors: Recep Durna, Alp Kağan Açan, Oğuz Ziya Tikenoğulları      %
%   Software License: GNU GPL v3.0   Copyright (C) 2023 Hüsnü Dal      %
%======================================================================%

%========== Stretch based data driven formulation ==========%

%============ Define the objective function ============%
control_pts = sym('c',[n_ctrl_pts n_free_ener]);  %control point vertices

%============ Setup the optimization problem ============%
w_1=1/4;
w_2=1/4;
w_3=1/4;
w_4=1/4;

% Calculate errors
[err_ut, err_et,err_ps,err_be, ~, ~,~,~,~] = ...
   object_str(control_pts,degree,lam_ut,lam_et,lam_ps,lam1_be,lam2_be,P_ut,P_et,P_ps,P1_be,P2_be);

err_tot = w_1*err_ut +w_2*err_et + w_3*err_ps + w_4*err_be;
err_tot = matlabFunction(err_tot,'vars',{control_pts}) ; 
% Initial guess of the control pts C 
init_control_pts=ones(n_ctrl_pts,n_free_ener); 
LB = zeros(n_ctrl_pts,n_free_ener);
UB = inf*ones(n_ctrl_pts,n_free_ener);
[A,b]     = convexity_constraint(n_free_ener,n_ctrl_pts);
% options = optimset('PlotFcns',@optimplotfval,'MaxFunEvals',2e4,'MaxIter',2e4, ...
%     'TolX',1e-10);
 opt_control_pts = fmincon(err_tot,init_control_pts,A,b,[],[],LB,UB,[],[]);
%   opt_control_pts = fmincon(err_tot,init_control_pts,[],[],[],[],LB,UB,[],options);

%============ Results at the end of optimization ============%
[err_ut, err_et,err_ps,err_be, Pdata_ut, Pdata_et,Pdata_ps,Pdata1_be, Pdata2_be] = ...
   object_str(opt_control_pts,degree,lam_ut,lam_et,lam_ps,lam1_be,lam2_be,P_ut,P_et,P_ps,P1_be,P2_be);

[~, ~,~,~, Pdata_ut_ext, Pdata_et_ext,Pdata_ps_ext,Pdata1_be_ext, Pdata2_be_ext] = ...
    object_str(opt_control_pts,degree,lam_ut_ext,lam_et_ext,lam_ps_ext,lam1_be_ext,lam2_be_ext,0,0,0,0,0);


Pdata_ut = double(Pdata_ut);
Pdata_et = double(Pdata_et);
Pdata_ps = double(Pdata_ps);
Pdata1_be = double(Pdata1_be);
Pdata2_be = double(Pdata2_be);

Pdata_ut_ext = double(Pdata_ut_ext);
Pdata_et_ext = double(Pdata_et_ext);
Pdata_ps_ext = double(Pdata_ps_ext);
Pdata1_be_ext = double(Pdata1_be_ext);
Pdata2_be_ext = double(Pdata2_be_ext);