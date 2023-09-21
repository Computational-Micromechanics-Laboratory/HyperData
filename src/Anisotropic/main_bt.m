%======================================================================%
%           Computational Micromechanics Laboratory (CMML)             %
%              Middle East Technical University (METU)                 %
%----------------------------------------------------------------------%
%                        Hyper-Data Toolbox - v1.0                     %
%----------------------------------------------------------------------%
%   Authors: Recep Durna, Alp Kağan Açan, Oğuz Ziya Tikenoğulları      %
%   Software License: GNU GPL v3.0   Copyright (C) 2023 Hüsnü Dal      %
%======================================================================%

% Calculate metrics of deformation for BT.
switch dispmodel
    case "io"
        [I1, E4, In, F, Ifib, Efib] = biax_defo_io(lam_1,lam_2);
    case "fs"
        [I1, I4f, F, Ifib] = biax_defo_fs(lam_1,lam_2);
end

%============ Define the objective function ============%
control_pts = sym('c',[n_ctrl_pts n_free_ener]);  % Control point vertices
% Error functions
switch dispmodel
    case "io"
        [err_11, err_22, ~,~,~,~] = ...
            object_ET(control_pts,degree,len,F,P_exp1,P_exp2,I1,E4,Efib);
    case "fs"
        [err_11, err_22, ~,~,~,~] = ...
            object_ET(control_pts,degree,len,F,P_exp1,P_exp2,I1,I4f,Ifib);
end

%============ Setup the optimization problem ============%
% Weight of data in each direction
w_1=0.5;
w_2=0.5;

% Total error
err_tot = w_1*err_11 + w_2*err_22  ;
err_tot = matlabFunction(err_tot,'vars',{control_pts}) ; 

% Initial guess of the control pts C 
init_control_pts=ones(n_ctrl_pts,n_free_ener);
% Lower and upper bounds
LB = zeros(n_ctrl_pts,n_free_ener);
UB = inf*ones(n_ctrl_pts,n_free_ener);
% Polyconvexity and normalization constraints
[A,b]     = convexity_constraint(n_free_ener,n_ctrl_pts);
[Aeq,beq] = zero_stress_constraint(n_free_ener,n_ctrl_pts);
% Optimization
opt_control_pts = fmincon(err_tot,init_control_pts,A,b,Aeq,beq,LB,UB);

%============ Results at the end of optimization ============%
% Stresses and first derivatives or B-splines
switch dispmodel
    case "io"
        [~, ~, P_11, P_22, p1, p4] = ...
            object_ET(opt_control_pts,degree,len,F,P_exp1,P_exp2,I1,E4,Efib);
    case "fs"
        [~, ~, P_11, P_22, p1, p4] = ...
            object_ET(opt_control_pts,degree,len,F,P_exp1,P_exp2,I1,I4f,Ifib);
end

% Quality of fit metrics
bnds_sz = size(bounds);
for i = 1:bnds_sz(1)
    qof1_ind(i,:) = qof_function(lam_1(bounds(i,1):bounds(i,2)), P_11(bounds(i,1):bounds(i,2)), P_exp1(bounds(i,1):bounds(i,2))); 
    qof2_ind(i,:) = qof_function(lam_2(bounds(i,1):bounds(i,2)), P_22(bounds(i,1):bounds(i,2)), P_exp2(bounds(i,1):bounds(i,2)));
    qof_tot_ind(i,:) = qof1_ind(i,:)+qof2_ind(i,:);
end
qof1 = zeros(1,4);
qof2 = zeros(1,4);
qof_tot = zeros(1,4);
for i = 1:bnds_sz(1)
    qof1 = qof1 + qof1_ind(i,:);
    qof2 = qof2 + qof2_ind(i,:);
    qof_tot = qof_tot + qof_tot_ind(i,:);
end

% Stresses and first derivatives or B-splines / extended
switch dispmodel
    case "io"
        [I1, E4, In, F, Ifib, Efib] = biax_defo_io(lam_1_ext,lam_2_ext);
    case "fs"
        [I1, I4f, F, Ifib] = biax_defo_fs(lam_1_ext,lam_2_ext);
end
switch dispmodel
    case "io"
        [~, ~, P_11, P_22, p1, p4] = ...
            object_ET(opt_control_pts,degree,len_ext,F,0,0,I1,E4,Efib);
    case "fs"
        [~, ~, P_11, P_22, p1, p4] = ...
            object_ET(opt_control_pts,degree,len_ext,F,0,0,I1,I4f,Ifib);
end


