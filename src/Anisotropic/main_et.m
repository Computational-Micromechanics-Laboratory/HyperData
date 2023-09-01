%======================================================================%
%           Computational Micromechanics Laboratory (CMML)             %
%              Middle East Technical University (METU)                 %
%----------------------------------------------------------------------%
%                        Hyper-Data Toolbox - v1.0                     %
%----------------------------------------------------------------------%
%   Authors: Recep Durna, Alp Kağan Açan, Oğuz Ziya Tikenoğulları      %
%   Software License: GNU GPL v3.0   Copyright (C) 2023 Hüsnü Dal      %
%======================================================================%

% Calculate metrics of deformation for the special case of ET.
switch dispmodel
    case "io"
        [I1_ET, E4_ET, In_ET, F_ET, Ifib_ET, Efib_ET] = biax_defo_io(lam_et,lam_et);
    case "fs"
        [I1_ET, I4f_ET, F_ET, Ifib_ET] = biax_defo_fs(lam_et,lam_et);
end

%============ Define the objective function ============%
control_pts = sym('c',[n_ctrl_pts n_free_ener]);  % Control point vertices
% Error functions
switch dispmodel
    case "io"
        [err_11, err_22, ~,~,~,~] = ...
            object_ET(control_pts,degree,len,F_ET,P_et1,P_et2,I1_ET,E4_ET,Efib_ET);
    case "fs"
        [err_11, err_22, ~,~,~,~] = ...
            object_ET(control_pts,degree,len,F_ET,P_et1,P_et2,I1_ET,I4f_ET,Ifib_ET);
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
            object_ET(opt_control_pts,degree,len,F_ET,P_et1,P_et2,I1_ET,E4_ET,Efib_ET);
    case "fs"
        [~, ~, P_11, P_22, p1, p4] = ...
            object_ET(opt_control_pts,degree,len,F_ET,P_et1,P_et2,I1_ET,I4f_ET,Ifib_ET);
end

% Quality of fit metrics
qof1 = qof_function(lam_et, P_11, P_et1); 
qof2 = qof_function(lam_et, P_22, P_et2);
qof_tot = qof1+qof2;



