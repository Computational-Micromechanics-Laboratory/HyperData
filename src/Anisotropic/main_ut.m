%======================================================================%
%           Computational Micromechanics Laboratory (CMML)             %
%              Middle East Technical University (METU)                 %
%----------------------------------------------------------------------%
%                        Hyper-Data Toolbox - v1.0                     %
%----------------------------------------------------------------------%
%   Authors: Recep Durna, Alp Kağan Açan, Oğuz Ziya Tikenoğulları      %
%   Software License: GNU GPL v3.0   Copyright (C) 2023 Hüsnü Dal      %
%======================================================================%

% Calculate metrics of deformation for the special case of UT.
[I1_ut1,I1_ut2,F_ut1,F_ut2,E4_ut1,E4_ut2,In_ut1,In_ut2,Ifib_ut1,Ifib_ut2,...
    Efib_ut1,Efib_ut2] = ut_defo(lam_ut1,lam_ut2);

%============ Define the objective function ============%
control_pts = sym('c',[n_ctrl_pts n_free_ener]);  % Control point vertices
% Error functions
switch dispmodel
    case "io"
        [err_11, ~, ~, ~, ~, ~] = ...
            object_UT(control_pts,degree,n_ut1,F_ut1,P_ut1,I1_ut1,I1_ut1,E4_ut1,E4_ut1,Efib_ut1);
        [~, err_22, ~, ~, ~, ~] = ...
            object_UT(control_pts,degree,n_ut2,F_ut2,P_ut2,I1_ut2,I1_ut2,E4_ut2,E4_ut2,Efib_ut2);
    case "fs"
%         [err_11,  ~, ~,~,~,~] = ...
%             object_UT(control_pts,degree,n_ut1,F_ut1,P_ut1,I1_ut1,I1_ut1,I4f_ut1,I4f_ut1,Ifib_ut1);
%         [~, err_22 ,~, ~,~,~,~] = ...
%             object_UT(control_pts,degree,n_ut2,F_ut2,P_ut2,I1_ut2,I1_ut2,I4f_ut2,I4f_ut2,Ifib_ut2);
        [err_11,  ~, ~,~,~,~] = ...
            object_UT(control_pts,degree,n_ut1,F_ut1,P_ut1,I1_ut1,I1_ut1,E4_ut1,E4_ut1,Ifib_ut1);
        [~, err_22 ,~, ~,~,~,~] = ...
            object_UT(control_pts,degree,n_ut2,F_ut2,P_ut2,I1_ut2,I1_ut2,E4_ut2,E4_ut2,Ifib_ut2);
end

%============ Setup the optimization problem ============%
% Weight of data for each experiment
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
% Polyconvexity and normalization?? constraints
[A,b]     = convexity_constraint(n_free_ener,n_ctrl_pts);
[Aeq,beq] = zero_stress_constraint(n_free_ener,n_ctrl_pts);
% Optimization
 opt_control_pts = fmincon(err_tot,init_control_pts,A,b,Aeq,beq,LB,UB);

%============ Results at the end of optimization ============%
% Stresses and first derivatives or B-splines
 switch dispmodel
    case "io"
        [err_11, ~, P11, ~, ~, ~,~,~] = ...
            object_UT(opt_control_pts,degree,n_ut1,F_ut1,P_ut1,I1_ut1,I1_ut1,E4_ut1,E4_ut1,Efib_ut1); %CHECK
        [~, err_22, ~, P22, p1_1, p4_1,p1_2,p4_2] = ...
            object_UT(opt_control_pts,degree,n_ut2,F_ut2,P_ut2,I1_ut2,I1_ut2,E4_ut2,E4_ut2,Efib_ut2);
    case "fs"
%         [err_11, ~, P11, ~, ~, ~,~,~] = ...
%             object_UT(opt_control_pts,degree,n_ut1,F_ut1,P_ut1,I1_ut1,I1_ut1,I4f_ut1,I4f_ut1,Ifib_ut1);
%         [~, err_22, ~, P22, p1_1, p4_1,p1_2,p4_2] = ...
%             object_UT(opt_control_pts,degree,n_ut2,F_ut2,P_ut2,I1_ut2,I1_ut2,I4f_ut2,I4f_ut2,Ifib_ut2);
        [err_11, ~, P11, ~, ~, ~,~,~] = ...
            object_UT(opt_control_pts,degree,n_ut1,F_ut1,P_ut1,I1_ut1,I1_ut1,E4_ut1,E4_ut1,Ifib_ut1);
        [~, err_22, ~, P22, p1_1, p4_1,p1_2,p4_2] = ...
            object_UT(opt_control_pts,degree,n_ut2,F_ut2,P_ut2,I1_ut2,I1_ut2,E4_ut2,E4_ut2,Ifib_ut2);
end
% Quality of fit metric
qof1 = qof_function(lam_ut1, P11, P_ut1); 
qof2 = qof_function(lam_ut2, P22, P_ut2);
qof_tot = qof1+qof2;
