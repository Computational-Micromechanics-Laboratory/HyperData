%======================================================================%
%           Computational Micromechanics Laboratory (CMML)             %
%              Middle East Technical University (METU)                 %
%----------------------------------------------------------------------%
%                        Hyper-Data Toolbox - v1.0                     %
%----------------------------------------------------------------------%
%   Authors: Recep Durna, Alp Kağan Açan, Oğuz Ziya Tikenoğulları      %
%   Software License: GNU GPL v3.0   Copyright (C) 2023 Hüsnü Dal      %
%======================================================================%

function [A,b]=zero_stress_constraint(n_free,n_ctrl_pts)

% this function generates the A and b matrices for equality constraint:
% such that A*x = b. 

% to have zero-stress at the undeformed configuration, 
% we need to have the first control point of psi4f, psi4s etc equal to zero
% ps1 does not have to start from zero 

% for example for 4 control points and 2 free energy functions, 
% constraints look like:
% 1- control_pts(1) = (free) 
% 2- control_pts(5) = 0


% gives us "n_free=2" inequality constraints. A is 2x8 and b is 2x1
% A= [0 0 0 0 0 0 0 0; 0 0 0 0 1 0 0 0] 
% b= [0; 0]

A=zeros( n_free , n_free*n_ctrl_pts );

for i=2:n_free  %(start from 2 because psi1 doesn't have a constraint)
    A(i, n_ctrl_pts*(i-1)+1)= 1;
end

b=zeros(n_free,1);

end