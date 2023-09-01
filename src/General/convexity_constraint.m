%======================================================================%
%           Computational Micromechanics Laboratory (CMML)             %
%              Middle East Technical University (METU)                 %
%----------------------------------------------------------------------%
%                        Hyper-Data Toolbox - v1.0                     %
%----------------------------------------------------------------------%
%   Authors: Recep Durna, Alp Kağan Açan, Oğuz Ziya Tikenoğulları      %
%   Software License: GNU GPL v3.0   Copyright (C) 2023 Hüsnü Dal      %
%======================================================================%

function [A,b]=convexity_constraint(n_free,n_ctrl_pts)

% this function generates the A and b matrices for inequality constraint:
% such that A*x <= b. 

% control points of psi_1 psi_4 etc need to be in increasing order, to
% satisfy the convexity of the free energy function psi. 

% for example for 4 control points, constraints look like:
% 1- control_pts(1) - control_pts(2) <=0 
% 2- control_pts(2) - control_pts(3) <=0 
% 3- control_pts(3) - control_pts(4) <=0 

%gives us 4-1=3 inequality constraints. A is 3x4 and b is 3x1
% A= [1 -1 0 0; 0 1 -1 0; 0 0 1 -1] 
% b= [0; 0; 0]

A=zeros( n_free*n_ctrl_pts-n_free , n_free*n_ctrl_pts );

k=0; l=0;

for i=1:n_free 
    for j=1:n_ctrl_pts-1
        
        k=k+1;
        l=l+1;
        
        A(k,l)=1;
        A(k,l+1)=-1;
    end
    l=l+1; 
end

b=zeros(n_free*n_ctrl_pts-n_free,1);

end