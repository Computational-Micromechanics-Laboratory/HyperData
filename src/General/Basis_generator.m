%======================================================================%
%           Computational Micromechanics Laboratory (CMML)             %
%              Middle East Technical University (METU)                 %
%----------------------------------------------------------------------%
%                        Hyper-Data Toolbox - v1.0                     %
%----------------------------------------------------------------------%
%   Authors: Recep Durna, Alp Kağan Açan, Oğuz Ziya Tikenoğulları      %
%   Software License: GNU GPL v3.0   Copyright (C) 2023 Hüsnü Dal      %
%======================================================================%

function [N,dN] = Basis_generator(i,k,x,knots)
% Generate basis functions according to cox-de boor recursion formula
tol = 1e-8;

% knots(1)    =knots(1)-1e-6;
% knots(end)  =knots(end)+1e-6;

if (k==0)
    N = ones(size(x));
    N(x < knots(i)) = 0 ;
    N(x >= knots(i+1)) = 0;
else
    N= (x-knots(i))/(knots(i+k) - knots(i)+tol) .* Basis_generator(i,k-1,x,knots)...
        + (knots(i+k+1)-x)/(knots(i+k+1) - knots(i+1)+tol) .* Basis_generator(i+1,k-1,x,knots);
    dN=(k)/(knots(i+k) - knots(i)+tol) .* Basis_generator(i,k-1,x,knots)...
        - (k)/(knots(i+k+1) - knots(i+1)+tol) .* Basis_generator(i+1,k-1,x,knots);
end
end

