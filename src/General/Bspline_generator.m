%======================================================================%
%           Computational Micromechanics Laboratory (CMML)             %
%              Middle East Technical University (METU)                 %
%----------------------------------------------------------------------%
%                        Hyper-Data Toolbox - v1.0                     %
%----------------------------------------------------------------------%
%   Authors: Recep Durna, Alp Kağan Açan, Oğuz Ziya Tikenoğulları      %
%   Software License: GNU GPL v3.0   Copyright (C) 2023 Hüsnü Dal      %
%======================================================================%

function [S,dS] = Bspline_generator(control_pts,limits,k,x)

% knots are clamped at the ends .
%
% S is the Bspline, 
% dS is the derivative of the bspline
% IS is the integral of the Bspline. according to:
%   de boor et al (1976), On Calculating with B-Splines II. Integration

%k is the degree
%n is the number of control points
%m is the number of knot points =k+n+1
%n-k is the number of Bezier curves
%xmax is the upper bound of the domain such that Domain=[0,xmax]

tol=1e-6;
len=length(x);

xmin= limits(1); 
xmax= limits(2);
if xmin == xmax
    xmax = xmax + tol;
end

n= size(control_pts,1);

m=k+n+1;
knots= zeros(m,1);
n_nonzeroknots = m - 2*(k+1);
knots(1:(k+1))= xmin;
knots((k+1):(k+2+n_nonzeroknots)) = linspace(xmin,xmax,n_nonzeroknots+2);
knots((end-k):end)= xmax+tol;

S=zeros(len,1);
dS=zeros(len,1);
% IS=zeros(len,1);

for i=1:(n)
    [N_i,dN_i] = Basis_generator(i,k,x,knots);
    %     [IN_i, ~ ] = Basis_generator(i,k+1,x,knots);
    
    S_i = N_i  * control_pts(i,:);
    S   = S    + S_i;
    
    dS_i= dN_i * control_pts(i,:);
    dS  = dS   + dS_i;
    
    %     for j=1:i
    %         control_pts_for_integral = knots
    %     end
    %     IS_i= IN_i * ( ,.. ) ;
    %     IS  = IS   + IS_i;
end

end

