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
% Symbolic expressions will be used in data_driven function to optimize the
% parameters
function [P11, P22] = data_driven3_function(F,lambda,lambda_inv,b,s,p1,p2)
% %----------------------------------------------------------%

% ----------- stretch part ------------%
lambda1 = lambda(:,1);
lambda2 = lambda(:,2);
lambda3 = lambda(:,3);

%------ 1/stretch part-------------------% 
lambda1_inv = lambda_inv(:,1);                  
lambda2_inv = lambda_inv(:,2);
lambda3_inv = lambda_inv(:,3);

% ---------------------------------------------%

% Derivatives of free energy function wrt. lambda1, lambda2 and lambda3

    wdL1 = p1(:,1);   % derivative w.r.t. lambda1
    wdL2 = p1(:,2);
    wdL3 = p1(:,3);
% Derivatives of free energy function wrt. 1/lambda1, 1/lambda2 and 1/lambda3 

    wdL1_inv = p2(:,1); % derivative w.r.t. nu1 (nu1=1/lambda1)
    wdL2_inv = p2(:,2);
    wdL3_inv = p2(:,3);
%-------------------------------------------------------------------%
% Computation of isochoric part of Kirchhoff stress         %GÜNCELLEME FAD
taubar = sym(zeros(3,3,s));
for m = 1:s
    taubar(1,1,m) =wdL1(m)*lambda1(m)- wdL1_inv(m)*lambda1_inv(m);
    taubar(2,2,m) =wdL2(m)*lambda2(m)- wdL2_inv(m)*lambda2_inv(m);
    taubar(3,3,m) =wdL3(m)*lambda3(m)- wdL3_inv(m)*lambda3_inv(m);
end

% Calculation of inverse of deformation gradient    
invF = zeros(3,3,s);
for m = 1:s
    invF(:,:,m) = inv(F(:,:,m));
end

% Computing Hydrostatic pressure term
phydr = sym(zeros(s,1));
for m = 1:s
    phydr(m,1) = -taubar(3,3,m);
end

% Kirchhoff stress,(Tau), is calculated using pressure term and Taubar
tau = sym(zeros(3,3,s));
for m = 1:s
    tau(:,:,m) = phydr(m,1)*eye(3) + taubar(:,:,m);
end

% Computing the 1st Piola stress and P11 terms
P = sym(zeros(3,3,s));
P11 = sym(zeros(s,1));
P22 = sym(zeros(s,1));
for m = 1:s
    P(:,:,m) = tau(:,:,m)*invF(:,:,m); % 1st Piola stress
    P11(m,1) = P(1,1,m);               % P11 term
    P22(m,1) = P(2,2,m);               % P22 term
end
 
% % Mean-square error expression
% err_11 = 1/s * sum((P11-P_exp).^2);
% err_22 = 1/s * sum((P22-P_exp).^2);
