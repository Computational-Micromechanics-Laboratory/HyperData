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
% Symbolic expressions will be used in data_driven2 function to optimize the
% parametersobject_modinv
function [P11,P22,wdI1,wdI2] = data_driven2_function(F,I1,I2,b,s,p1,p2)

% Calculate Extended Eight Chain like parameretes lambda_n and nu_n
% and their derivatives with respect to I1 and I2
lambda_chn=zeros(s,1);
dlambda_chn=zeros(s,1);
nu_chn=zeros(s,1);
dnu_chn=zeros(s,1);

for t=1:s
    lambda_chn(t,1)=(I1(t,1)/3)^(1/2);
    nu_chn(t,1)=(I2(t,1)/3)^(1/3);
    % dlambda_chn(t,1)=(1/(12*I1(t,1)))^(1/2);
    % dnu_chn(t,1)=((1/(81*((I2(t,1))^2)))^(1/3));
    dlambda_chn(t,1)=(1/6)*(1/lambda_chn(t,1));
    dnu_chn(t,1)=(1/9)*(1/nu_chn(t,1))^2;
end

% Derivatives of free energy function wrt. I1 and I2
wdI1 = p1.*dlambda_chn;
wdI2 = p2.*dnu_chn;
    
% Computation of isochoric part of Kirchhoff stress 
taubar = sym(zeros(3,3,s));
for m = 1:s
    taubar(:,:,m) = 2*(wdI1(m,1)+I1(m,1)*wdI2(m,1))*b(:,:,m,1)-...
                    2*wdI2(m,1)*(b(:,:,m))^2;
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

