%======================================================================%
%           Computational Micromechanics Laboratory (CMML)             %
%              Middle East Technical University (METU)                 %
%----------------------------------------------------------------------%
%                        Hyper-Data Toolbox - v1.0                     %
%----------------------------------------------------------------------%
%   Authors: Recep Durna, Alp Kağan Açan, Oğuz Ziya Tikenoğulları      %
%   Software License: GNU GPL v3.0   Copyright (C) 2023 Hüsnü Dal      %
%======================================================================%

function [F,I1,I4s,I4f,I8fs] = shear_defo(shear,mode)
% Calculate  the metrics ; knowing shear for triaxial shear.
global kappaf kappas

len=length(shear);
% Initialization
F= zeros(3,3,len);
C=zeros(3,3,len);
I1=zeros(len,1);
I4f = zeros(len,1);
I4s = zeros(len,1);
I8fs = zeros(len,1);

for i=1:len
    % Deformation gradient
    if mode == 12
        F(:,:,i) = eye(3) +  shear(i) * [ 0 0 0; 1 0 0 ;0 0 0];
    elseif mode == 21
        F(:,:,i) = eye(3) +  shear(i) * [ 0 1 0; 0 0 0 ;0 0 0];
    elseif mode == 13
        F(:,:,i) = eye(3) +  shear(i) * [ 0 0 0; 0 0 0 ;1 0 0];
    elseif mode == 31
        F(:,:,i) = eye(3) +  shear(i) * [ 0 0 1; 0 0 0 ;0 0 0];
    elseif mode == 23
        F(:,:,i) = eye(3) +  shear(i) * [ 0 0 0; 0 0 0 ;0 1 0];
    elseif mode == 32
        F(:,:,i) = eye(3) +  shear(i) * [ 0 0 0; 0 0 1 ;0 0 0];
    end
    % Right Cauchy-Green tensor
    C(:,:,i)= F(:,:,i)'*F(:,:,i);
    % First invariant
    I1(i)   = C(1,1,i)+ C(2,2,i)+ C(3,3,i);
    % Anisotropic invariants
    I4f(i)  = kappaf*I1(i) + (1-3*kappaf)*C(1,1,i);
    I4s(i)  = kappas*I1(i) + (1-3*kappas)*C(2,2,i);
    % Interaction term
    I8fs(i) = 0.5*(C(1,2,i)+C(2,1,i));
end

end