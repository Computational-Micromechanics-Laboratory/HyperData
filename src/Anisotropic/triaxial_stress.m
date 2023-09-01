%======================================================================%
%           Computational Micromechanics Laboratory (CMML)             %
%              Middle East Technical University (METU)                 %
%----------------------------------------------------------------------%
%                        Hyper-Data Toolbox - v1.0                     %
%----------------------------------------------------------------------%
%   Authors: Recep Durna, Alp Kağan Açan, Oğuz Ziya Tikenoğulları      %
%   Software License: GNU GPL v3.0   Copyright (C) 2023 Hüsnü Dal      %
%======================================================================%

function [Pcalc,Sigma] = ...
    triaxial_stress(psi1, psif,psis,psifs, F, len,mode)
% calculate stresses according to Eriksson 2013 model. 
% assume one fiber family in 11 direction perfectly aligned with dispersion

global kappaf kappas fib

if (len ~= size(F,3))
    error('size of I1 is not equal to F')
end
         
f0 = fib(:,1); % f direction in Lagrangian configuration
s0 = fib(:,2); % s direction in Lagrangian configuration
% n direction is [0 0 1], however it is not used
% in this code.

% Initialization
f0_tens=zeros(3,3);
s0_tens=zeros(3,3);
fs_tens=zeros(3,3);

% Dyadic products
for i=1:3
    for j=1:3
        f0_tens(i,j) = f0(i)*f0(j);
        s0_tens(i,j) = s0(i)*s0(j);
        fs_tens(i,j) = 0.5*(f0(i)*s0(j)+s0(i)*f0(j));  % In articles, it is denoted by H_{fs} = sym(f_0 \dyad s_0)
    end
end

I       = eye(3,3);
% Generalized structure tensors
Hf       = kappaf*eye(3) + (1-3*kappaf)*f0_tens ;
Hs       = kappas*eye(3) + (1-3*kappas)*s0_tens ;

% Symbolic initialization
S       = sym('S', [3 3 len]);
Sigma   = sym('Sigma', [3 3 len]);
P       = sym('P', [3 3 len]);

% Stress calculations
for i=1:len
    if mode == 12
        S(:,:,i)=2*psi1(i)*I + 2*psif(i)*Hf + 2*psifs(i)*fs_tens;    %fs 
    elseif mode == 21
        S(:,:,i)=2*psi1(i)*I + 2*psis(i)*Hs + 2*psifs(i)*fs_tens;    %sf
    elseif mode == 13
        S(:,:,i)=2*psi1(i)*I + 2*psif(i)*Hf;                         %fn
    elseif mode == 31
        S(:,:,i)=2*psi1(i)*I ;                                       %nf
    elseif mode == 23
        S(:,:,i)=2*psi1(i)*I + 2*psis(i)*Hs;                         %sn
    elseif mode == 32
        S(:,:,i)=2*psi1(i)*I;                                        %ns
    end
    % subtract the hydrostatic pressures??

    % because of conservation of volume Kirchhoff=Cauchy:
    Sigma(:,:,i)= F(:,:,i)*S(:,:,i)*F(:,:,i)';
    % PK1 stress:
     P(:,:,i)    = Sigma(:,:,i)*inv(F(:,:,i));
end

% Return the stress components as arrays:
Pcalc = sym(zeros(len,1));
for i = 1:len
    if mode == 12
        Pcalc(i,:) = Sigma(2,1,i) ;
    elseif mode == 21
         Pcalc(i,:) = Sigma(1,2,i) ;
    elseif mode == 13 
         Pcalc(i,:) = Sigma(3,1,i) ;
    elseif mode == 31
         Pcalc(i,:) = Sigma(1,3,i) ;  
    elseif mode == 23 
         Pcalc(i,:) = Sigma(3,2,i) ;
    elseif mode == 32
         Pcalc(i,:) = Sigma(2,3,i) ;       
    end
end

end
