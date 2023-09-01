%======================================================================%
%           Computational Micromechanics Laboratory (CMML)             %
%              Middle East Technical University (METU)                 %
%----------------------------------------------------------------------%
%                        Hyper-Data Toolbox - v1.0                     %
%----------------------------------------------------------------------%
%   Authors: Recep Durna, Alp Kağan Açan, Oğuz Ziya Tikenoğulları      %
%   Software License: GNU GPL v3.0   Copyright (C) 2023 Hüsnü Dal      %
%======================================================================%

function [sigma11, sigma22, sigma33, P11, P22, P33] = ...
    fs_stress(psi1, psi4, F, len)
% Calculate stresses according to Eriksson 2013 model. 
% assume one fiber family in 11 direction perfectly aligned with dispersion
% assume no stress boundary condition in 3-direction; such as biaxial
% loading or uniaxial loading

global fib kappafib n_free_ener tnsonlycond

% f=[1 0 0]';
% s=[0 1 0]';

if (len ~= size(F,3))
    error('size of I1 is not equal to F')
end

% fof=zeros(3,3);
% for i=1:3
%     for j=1:3
%         fof(i,j)=f(i)*f(j);
%     end
% end

% Dyadic product
fof=zeros(3,3,n_free_ener-1);
for i=1:3
    for j=1:3
        for k = 1:n_free_ener-1
            fof(i,j,k)=fib(i,k)*fib(j,k);
        end
    end
end

I       = eye(3,3);
% Generalized structure tensors
for k = 1:n_free_ener-1
    H(:,:,k)       = kappafib(k)*I + (1-3*kappafib(k))*fof(:,:,k);
end
% Symbolic initialization
S       = sym('S', [3 3 len]);
Sigma   = sym('Sigma', [3 3 len]);
P       = sym('P', [3 3 len]);

I4=ones(len,n_free_ener-1);
if tnsonlycond
    for i=1:len
	    I4(i,:) = F(1,1,i) * F(1,1,i) * fib(1,:).^2 + F(2,2,i) * F(2,2,i) * fib(2,:).^2 + F(3,3,i) * F(3,3,i) * fib(3,:).^2;
    end
end

% Stress calculations
for i=1:len
%     S(:,:,i)=2*psi1(i)*I + 2*psi4(i)*H;% -pC^(-1) is not necessary, 

    S(:,:,i)=2*psi1(i)*I;
    for k = 1:n_free_ener-1
        if I4(i,k) >= 1
            S(:,:,i)= S(:,:,i)+2*psi4(i,k)*H(:,:,k);
        end
    end

    % because of conservation of volume Kirchhoff=Cauchy:
    Sigma(:,:,i)= F(:,:,i)*S(:,:,i)*F(:,:,i)';
    % subtract the hydrostatic pressures.
    Sigma(:,:,i)= Sigma(:,:,i) - Sigma(3,3,i)*I;

    P(:,:,i)    = Sigma(:,:,i)*inv(F(:,:,i)');
    S(:,:,i)    = inv(F(:,:,i))*Sigma(:,:,i)*inv(F(:,:,i)');

end

% Return the stress components as arrays:
sigma11=reshape(Sigma(1,1,:),len,1);
sigma22=reshape(Sigma(2,2,:),len,1);
sigma33=reshape(Sigma(3,3,:),len,1);

P11=reshape(P(1,1,:),len,1);
P22=reshape(P(2,2,:),len,1);
P33=reshape(P(3,3,:),len,1);

end

