%======================================================================%
%           Computational Micromechanics Laboratory (CMML)             %
%              Middle East Technical University (METU)                 %
%----------------------------------------------------------------------%
%                        Hyper-Data Toolbox - v1.0                     %
%----------------------------------------------------------------------%
%   Authors: Recep Durna, Alp Kağan Açan, Oğuz Ziya Tikenoğulları      %
%   Software License: GNU GPL v3.0   Copyright (C) 2023 Hüsnü Dal      %
%======================================================================%

function [I1, I4f, F, Ifib] = biax_defo_fs(lambda1, lambda2)
% Calculate  the metrics ; knowing lambda(stretch) for biaxial stretch test.
global n_free_ener fib kappafib

if (length(lambda1) ~= length(lambda2))
    error('Array length of stretch data in 1 and 2 directions are not equal')
end

len=length(lambda1);
% Initialize
F= zeros(3,3,len);
C= zeros(3,3,len);
I1= zeros(len,1);
I4f= zeros(len,1);
Ifib_n= zeros(len,n_free_ener-1);
Ifib= zeros(len,n_free_ener-1);

for i=1:len
    % Deformation gradient
    F(:,:,i) = [lambda1(i)  0           0;...
                0           lambda2(i)  0;...
                0           0           1/(lambda1(i)*lambda2(i))];
    % Right Cauchy-Green tensor
    C(:,:,i) = F(:,:,i)'*F(:,:,i);
    % First invariant
    I1(i)= C(1,1,i)+C(2,2,i)+C(3,3,i);
    % Anisotropic invariants
    for p = 1:n_free_ener-1
        for j = 1:3
            for k = 1:3
                Ifib_n(i,p) = Ifib_n(i,p)+fib(j,p)*C(j,k,i)*fib(k,p);
            end
        end
    end

    for p = 1:n_free_ener-1
        Ifib(i,p)= kappafib(p)*I1(i) + (1-3*kappafib(p))*Ifib_n(i,p);
    end
    
end

end

