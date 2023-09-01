%======================================================================%
%           Computational Micromechanics Laboratory (CMML)             %
%              Middle East Technical University (METU)                 %
%----------------------------------------------------------------------%
%                        Hyper-Data Toolbox - v1.0                     %
%----------------------------------------------------------------------%
%   Authors: Recep Durna, Alp Kağan Açan, Oğuz Ziya Tikenoğulları      %
%   Software License: GNU GPL v3.0   Copyright (C) 2023 Hüsnü Dal      %
%======================================================================%

function [I1, E4, In, F, Ifib, Efib] = biax_defo_io(lambda1, lambda2)
% Calculate  the metrics ; knowing lambda(stretch) for biaxial stretch test.
global fib n_free_ener
global kip
global kop

if (length(lambda1) ~= length(lambda2))
    error('Array length of stretch data in 1 and 2 directions are not equal')
end

len=length(lambda1);

% Initialize
F= zeros(3,3,len);
C= zeros(3,3,len);
I1= zeros(len,1);
I4= zeros(len,1);
In= zeros(len,1);
E4 = zeros(len,1);
E6 = zeros(len,1);
Ifib = zeros(len,n_free_ener-1);
Efib = zeros(len,n_free_ener-1);

A = 2*kip*kop;
B = 2*kop*(1-2*kip);


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
    for j = 1:3
        for k = 1:3
            for p = 1:n_free_ener-1
                Ifib(i,p) = Ifib(i,p) + fib(j,p)*C(j,k,i)*fib(k,p);
            end
        end
    end
    
    In(i)= C(3,3,i);

    Efib(i,:) = A*I1(i) + B*Ifib(i,:) + (1-3*A-B)*In(i) -1;

%     E4(i) = A*I1(i) + B*I4(i) + (1-3*A-B)*In(i) -1;
%     E6(i) = A*I1(i) + B*I4(i) + (1-3*A-B)*In(i) -1;
end

end

