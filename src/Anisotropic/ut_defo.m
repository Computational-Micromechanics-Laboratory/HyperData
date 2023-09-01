%======================================================================%
%           Computational Micromechanics Laboratory (CMML)             %
%              Middle East Technical University (METU)                 %
%----------------------------------------------------------------------%
%                        Hyper-Data Toolbox - v1.0                     %
%----------------------------------------------------------------------%
%   Authors: Recep Durna, Alp Kağan Açan, Oğuz Ziya Tikenoğulları      %
%   Software License: GNU GPL v3.0   Copyright (C) 2023 Hüsnü Dal      %
%======================================================================%

function [I1_ut1,I1_ut2,F_ut1,F_ut2,E4_ut1,E4_ut2,In_ut1,In_ut2, ...
    Ifib_ut1,Ifib_ut2,Efib_ut1,Efib_ut2] = ut_defo(lam_ut1,lam_ut2)

global fib kappafib n_free_ener dispmodel
global kip kop


A = 2*kip*kop;
B = 2*kop*(1-2*kip);

% Initialization
F_ut1 = zeros(3,3,length(lam_ut1));
F_ut2 = zeros(3,3,length(lam_ut2));
C_ut1= zeros(3,3,length(lam_ut1));
I1_ut1= zeros(length(lam_ut1),1);
I4_ut1= zeros(length(lam_ut1),1);
In_ut1= zeros(length(lam_ut1),1);
E4_ut1 = zeros(length(lam_ut1),1);
Ifib_ut1 = zeros(length(lam_ut1),n_free_ener-1);
Efib_ut1 = zeros(length(lam_ut1),n_free_ener-1);

for i = 1:length(lam_ut1)
    % Deformation gradient
    F_ut1(:,:,i) = [lam_ut1(i) 0 0;0 1/sqrt(lam_ut1(i)) 0; 0 0 1/sqrt(lam_ut1(i))];
    % Right Cauchy-Green tensor
    C_ut1(:,:,i) = F_ut1(:,:,i)'*F_ut1(:,:,i);
    % First invariant
    I1_ut1(i)= C_ut1(1,1,i)+C_ut1(2,2,i)+C_ut1(3,3,i);
%     I4_ut1(i)= C_ut1(1,1,i)*(cos(g))^2 + C_ut1(1,1,i) * (sin(g))^2 ;
    % Anisotropic invariants
    for j = 1:3
        for k = 1:3
            for p = 1:n_free_ener-1
                Ifib_ut1(i,p) = Ifib_ut1(i,p) + fib(j,p)*C_ut1(j,k,i)*fib(k,p);
            end
        end
    end

    In_ut1(i)= C_ut1(3,3,i);

    switch dispmodel
        case 'io'
            Efib_ut1(i,:) = A*I1_ut1(i) + B*Ifib_ut1(i,:) + (1-3*A-B)*In_ut1(i) -1;
%             E4_ut1(i) = A*I1_ut1(i) + B*I4_ut1(i) + (1-3*A-B)*In_ut1(i) -1;
        case 'fs'
            for p = 1:n_free_ener-1
                Ifib_ut1(i,p)= kappafib(p)*I1_ut1(i) + (1-3*kappafib(p))*Ifib_ut1(i,p);
            end
    end

end

% Initialization
C_ut2= zeros(3,3,length(lam_ut2));
I1_ut2= zeros(length(lam_ut2),1);
I4_ut2= zeros(length(lam_ut2),1);
In_ut2= zeros(length(lam_ut2),1);
E4_ut2 = zeros(length(lam_ut2),1);
Ifib_ut2 = zeros(length(lam_ut2),n_free_ener-1);
Efib_ut2 = zeros(length(lam_ut2),n_free_ener-1);

for i = 1:length(lam_ut2)
    % Deformation gradient
    F_ut2(:,:,i) = [1/sqrt(lam_ut2(i)) 0 0;0 lam_ut2(i) 0; 0 0 1/sqrt(lam_ut2(i))];
    % Right Cauchy-Green tensor
    C_ut2(:,:,i) = F_ut2(:,:,i)'*F_ut2(:,:,i);
    % First invariant
    I1_ut2(i)= C_ut2(1,1,i)+C_ut2(2,2,i)+C_ut2(3,3,i);
%     I4_ut2(i)= C_ut2(1,1,i)*(cos(g))^2 + C_ut2(1,1,i) * (sin(g))^2 ;

    % Anisotropic invariants
    for j = 1:3
        for k = 1:3
            for p = 1:n_free_ener-1
                Ifib_ut2(i,p) = Ifib_ut2(i,p) + fib(j,p)*C_ut2(j,k,i)*fib(k,p);
            end
        end
    end

    In_ut2(i)= C_ut2(3,3,i);

    switch dispmodel
        case 'io'
            Efib_ut2(i,:) = A*I1_ut2(i) + B*Ifib_ut2(i,:) + (1-3*A-B)*In_ut2(i) -1;
%             E4_ut2(i) = A*I1_ut2(i) + B*I4_ut2(i) + (1-3*A-B)*In_ut2(i) -1;
        case 'fs'
            for p = 1:n_free_ener-1
                Ifib_ut2(i,p)= kappafib(p)*I1_ut2(i) + (1-3*kappafib(p))*Ifib_ut2(i,p);
            end
    end

    
end



