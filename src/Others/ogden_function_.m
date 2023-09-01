% This function will calculate the error function NUMERICALLY accoring 
% to the parameters obtained through optimization
function [P11, taubar, tau, err_11, err_22, P22] = ogden_function_(F,s,P_exp,matset)

% material parameters 
    a1  = matset(1);
    a2  = matset(2);
    a3  = matset(3);
    mu1 = matset(4);
    mu2 = matset(5);
    mu3 = matset(6);
 
% Derivatives if free energy function wrt. I1 and I2
wdL = zeros(s,3,1);

a = [a1 a2 a3];
mu = [mu1 mu2 mu3];
na = [1 0 0; 0 1 0; 0 0 1];

    for m = 1:s
        for i = 1:3
            for j = 1:3                
            wdL(m,i,1) = wdL(m,i,1) + mu(j)*((F(i,i,m))^(a(j)-1));    
            end
        end
    end

% Computaion of isochoric part of Kirchhoff stress
taubar = zeros(3,3,s);
for m = 1:s
    for i = 1:3
    taubar(:,:,m) = taubar(:,:,m) +wdL(m,i,1)*F(i,i,m)*na(i,:)'*na(i,:);
    end
end

% Calculation of inverse of deformation gradient    
invF = zeros(3,3,s);
for m = 1:s
    invF(:,:,m) = inv(F(:,:,m));
end

% Computing Hydrostatic pressure term
phydr = zeros(s,1);
for m = 1:s
    phydr(m,1) = -taubar(3,3,m);
end

% Kirchhoff stress,(Tau), is calculated using pressure term and Taubar
tau = zeros(3,3,s);
for m = 1:s
    tau(:,:,m) = phydr(m,1)*eye(3) + taubar(:,:,m);
end

% Computing the 1st Piola stress and P11 terms
P = zeros(3,3,s);
P11 = zeros(s,1);
P22 = zeros(s,1);
for m = 1:s
    P(:,:,m) = tau(:,:,m)*invF(:,:,m); % 1st Piola stress
    P11(m,1) = P(1,1,m);               % P11 term
    P22(m,1) = P(2,2,m);               % P22 term
end
 
% Mean-square error expression
err_11 = sum((P11-P_exp).^2);
err_22 = sum((P22-P_exp).^2);