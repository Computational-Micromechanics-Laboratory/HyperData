% This function will calculate the error function NUMERICALLY accoring 
% to the parameters obtained through optimization
function [P11, taubar, tau, err_11, err_22, P22] = abm_function_(F,b,s,P_exp,I1,I2,matset)

mu = matset(1);  % Regular eight-chain parameter mu
N  = matset(2);  % Regular eight-chain parameter N
mu2 = matset(3); % Constraining parameter related to area contraction
nu = (I2/3).^(1/3); % Area contraction
    
% Computaion of isochoric part of Kirchhoff stress
taubar = zeros(3,3,s);

for m = 1:s
    ln2 = (F(1,1,m)^2 + F(2,2,m)^2 + F(3,3,m)^2)/3.0;
    taubar(:,:,m) = (mu/3.0)*((3.0*N-ln2)/(N-ln2))*b(:,:,m);
end

for m=1:s
    taubar(:,:,m) = taubar(:,:,m) + (2/9)*(mu2/((nu(m))^(2)))*(I1(m)*b(:,:,m)-b(:,:,m)*b(:,:,m));
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