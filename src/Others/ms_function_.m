% This function will calculate the error function NUMERICALLY accoring 
% to the parameters obtained through optimization
function [P11, taubar, tau, err_11, err_22, P22] = ms_function_(F,s,P_exp,matset)

xr=    [0,                 0,                 1                ;...
        0,                 1,                 0                ;...
        1,                 0,                 0                ;...
        0,                 0.707106781187000, 0.707106781187000;...
        0,                -0.707106781187000, 0.707106781187000;...
        0.707106781187000, 0,                 0.707106781187000;...
       -0.707106781187000, 0,                 0.707106781187000;...
        0.707106781187000, 0.707106781187000, 0                ;...
       -0.707106781187000, 0.707106781187000, 0                ;...
        0.836095596749000, 0.387907304067000, 0.387907304067000;...
       -0.836095596749000, 0.387907304067000, 0.387907304067000;...
        0.836095596749000,-0.387907304067000, 0.387907304067000;...
       -0.836095596749000,-0.387907304067000, 0.387907304067000;...
        0.387907304067000, 0.836095596749000, 0.387907304067000;...
       -0.387907304067000, 0.836095596749000, 0.387907304067000;...
        0.387907304067000,-0.836095596749000, 0.387907304067000;...
       -0.387907304067000,-0.836095596749000, 0.387907304067000;...
        0.387907304067000, 0.387907304067000, 0.836095596749000;...
       -0.387907304067000, 0.387907304067000, 0.836095596749000;...
        0.387907304067000,-0.387907304067000, 0.836095596749000;...
       -0.387907304067000,-0.387907304067000, 0.836095596749000];
% weight factors
   w= 2 * [0.0265214244093000;...
           0.0265214244093000;...
           0.0265214244093000;...
           0.0199301476312000;...
           0.0199301476312000;...
           0.0199301476312000;...
           0.0199301476312000;...
           0.0199301476312000;...
           0.0199301476312000;...
           0.0250712367487000;...
           0.0250712367487000;...
           0.0250712367487000;...
           0.0250712367487000;...
           0.0250712367487000;...
           0.0250712367487000;...
           0.0250712367487000;...
           0.0250712367487000;...
           0.0250712367487000;...
           0.0250712367487000;...
           0.0250712367487000;...
           0.0250712367487000];
       
% material parameters 
 mu = matset(1); 
 N  = matset(2);
 p  = matset(3);
 U  = matset(4);
 q  = matset(5);
       
% ----------------------- PART I -------------------------%
% ----------- compute non-affine stretch part ------------%
t = zeros(21,3,s);
for i=1:21
    for j = 1:3
        for k = 1:3
            for m = 1:s
                t(i,j,m) = t(i,j,m) + F(j,k,m) * xr(i,k);
            end
        end
    end
end
      
% compute orientation stretches
lam = zeros(21,s);
for i=1:21
    for m = 1:s
        lam(i,m) = sqrt(t(i,1,m)*t(i,1,m)+t(i,2,m)*t(i,2,m)+t(i,3,m)*t(i,3,m));
    end
end

% average network stretch
lambda= zeros(s,1);
for i=1:21
    for m = 1:s
        lambda(m) = lambda(m) + lam(i,m)^p*w(i);
    end
end

    lambda = lambda.^(1/p);
        
% compute Kirchhoff stresses
xh = zeros(3,3,s);
for k=1:21
    for i = 1:3
        for j = 1:3
            for m = 1:s
                xh(i,j,m) = xh(i,j,m) + lam(k,m)^(p-2) * t(k,i,m) * t(k,j,m) * w(k) ;
            end
        end
    end
end

    tauf = mu*lambda.*(3*N-lambda.^2)./(N-lambda.^2);
    
    tau1 = zeros(3,3,s);
    
    for m=1:s
        tau1(:,:,m) = tauf(m) * lambda(m)^(1-p) .* xh(:,:,m);
    end

        
% -------------------- PART II ------------------------%
% ----------- compute non-affine tube part ------------%
invF = zeros(3,3,s);

for m = 1:s
    invF(:,:,m) = inv(F(:,:,m));
end

xn = zeros(21,3,s);
for i=1:21
    for j = 1:3
        for k = 1:3
            for m = 1:s   
                xn(i,j,m) = xn(i,j,m) + invF(j,k,m) * xr(i,k);
            end
        end
    end
end
      
% compute orientation stretches
nu = zeros(21,s);
for i=1:21
    for m = 1:s
        nu(i,m) = sqrt(xn(i,1,m)*xn(i,1,m)+xn(i,2,m)*xn(i,2,m)+xn(i,3,m)*xn(i,3,m));
    end
end
% compute Kirchhoff stresses
xk = zeros(3,3,s);
for k=1:21
    for i = 1:3
        for j = 1:3
            for m = 1:s
                xk(i,j,m) = xk(i,j,m) + nu(k,m)^(q-2) * xn(k,i,m) * xn(k,j,m) * w(k) ;
            end
        end
    end
end
        
tau2 = - (mu*N*U*q)*xk;
        
% -------------------- Total Kirchoff Stress ------------------------%
 
taubar = tau1+tau2;

phydr = zeros(s,1);
for m = 1:s
    phydr(m,1) = -taubar(3,3,m);
end

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