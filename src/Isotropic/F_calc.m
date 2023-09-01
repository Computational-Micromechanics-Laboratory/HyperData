function [F_ut, F_et, F_ps] = F_calc(lam_ut,lam_et,lam_ps,n_ut,n_et,n_ps)

F_ut = zeros(3,3,n_ut);
F_et = zeros(3,3,n_et);
F_ps = zeros(3,3,n_ps);

for m = 1:n_ut
    F_ut(:,:,m) = [lam_ut(m) 0 0;0 1/sqrt(lam_ut(m)) 0; 0 0 1/sqrt(lam_ut(m))];
end

for m = 1:n_et
    F_et(:,:,m) = [lam_et(m) 0 0;0 lam_et(m) 0; 0 0 1/lam_et(m)^2];
end

for m = 1:n_ps
    F_ps(:,:,m) = [lam_ps(m) 0 0;0 1 0; 0 0 1/lam_ps(m)];
end