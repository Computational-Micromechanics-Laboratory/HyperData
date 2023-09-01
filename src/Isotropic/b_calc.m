function [b_ut, b_et, b_ps] = b_calc(F_ut,F_et,F_ps,n_ut,n_et,n_ps)

b_ut = zeros(3,3,n_ut);
b_et = zeros(3,3,n_et);
b_ps = zeros(3,3,n_ps);

for m = 1:n_ut
   b_ut(:,:,m) = F_ut(:,:,m)*transpose(F_ut(:,:,m));
end
for m = 1:n_et
   b_et(:,:,m) = F_et(:,:,m)*transpose(F_et(:,:,m));
end
for m = 1:n_ps
   b_ps(:,:,m) = F_ps(:,:,m)*transpose(F_ps(:,:,m));
end