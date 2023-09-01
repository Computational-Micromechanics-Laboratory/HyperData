function [I1_ut, I1_et, I1_ps] = I1_calc(b_ut,b_et,b_ps,n_ut,n_et,n_ps)

I1_ut = zeros(n_ut,1);
I1_et = zeros(n_et,1);
I1_ps = zeros(n_ps,1);

for m = 1:n_ut
   I1_ut(m,1) = trace(b_ut(:,:,m));
end

for m = 1:n_et
   I1_et(m,1) = trace(b_et(:,:,m));
end

for m = 1:n_ps
   I1_ps(m,1) = trace(b_ps(:,:,m));
end