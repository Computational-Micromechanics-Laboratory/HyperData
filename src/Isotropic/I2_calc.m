function [I2_ut, I2_et, I2_ps] = I2_calc(b_ut,b_et,b_ps,n_ut,n_et,n_ps)

I2_ut = zeros(n_ut,1);
I2_et = zeros(n_et,1);
I2_ps = zeros(n_ps,1);

for m = 1:n_ut
   I2_ut(m,1) = (1/2)*((trace(b_ut(:,:,m)))^2 - trace(b_ut(:,:,m)*b_ut(:,:,m)));
end

for m = 1:n_et
    I2_et(m,1) = (1/2)*((trace(b_et(:,:,m)))^2 - trace(b_et(:,:,m)*b_et(:,:,m)));
end

for m = 1:n_ps
   I2_ps(m,1) = (1/2)*((trace(b_ps(:,:,m)))^2 - trace(b_ps(:,:,m)*b_ps(:,:,m)));
end
