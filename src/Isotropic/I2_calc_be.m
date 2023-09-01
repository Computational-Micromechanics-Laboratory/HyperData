function I2_be = I2_calc_be(b_be,n_be)

I2_be = zeros(n_be,1);

for m = 1:n_be
   I2_be(m,1) = (1/2)*((trace(b_be(:,:,m)))^2 - trace(b_be(:,:,m)*b_be(:,:,m)));
end
