function I1_be = I1_calc_be(b_be,n_be)

I1_be = zeros(n_be,1);

for m = 1:n_be
   I1_be(m,1) = trace(b_be(:,:,m));
end