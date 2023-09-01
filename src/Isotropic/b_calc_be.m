function b_be = b_calc_be(F_be,n_be)

b_be = zeros(3,3,n_be,1);

for m = 1:n_be
   b_be(:,:,m) = F_be(:,:,m)*transpose(F_be(:,:,m));
end
