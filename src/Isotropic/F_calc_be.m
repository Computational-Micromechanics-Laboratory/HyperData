function F_be = F_calc_be(lam1_be,lam2_be,n_be)

F_be = zeros(3,3,n_be);

for m = 1:n_be
    F_be(:,:,m) = [lam1_be(m) 0 0;0 lam2_be(m) 0; 0 0 1/(lam1_be(m)*lam2_be(m))];
end
