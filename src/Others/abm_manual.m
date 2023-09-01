function [x] = abm_manual(err_ut,err_et,err_ps,err_be,w_ut,w_et,w_ps,w_be,experiment)

err_tot = w_ut*err_ut + w_et*err_et + w_ps*err_ps + w_be*err_be;
err_tot = matlabFunction(err_tot);
err_tot = @(c)err_tot(c(1),c(2),c(3));

 x0 = [0.32  60 2.32];
 LB = [0     0   0];
 UB = [ inf  inf inf];

 x = fmincon(err_tot,x0,[],[],[],[],LB,UB);
