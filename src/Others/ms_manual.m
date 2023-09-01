function [x] = ms_manual(err_ut,err_et,err_ps,err_be,w_ut,w_et,w_ps,w_be,experiment)

err_tot = w_ut*err_ut + w_et*err_et + w_ps*err_ps + w_be*err_be;
err_tot = matlabFunction(err_tot);
err_tot = @(c)err_tot(c(1),c(2),c(3),c(4),c(5));

  x0 = [0.292 22.01 1.472 0.744 0.1086];
  LB = [0 0 0 0 0];
  UB = [inf inf inf inf inf];

 x = fmincon(err_tot,x0,[],[],[],[],LB,UB);