function [x] = ogden_manual(err_ut,err_et,err_ps,err_be,w_ut,w_et,w_ps,w_be,experiment)


err_tot = w_ut*err_ut + w_et*err_et + w_ps*err_ps + w_be*err_be;
err_tot = matlabFunction(err_tot);
err_tot = @(c)err_tot(c(1),c(2),c(3),c(4),c(5),c(6));

  x0 = [1.3 5 -2 0.63 1.2e-3 -1e-2];
  LB = [-inf -inf -inf -inf -inf -inf];
  UB = [inf inf inf inf inf inf];

 x = fmincon(err_tot,x0,[],[],[],[],LB,UB);
%  %% Multistart version
% 
% problem = createOptimProblem('fmincon','objective',...
%     err_tot,'x0',x0,'lb',LB,'ub',UB);
% ms = MultiStart;
% x = run(ms,problem,1000);