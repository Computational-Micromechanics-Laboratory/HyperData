%======================================================================%
%           Computational Micromechanics Laboratory (CMML)             %
%              Middle East Technical University (METU)                 %
%----------------------------------------------------------------------%
%                        Hyper-Data Toolbox - v1.0                     %
%----------------------------------------------------------------------%
%   Authors: Recep Durna, Alp Kağan Açan, Oğuz Ziya Tikenoğulları      %
%   Software License: GNU GPL v3.0   Copyright (C) 2023 Hüsnü Dal      %
%======================================================================%

% This function will calculate the error function SYMBOLICALLY and the 
% Symbolic expressions will be used in datadriven function to optimize the
% parameters
function [err_11, err_22, P11, P22, p1, pfib] = ...
    object_ET(control_pts,degree,len,F,P_et1,P_et2,I1,I4,Ifib)
global dispmodel n_free_ener
% Calculations according to:
% niestrawska et al 2016 microstructure 
% holzapfel et al 2015 modelling 

% Calculate Ixmax and Ixmin values. If stretching tests are done in various
% ranges of Ix, then use the smallest value of Ix for Ixmin and largest one
% for the Ixmax. Ix can be I1 I4 E4 etc.
I1max= max(I1); I1min= min(I1);
I4max= max(I4); I4min= min(I4); 

for k = 1:n_free_ener-1
    Ifibmax(k) = max(Ifib(:,k));
    Ifibmin(k) = min(Ifib(:,k));
    % Generate B-splines
    pfib(:,k)=Bspline_generator(control_pts(:,k+1),[Ifibmin(k),Ifibmax(k)],degree,Ifib(:,k));
end
%generate bsplines
p1=Bspline_generator(control_pts(:,1),[I1min,I1max],degree,I1);
% p4=Bspline_generator(control_pts(:,2),[I4min,I4max],degree,I4);


% Run the constitutive model
switch dispmodel
    case "io"
        [P11, P22] = io_stress(p1, pfib, F, len);
    case "fs"
        [sigma11, sigma22, sigma33, P11, P22, P33] = fs_stress(p1, pfib, F, len);
end

% Calculate mean square error
err_11 = 1/len * sum((P11-P_et1).^2);
err_22 = 1/len * sum((P22-P_et2).^2);

end