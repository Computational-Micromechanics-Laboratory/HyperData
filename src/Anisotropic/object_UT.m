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
function [err_11, err_22, P11, P22, p1_1, p4_1,p1_2,p4_2] = ...
   object_UT(control_pts,degree,len,F,P_ut,I1_1,I1_2,E4_1,E4_2,Efib)

global dispmodel n_free_ener
% Calculations according to:
% niestrawska et al 2016 microstructure 
% holzapfel et al 2015 modelling 

% Calculate Ixmax and Ixmin values. If stretching tests are done in various
% ranges of Ix, then use the smallest value of Ix for Ixmin and largest one
% for the Ixmax. Ix can be I1 I4 E4 etc.

I1max= max([I1_1;I1_2]); I1min= min([I1_1;I1_2]);
E4max= max([E4_1;E4_2]); E4min= min([E4_1;E4_2]);
for k = 1:n_free_ener-1
    Efibmax(k) = max(Efib(:,k));
    Efibmin(k) = min(Efib(:,k));
    % Generate B-splines
    pfib(:,k)=Bspline_generator(control_pts(:,k+1),[Efibmin(k),Efibmax(k)],degree,Efib(:,k));
end
% generate bsplines
p1_1=Bspline_generator(control_pts(:,1),[I1min,I1max],degree,I1_1);
p1_2=Bspline_generator(control_pts(:,1),[I1min,I1max],degree,I1_2);
p4_1=Bspline_generator(control_pts(:,2),[E4min,E4max],degree,E4_1);
p4_2=Bspline_generator(control_pts(:,2),[E4min,E4max],degree,E4_2);

% Run the constitutive model
% [P11, ~] = Holzapfel_model_singlefamily(p1_1, p4_1, F, len);
% [~, P22] = Holzapfel_model_singlefamily(p1_2, p4_2, F, len);
switch dispmodel
    case "io"
%         [P11, ~] = Holzapfel_model(p1_1, p4_1, F, len);
%         [~, P22] = Holzapfel_model(p1_2, p4_2, F, len);
        [P11, ~] = io_stress(p1_1, pfib, F, len);
        [~, P22] = io_stress(p1_2, pfib, F, len);
    case "fs"
%         [sigma11, ~, ~, P11, ~, ~] = Eriksson13(p1_1, p4_1, F, len);
%         [~, sigma22, ~, ~, P22, ~] = Eriksson13(p1_2, p4_2, F, len);
        [sigma11, ~, ~, P11, ~, ~] = fs_stress(p1_1, pfib, F, len);
        [~, sigma22, ~, ~, P22, ~] = fs_stress(p1_2, pfib, F, len);
end
% Calculate mean square error
err_11 = 1/len * sum((P11-P_ut).^2);
err_22 = 1/len * sum((P22-P_ut).^2);

end
