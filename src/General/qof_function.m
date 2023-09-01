%======================================================================%
%           Computational Micromechanics Laboratory (CMML)             %
%              Middle East Technical University (METU)                 %
%----------------------------------------------------------------------%
%                        Hyper-Data Toolbox - v1.0                     %
%----------------------------------------------------------------------%
%   Authors: Recep Durna, Alp Kağan Açan, Oğuz Ziya Tikenoğulları      %
%   Software License: GNU GPL v3.0   Copyright (C) 2023 Hüsnü Dal      %
%======================================================================%

function qof_m = qof_function(lam,Pfit,Pexp)
    % This function is for quality of fit calculations
    n=length(lam);
    qof1 = 0;
    qof2 = 0;
    qof3 = 0;
    qof = 0;
    for i=1:round(n/3)
        if Pexp(i)~=0
            qof1 = qof1 + (Pfit(i)-Pexp(i))^2/Pexp(i);
        end
    end
    for i=(round(n/3)+1):round(2*n/3)
        if Pexp(i)~=0
            qof2 = qof2 + (Pfit(i)-Pexp(i))^2/Pexp(i);
        end
    end
    for i=(round(2*n/3)+1):n
        if Pexp(i)~=0
            qof3 = qof3 + (Pfit(i)-Pexp(i))^2/Pexp(i);
        end
    end
    for i=1:n
        if Pexp(i)~=0
            qof = qof + (Pfit(i)-Pexp(i))^2/Pexp(i);
        end
    end
    qof_m = [qof, qof1, qof2, qof3];

end