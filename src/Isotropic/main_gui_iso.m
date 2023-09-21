%======================================================================%
%           Computational Micromechanics Laboratory (CMML)             %
%              Middle East Technical University (METU)                 %
%----------------------------------------------------------------------%
%                        Hyper-Data Toolbox - v1.0                     %
%----------------------------------------------------------------------%
%   Authors: Recep Durna, Alp Kağan Açan, Oğuz Ziya Tikenoğulları      %
%   Software License: GNU GPL v3.0   Copyright (C) 2023 Hüsnü Dal      %
%======================================================================%

% Size of data info array
sz_data = size(Data_Info);
% Call script to read data
% Initialize
data_pos = {};
lam_ut = [];
P_ut = [];
n_ut = 0;
lam_et = [];
P_et = [];
n_et = 0;
lam_ps = [];
P_ps = [];
n_ps = 0;
lam1_be = [];
P1_be = [];
lam2_be = [];
P2_be = [];
n_be = 0;

lam_ut_ext = [];
n_ut_ext = 0;
lam_et_ext = [];
n_et_ext = 0;
lam_ps_ext = [];
n_ps_ext = 0;
lam1_be_ext = [];
lam2_be_ext = [];
n_be_ext = 0;

for i = 1:sz_data(1);
    type = string(Data_Info(i,1));
    data_1 = dlmread(string(Data_Info(i,2)));
    data_pos(i,1) = {type};
    switch type
        case "UT"
            P_1       = data_1(:,2);
            lam_1      = data_1(:,1);  %stretch
            len         = length(lam_1);
            lam_ut = [lam_ut; lam_1];
            P_ut   = [P_ut; P_1];
            data_pos(i,2) = {[n_ut+1 n_ut+len]};
            n_ut   = length(lam_ut);

            if len < 100
                lam_1_ext = linspace(lam_1(1),lam_1(end),100)';
            else
                lam_1_ext = lam_1;
            end
            len_ext = length(lam_1_ext);
            lam_ut_ext = [lam_ut_ext; lam_1_ext];
            data_pos(i,3) = {[n_ut_ext+1 n_ut_ext+len_ext]};
            n_ut_ext   = length(lam_ut_ext);
        case "ET"
            P_1       = data_1(:,2);
            lam_1      = data_1(:,1);  %stretch
            len         = length(lam_1);
            lam_et = [lam_et; lam_1];
            P_et   = [P_et; P_1];
            data_pos(i,2) = {[n_et+1 n_et+len]};
            n_et   = length(lam_et);
            
            if len < 100
                lam_1_ext = linspace(lam_1(1),lam_1(end),100)';
            else
                lam_1_ext = lam_1;
            end
            len_ext = length(lam_1_ext);
            lam_et_ext = [lam_et_ext; lam_1_ext];
            data_pos(i,3) = {[n_et_ext+1 n_et_ext+len_ext]};
            n_et_ext   = length(lam_et_ext);
        case "PS"
            P_1       = data_1(:,2);
            lam_1      = data_1(:,1);  %stretch
            len         = length(lam_1);
            lam_ps = [lam_ps; lam_1];
            P_ps   = [P_ps; P_1];
            data_pos(i,2) = {[n_ps+1 n_ps+len]};
            n_ps   = length(lam_ps);

            if len < 100
                lam_1_ext = linspace(lam_1(1),lam_1(end),100)';
            else
                lam_1_ext = lam_1;
            end
            len_ext = length(lam_1_ext);
            lam_ps_ext = [lam_ps_ext; lam_1_ext];
            data_pos(i,3) = {[n_ps_ext+1 n_ps_ext+len_ext]};
            n_ps_ext   = length(lam_ps_ext);
        case "BT"
            P_1       = data_1(:,3);
            P_2       = data_1(:,4);
            lam_1      = data_1(:,1);  %stretch
            lam_2      = data_1(:,2);  %stretch
            len         = length(lam_1);
            lam1_be = [lam1_be; lam_1];
            lam2_be = [lam2_be; lam_2];
            P1_be   = [P1_be; P_1];
            P2_be   = [P2_be; P_2];
            data_pos(i,2) = {[n_be+1 n_be+len]};
            n_be   = length(lam1_be);

            if len < 100
                lam_1_ext = linspace(lam_1(1),lam_1(end),100)';
                lam_2_ext = linspace(lam_2(1),lam_2(end),100)';
            else
                lam_1_ext = lam_1;
                lam_2_ext = lam_2;
            end
            len_ext = length(lam_1_ext);
            lam1_be_ext = [lam1_be_ext; lam_1_ext];
            lam2_be_ext = [lam2_be_ext; lam_2_ext];
            data_pos(i,3) = {[n_be_ext+1 n_be_ext+len_ext]};
            n_be_ext   = length(lam1_be_ext);
    end
end

% Call respective function based on the formulation base
switch form_type
    case "Invariant-based"
        main_inv
    case "Modified invariant-based"
        main_modinv
    case "Principal stretch-based"
        main_str
end

% Write the identified parameters into a file
dataloc1 = string(Data_Info(1,2));
outdata

