%======================================================================%
%           Computational Micromechanics Laboratory (CMML)             %
%              Middle East Technical University (METU)                 %
%----------------------------------------------------------------------%
%                        Hyper-Data Toolbox - v1.0                     %
%----------------------------------------------------------------------%
%   Authors: Recep Durna, Alp Kağan Açan, Oğuz Ziya Tikenoğulları      %
%   Software License: GNU GPL v3.0   Copyright (C) 2023 Hüsnü Dal      %
%======================================================================%

% Define global variables
global dispmodel kip kop kappaf kappas fib kappafib n_free_ener

%%%%  Read experiment data  %%%% 
switch type
    % Equibiaxial Loading
    case "ET"
        ET_data     = dlmread(dataloc1);
        P_et2       = ET_data(:,3);
        P_et1       = ET_data(:,2);
        lam_et      = ET_data(:,1);  %stretch
        len         = length(lam_et);

        if len < 100
            lam_et_ext = linspace(lam_et(1),lam_et(end),100)';
        else
            lam_et_ext = lam_et;
        end
        len_ext = length(lam_et_ext);
    % Uniaxial Loading
    case "UT"
        P1_exp = dlmread(dataloc1);
        P2_exp = dlmread(dataloc2);
        %%%%        Uniaxial Loading in 1-direction       %%%%
        P_ut1   = P1_exp(:,2);
        lam_ut1 = P1_exp(:,1);
        n_ut1   = length(lam_ut1);
        %%%%        Uniaxial Loading in 2-direction       %%%%
        P_ut2   = P2_exp(:,2);
        lam_ut2 = P2_exp(:,1);
        n_ut2   = length(lam_ut2);

        if n_ut1 < 100
            lam_ut1_ext = linspace(lam_ut1(1),lam_ut1(end),100)';
        else
            lam_ut1_ext = lam_ut1;
        end
        n_ut1_ext = length(lam_ut1_ext);

        if n_ut2 < 100
            lam_ut2_ext = linspace(lam_ut2(1),lam_ut2(end),100)';
        else
            lam_ut2_ext = lam_ut2;
        end
        n_ut2_ext = length(lam_ut2_ext);
    % Triaxial Shear
    case "TS"
        fs_data     = dlmread(datalocfs);
        sf_data     = dlmread(datalocsf);
        fn_data     = dlmread(datalocfn);
        nf_data     = dlmread(datalocnf);
        sn_data     = dlmread(datalocsn);
        ns_data     = dlmread(datalocns);
        % experimental stress and shear values 
        shear   = fs_data(:,1); % the shear amount of each dataset is identical.
        Pexp_fs = fs_data(:,2); 
        Pexp_sf = sf_data(:,2);
        Pexp_fn = fn_data(:,2); 
        Pexp_nf = nf_data(:,2); 
        Pexp_sn = sn_data(:,2); 
        Pexp_ns = ns_data(:,2); 
        
        len     = length(shear);  % Length of the dataset

        if len < 100
            shear_ext = linspace(shear(1),shear(end),100)';
        else
            shear_ext = shear;
        end
        len_ext = length(shear_ext);
    % Biaxial Tension
    case "BT"
        P_exp1 = [];
        P_exp2 = [];
        lam_1 = [];
        lam_2 = [];
        lam_1_ext = [];
        lam_2_ext = [];
        bounds = [];
        bounds_ext = [];
        len = 0;
        len_ext = 0;
        for i = 1:length(Loc_Data)
            data_1     = dlmread(string(Loc_Data(i)));
            P_exp1       = [P_exp1; data_1(:,3)];
            P_exp2       = [P_exp2; data_1(:,4)];
            lam_1      = [lam_1; data_1(:,1)];  %stretch
            lam_2      = [lam_2; data_1(:,2)];  %stretch
            bounds(i,:) = [len+1 length(lam_1)];
            len = length(lam_1);

            if length(data_1(:,1)) < 100
                lam1_be_ext = linspace(data_1(1,1),data_1(end,1),100)';
                lam2_be_ext = linspace(data_1(1,2),data_1(end,2),100)';
            else
                lam1_be_ext = data_1(:,1);
                lam2_be_ext = data_1(:,2);
            end
            lam_1_ext = [lam_1_ext; lam1_be_ext];
            lam_2_ext = [lam_2_ext; lam2_be_ext];
            bounds_ext(i,:) = [len_ext+1 length(lam_1_ext)];
            len_ext = length(lam_1_ext);
        end
end

% Call respective function depending on the experiment type
switch type
    case "ET"
        main_et
    case "UT"
        main_ut
    case "TS"
        main_ts
    case "BT"
        main_bt
end

% Register identified parameters into a file
outdata

