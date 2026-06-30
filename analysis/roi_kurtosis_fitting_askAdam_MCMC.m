% roi_kurtosis_fitting_askAdam_MCMC.m
% ----------------------------------------------------------------------
% Companion code for:
%   Sung D, et al. Microstructural and exchange imaging with oscillating
%   gradient spin-echo (OGSE) diffusion MRI. Magn Reson Med, 2025.
%   DOI: 10.1002/mrm.70300
%
% Author: Dongsuk Sung (dsung2@mgh.harvard.edu)
% Athinoula A. Martinos Center for Biomedical Imaging, MGH/Harvard Medical School
% Created: 2025-07-18
%
% Description:
%   Loads ROI masks (whole-brain, WM tracts, GM/cortical regions from
%   FreeSurfer/SynthSeg parcellations) for each subject, extracts
%   region-wise mean kurtosis across the OGSE frequency protocol, and
%   fits the exchange-time Karger model (tex, fn, aD, aC, Kinf) to the
%   frequency-dependent kurtosis using the askAdam solver (gacelle
%   toolbox) and/or MCMC.
%
% Dependencies (not included in this repo):
%   - gacelle toolbox   https://github.com/kschan0214/gacelle
%       (provides the askAdam minimisation solver, mcmc solver, and the
%        `utils` helper class used throughout this script)
%   - util_ogse         custom OGSE waveform utility class (see README.md)
%   - gpu_KurtModel.m   forward model class, included in ../models
%   - MATLAB Image Processing Toolbox (niftiread), Parallel Computing
%     Toolbox (parfor/gpuArray)
% ----------------------------------------------------------------------

%% Add path
clear; close all;
%% ==================== USER CONFIGURATION ====================
% Edit the paths below to match your local setup.
gacelle_dir   = '/path/to/gacelle';        % https://github.com/kschan0214/gacelle
util_ogse_dir = '/path/to/util_ogse';      % custom OGSE waveform utility class (see README)
project_dir   = '/path/to/OGSE_dataset';   % root of the OGSE BIDS dataset
% ==============================================================
addpath(genpath(gacelle_dir));
addpath(genpath(util_ogse_dir));
cd(project_dir);
addpath(genpath(project_dir));


%% Load data
do_plot = false;
for lr_combine = [0,1]
    subject_list = {'sub-ogse005-02','sub-ogse006','sub-ogse007','sub-ogse008','sub-ogse009','sub-ogse010',...
        'sub-ogse011','sub-ogse012','sub-ogse013','sub-ogse014','sub-ogse015'};
    for id = 1:numel(subject_list)
        subj_label = subject_list{id};
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf('Start Processing *** %s ****\n',subj_label);
        protocol = 'six_acqs'; % choose protocol: 1) all_acqs, 2) ten_acqs, 3) six_acqs, etc.
        switch protocol
            case 'all_acqs'

                acq={'30Hz_N1','40Hz_N1','50Hz_N1','60Hz_N1',...
                    '60Hz_N2','65Hz_N2','70Hz_N2','80Hz_N2','80Hz_N3','85Hz_N3','90Hz_N3'};
                N_list    = [1,1,1,1,2,2,2,2,3,3,3];
                freq_list = [30,40,50,60,60,65,70,80,80,85,90];                              % [Hz]
                Tm_list   = [7.4, 30.9, 24.7, 20.4, 16.5, 15.2, 14, 12.1, 11.5, 10.8, 10.1]; % Mixing Time [ms]
                list_fixN1 = find(N_list==1);
                list_fixN2 = find(N_list==2);
                list_fixN3 = find(N_list==3);
                list_fixT = [1,6,11];
                list_fixT2 = [1,7,9];
            case 'ten_acqs'

                acq={'40Hz_N1','50Hz_N1','60Hz_N1',...
                    '60Hz_N2','65Hz_N2','70Hz_N2','80Hz_N2','80Hz_N3','85Hz_N3','90Hz_N3'};
                N_list    = [1,1,1,2,2,2,2,3,3,3];
                freq_list = [40,50,60,60,65,70,80,80,85,90];                              % [Hz]
                Tm_list   = [30.9, 24.7, 20.4, 16.5, 15.2, 14, 12.1, 11.5, 10.8, 10.1]; % Mixing Time [ms]

            case 'six_acqs'
                acq={'30Hz_N1','40Hz_N1','50Hz_N1','60Hz_N1','65Hz_N2','90Hz_N3'};
                N_list    = [1,1,1,1,2,3];
                freq_list = [30,40,50,60,65,90];                 % [Hz]
                Tm_list   = [7.4, 30.9, 24.7, 20.4, 15.2, 10.1]; % Mixing Time [ms]
                list_fixN1 = [1,2,3,4];
                list_fixT = [1,5,6];
            case 'six_acqs2'
                acq={'30Hz_N1','40Hz_N1','50Hz_N1','60Hz_N1','60Hz_N2','80Hz_N3'};
                N_list    = [1,1,1,1,2,3];
                freq_list = [30,40,50,60,60,80];                              % [Hz]
                Tm_list   = [7.4, 30.9, 24.7, 20.4, 16.5, 11.5]; % Mixing Time [ms]
                list_fixN1 = [1,2,3,4];
                list_fixT = [1,5,6];
            case 'only_fixN'
                acq={'40Hz_N1','50Hz_N1','60Hz_N1'};
                N_list    = [1,1,1];
                freq_list = [40,50,60];                              % [Hz]
                Tm_list   = [30.9, 24.7, 20.4]; % Mixing Time [ms]
                list_fixN1 = [1,2,3];
                list_fixT = [1,2,3];
            case 'only_fixT'
                acq={'30Hz_N1','65Hz_N2','90Hz_N3'};
                N_list    = [1,2,3];
                freq_list = [30,65,90];                 % [Hz]
                Tm_list   = [7.4, 15.2, 10.1]; % Mixing Time [ms]
                list_fixN1 = [1,2,3];
                list_fixT = [1,2,3];
        end


        disp('Data loaded')
        %% Acquire ROI Masks

        dwi_proc_dirs = fullfile(project_dir,'bids/derivatives/mrtrix_tensor',subj_label);
        cd(dwi_proc_dirs);
        synthseg_dir=fullfile(project_dir,'bids/derivatives/synthseg',subj_label);
        mrtrix_dir = fullfile(project_dir,'bids/derivatives/mrtrix_tensor',subj_label);

        Nfreq = length(N_list);
        mk = cell(Nfreq,5);

        % Acquire ROI Masks
        brain_mask = double(niftiread(fullfile(project_dir,...
            'bids/derivatives/preprocessed',subj_label,sprintf('%s_dwi_brain_mask_gncorr_C2.nii.gz',subj_label))));
        wm_mask =  double(niftiread(fullfile(synthseg_dir,sprintf('%s_wm_mask_res.nii.gz',subj_label))));
        cortex_mask =  niftiread(fullfile(synthseg_dir,sprintf('%s_cortex_mask_res.nii.gz',subj_label)));
        csf_mask =  niftiread(fullfile(synthseg_dir,sprintf('%s_csf_mask_res.nii.gz',subj_label)));
        subnuclei_mask =  niftiread(fullfile(synthseg_dir,sprintf('%s_subnuclei_mask_res.nii.gz',subj_label)));
        gm_mask = cortex_mask+subnuclei_mask;
        gm_mask(gm_mask(:)>1) = 1;

        % WM regions
        template = round(niftiread(fullfile(mrtrix_dir,sprintf('%s_acq-90Hz_N3_dwi_JHU-labels-2mm.nii',subj_label))));
        wm_label = unique(nonzeros(template(:)));
        wm_region_mask = cell((numel(wm_label)+6)/2,1);
        if lr_combine == 1
            for wi = 1:(numel(wm_label)+6)/2
                if sum(ismember(1:6,wm_label(wi)))
                    wm_region_mask{wi} = double(template==wm_label(wi));
                else
                    wm_region_mask{wi} = double(template==wm_label(2*wi-7)) + double(template==wm_label(2*wi-6));
                end
            end
        else
            for wi = 1:numel(wm_label)
                wm_region_mask{wi} = double(template==wm_label(wi));
            end
        end

        % Larger WM regions (Tract-based)
        tract_label = {[31:38, 41:46], [6:10, 15:30], [3,4,5]};
        tract_mask = cell(3,1);
        for wi = 1:numel(wm_label)
            for ti = 1:length(tract_mask)
                if sum(ismember(tract_label{ti},wm_label(wi)))==1
                    if isempty(tract_mask{ti})
                        tract_mask{ti} = double(template==wm_label(wi));
                    else
                        tract_mask{ti} = tract_mask{ti} + double(template==wm_label(wi));
                    end
                end
            end
        end

        % GM Regions
        gm_label = [10:13,17,18,26,1001:1003,1005:1035, 49:54,58,2001:2003,2005:2035];%[12,51,1024,2024,1022,2022];
        gm_group_label = {10,11,12,13,17,18,26, ...
            1028, [1003,1027], [1018,1019,1020], [1012,1014], 1024, 1017, 1032, ...
            1029,1008,1031,1022,1025,[1009,1015,1030],1001,1007,1034,1006,1033,1016,...
            1011,1013,1005,1021,1026,1002,1023,1010,1035,...
            49,50,51,52,53,54,58, ...
            2028, [2003,2027], [2018,2019,2020], [2012,2014], 2024, 2017, 2032, ...
            2029,2008,2031,2022,2025,[2009,2015,2030],2001,2007,2034,2007,2033,2016,...
            2011,2013,2005,2021,2026,2002,2023,2010,2035};
        lobe_label = {[10:13,17,18,26], [1003,1012,1014,1017:1020,1024,1027,1028,1032],...
            [1008,1022,1025,1029,1031], [1001,1006,1007,1009,1015,1016,1030,1033,1034],...
            [1005,1011,1013,1021], [1002,1010,1023,1026],...
            [49:54,58], [2003,2012,2014,2017:2020,2024,2027,2028,2032],...
            [2008,2022,2025,2029,2031], [2001,2006,2007,2009,2015,2016,2030,2033,2034],...
            [2005,2011,2013,2021], [2002,2010,2023,2026]};



        if lr_combine == 1
            gm_region_mask = cell(numel(gm_label)/2,1);
            gm_group_mask = cell(length(gm_group_label)/2,1);
            lobe_mask = cell(length(lobe_label)/2,1);
            for gi = 1:numel(gm_region_mask)
                gm_region_mask_L = double(niftiread(fullfile(synthseg_dir,'all_labels',...
                    sprintf('%s_synthseg_%d_res.nii.gz',subj_label,gm_label(gi)))));
                gm_region_mask_R = double(niftiread(fullfile(synthseg_dir,'all_labels',...
                    sprintf('%s_synthseg_%d_res.nii.gz',subj_label,gm_label(numel(gm_label)/2+gi)))));
                gm_region_mask{gi} = gm_region_mask_L + gm_region_mask_R;

                for li = 1:length(gm_group_mask)
                    if sum(ismember(gm_group_label{li},gm_label(gi)))==1
                        if isempty(gm_group_mask{li})
                            gm_group_mask{li} = gm_region_mask_L + gm_region_mask_R;
                        else
                            gm_group_mask{li} = gm_group_mask{li} + gm_region_mask_L + gm_region_mask_R;
                        end
                    end
                end

                for li = 1:length(lobe_mask)
                    if sum(ismember(lobe_label{li},gm_label(gi)))==1
                        if isempty(lobe_mask{li})
                            lobe_mask{li} = gm_region_mask_L + gm_region_mask_R;
                        else
                            lobe_mask{li} = lobe_mask{li} + gm_region_mask_L + gm_region_mask_R;
                        end
                    end
                end
            end
        else
            gm_region_mask = cell(numel(gm_label),1);
            gm_group_mask = cell(length(gm_group_label),1);
            lobe_mask = cell(length(lobe_label),1);
            for gi = 1:numel(gm_label)
                gm_region_mask{gi} = double(niftiread(fullfile(synthseg_dir,'all_labels',...
                    sprintf('%s_synthseg_%d_res.nii.gz',subj_label,gm_label(gi)))));

                for li = 1:length(gm_group_mask)
                    if sum(ismember(gm_group_label{li},gm_label(gi)))==1
                        if isempty(gm_group_mask{li})
                            gm_group_mask{li} = gm_region_mask{gi};
                        else
                            gm_group_mask{li} = gm_group_mask{li} + gm_region_mask{gi};
                        end
                    end
                end

                for li = 1:length(lobe_mask)
                    if sum(ismember(lobe_label{li},gm_label(gi)))==1
                        if isempty(lobe_mask{li})
                            lobe_mask{li} = gm_region_mask{gi};
                        else
                            lobe_mask{li} = lobe_mask{li} + gm_region_mask{gi};
                        end
                    end
                end
            end
        end
        disp('All ROIs were acquired')
        %% Acquire median MK across each ROI
        if lr_combine==1
            Nlabel = 1+(numel(wm_label)+6)/2+numel(tract_mask) + 1+numel(gm_group_mask)+numel(lobe_mask);
        else
            Nlabel = 1+numel(wm_label)+numel(tract_mask) + 1+numel(gm_group_mask)+numel(lobe_mask);
        end
        md_roi = cell(length(acq),Nlabel);
        mk_roi = cell(length(acq),Nlabel);

        mean_md = zeros(Nfreq,4);
        std_md = zeros(Nfreq,4);
        median_md = zeros(Nfreq,4);
        iqr_md = zeros(Nfreq,4);

        mean_mk = zeros(Nfreq,4);
        std_mk = zeros(Nfreq,4);
        median_mk = zeros(Nfreq,4);
        iqr_mk = zeros(Nfreq,4);

        for i = 1:length(acq)
            %%% Load MD & MK %%%
            md = 1000*niftiread(fullfile(dwi_proc_dirs,sprintf('%s_acq-%s_dwi_adc.nii.gz',subj_label,acq{i})));
            mk = niftiread(fullfile(dwi_proc_dirs,sprintf('%s_acq-%s_dwi_mk.nii.gz',subj_label,acq{i})));

            %%% Remove outliers %%%
            lower_lim = mk<0;
            upper_lim = mk>3;
            mk(lower_lim) = NaN; %0;
            mk(upper_lim) = NaN; %3

            %%% Apply masks %%%
            % WM
            md_roi{i,1} = md.*wm_mask;
            mk_roi{i,1} = mk.*wm_mask;

            % WM Regions
            cnt = 2;
            if lr_combine==1
                nwm_label = (numel(wm_label)+6)/2;
            else
                nwm_label = numel(wm_label);
            end
            for wi = 1:nwm_label
                md_roi{i,cnt} = md.*wm_region_mask{wi};
                mk_roi{i,cnt} = mk.*wm_region_mask{wi};
                cnt = cnt + 1;
            end

            % Tract Regions
            for ti = 1:length(tract_mask)
                md_roi{i,cnt} = md.*tract_mask{ti};
                mk_roi{i,cnt} = mk.*tract_mask{ti};
                cnt = cnt + 1;
            end

            % Cortex
            md_roi{i,cnt} = md.*cortex_mask;
            mk_roi{i,cnt} = mk.*cortex_mask;
            cnt = cnt+1;

            % GM Regions in Group
            for hi = 1:(2-lr_combine)
                for gi = 1:length(gm_group_mask)/(2-lr_combine)
                    half_Nsample = (hi-1)*length(gm_group_mask)/2;
                    md_roi{i,cnt} = md.*gm_group_mask{half_Nsample+gi};
                    mk_roi{i,cnt} = mk.*gm_group_mask{half_Nsample+gi};
                    cnt = cnt + 1;
                end

                % Lobes
                for li = 1:length(lobe_mask)/(2-lr_combine)
                    half_Nsample = (hi-1)*length(lobe_mask)/2;
                    md_roi{i,cnt} = md.*lobe_mask{half_Nsample+li};
                    mk_roi{i,cnt} = mk.*lobe_mask{half_Nsample+li};
                    cnt = cnt + 1;
                end
            end

            %%% Stats %%%
            for j = 1:cnt-1
                % MD
                mean_md(i,j) = mean(nonzeros(md_roi{i,j}(:)),'omitnan');
                std_md(i,j) = std(nonzeros(md_roi{i,j}(:)),'omitnan');
                median_md(i,j) = median(nonzeros(md_roi{i,j}(:)),'omitnan');
                iqr_md(i,j) = iqr(nonzeros(md_roi{i,j}(:)));

                % MK
                mean_mk(i,j) = mean(nonzeros(mk_roi{i,j}(:)),'omitnan');
                std_mk(i,j) = std(nonzeros(mk_roi{i,j}(:)),'omitnan');
                median_mk(i,j) = median(nonzeros(mk_roi{i,j}(:)),'omitnan');
                iqr_mk(i,j) = iqr(nonzeros(mk_roi{i,j}(:)));
            end
        end

        disp('Done load data')

        %% Prepare input parameters
        % Calculate waveform related parameters: t, g(t), q(t)^2, and q(w)^2
        Nt = 1e3;
        tr = 1.59;
        ut = util_ogse();
        t_vec = zeros(numel(N_list),Nt,'single');
        gt_vec = zeros(numel(N_list),Nt,'single');
        qt2_vec = zeros(numel(N_list),Nt,'single');
        qt2_fun = cell(numel(N_list),1);
        qw2_int = zeros(numel(N_list),1,'single');
        true_freq_list = zeros(numel(N_list),1,'single');

        bval = 2;
        tau = 1000./freq_list;
        Gmax = (sqrt(bval)*40/0.0107) ./ sqrt( N_list.*((1/24)*tau.^3 - (2/3)*tr^2*tau + (16/15)*tr^3) + (1/30)*tr^3 );
        for ii = 1:numel(N_list)
            [t_vec(ii,:),gt_vec(ii,:)] = ut.costrap_waveform(N_list(ii),freq_list(ii),Tm_list(ii),Gmax(ii),tr,Nt); % t[ms]
            dt = t_vec(ii,2)-t_vec(ii,1);
            qt2_vec(ii,:) = (cumsum(gt_vec(ii,:))*dt).^2;
            qt2_fun{ii} = griddedInterpolant(t_vec(ii,:), qt2_vec(ii,:), 'linear', 'none');
            [w_tmp,qw2_tmp,bval_tmp] = ut.convQw(t_vec(ii,:), gt_vec(ii,:));
            qw2_int(ii) = 1./(2*pi*bval_tmp).*trapz(w_tmp, sqrt(abs(w_tmp)).*qw2_tmp);
            [~,loc] = findpeaks(qw2_tmp,w_tmp,'SortStr','descend');
            true_freq_list(ii) = 1000*abs(loc(1))/(2*pi);
        end

        disp('Data is prepared')

        %% Choose solver
        % close all
        model_obj = gpu_KurtModel(t_vec, gt_vec, qw2_int);
        ROI_list = {'WM_ROI','GM_ROI'};
        for roii = 1:numel(ROI_list)

            ROI = ROI_list{roii};
            lmax = 0;
            if lmax==0; Nparam=5; else; Nparam=6; end

            % Mask for input
            if lr_combine==1
                nWM = 1+(numel(wm_label)+6)/2+numel(tract_mask);
                nGM = 1+numel(gm_group_mask)+numel(lobe_mask);
            else
                nWM = 1+numel(wm_label)+numel(tract_mask);
                nGM = 1+numel(gm_group_mask)+numel(lobe_mask);
            end
            switch ROI
                case 'WM_ROI'
                    mask_sample = ones(1,nWM)>0;
                    startPoint  = [150; 0.7; 0.4; 5; 0.5];
                case 'GM_ROI'
                    mask_sample = ones(1,nGM)>0;
                    startPoint  = [20; 0.4;  0.5; 10;  0.5];
            end
            Nsample = sum(mask_sample(:));
            ind = find(mask_sample(:));

            solver = 'askAdam'; % choose between two options: 1) askAdam, 2) MCMC

            switch solver
                case 'askAdam'
                    %% Model fitting using askAdam
                    disp('Start askAdam fitting')
                    start_entire = tic;
                    % Fitting parameters
                    fitting = [];
                    fitting.iteration           = 10000;
                    fitting.initialLearnRate    = 0.001;
                    fitting.decayRate           = 0.0001;
                    fitting.convergenceValue    = 1e-10;         % NEXI: 1e-8
                    fitting.lossFunction        = 'l1';         % 'l1', 'l2', 'huber', 'mse'
                    fitting.tol                 = 1e-4;         % NEXI: 1e-4
                    fitting.patience            = 5;          % NEXI = 5
                    fitting.isDisplay           = false;
                    fitting.lmax                = 0;
                    fitting.start               = 'default';

                    %%%% fitting parameters for TV regularization %%%%
                    % fitting.lambda              = 0.001;
                    % fitting.TVmode              = '2D';
                    % fitting.voxelSize           = [2,2,2];


                    disp(subj_label);
                    % get all fitting algorithm parameters
                    fitting                 = model_obj.check_set_default(fitting);
                    % determine fitting parameters
                    model_obj               = model_obj.updateProperty(fitting);
                    fitting.modelParams     = model_obj.modelParams;
                    % set fitting boundary if no input from user
                    if isempty( fitting.ub); fitting.ub = model_obj.ub(1:numel(fitting.modelParams)); end
                    if isempty( fitting.lb); fitting.lb = model_obj.lb(1:numel(fitting.modelParams)); end
                    % setup fitting weights
                    w = model_obj.compute_optimisation_weights(mask_sample,fitting.lossFunction,fitting.lmax); % This is a customised funtion
                    w = squeeze(w)';
                    % askAdam optimisation main
                    askadamObj = askadam();

                    disp('MK fitting')
                    modelFWD = @model_obj.FWD;
                    switch ROI
                        case 'WM_ROI'
                            y = median_mk(:,1:nWM);
                        case 'GM_ROI'
                            y = median_mk(:,(nWM+1):end);
                    end

                    fprintf('Generate model-based data N=%d, Nacq=%s, Nparam=%d, lmax=%d, loss=%s\n',Nsample,protocol,Nparam,lmax,fitting.lossFunction)
                    fprintf('Estimation lower bound [%s]: [%s]\n',cell2str(model_obj.modelParams),num2str(model_obj.lb(:).',' %6.1f'));
                    fprintf('Estimation upper bound [%s]: [%s]\n',cell2str(model_obj.modelParams),num2str(model_obj.ub(:).',' %5.1f'));
                    ('---------------');
                    start = tic;
                    fprintf('Using default starting points for all voxels at [%s]: [%s]\n', cell2str(model_obj.modelParams),replace(num2str(model_obj.startPoint(:).',' %.2f'),' ',','));
                    pars0 = zeros(Nparam,size(y,2));

                    switch fitting.start
                        case 'likelihood'
                            disp('Estimate starting points based on likelihood ...')
                            % Generate ground truth parameter maps
                            Nsamples = 1e6;
                            [x_train, S_train] = model_obj.traindata(modelFWD, Nsamples, lmax);
                            pars0 = model_obj.likelihood(y, x_train, S_train, lmax);
                            fprintf('Likelihood method was used with number of sample = %d\n',Nsamples);
                        case 'default'
                            % use fixed points
                            fprintf('Using default starting points for all voxels at [%s]: [%s]\n', cell2str(model_obj.modelParams),replace(num2str(model_obj.startPoint(:).',' %.2f'),' ',','));
                            pars0 = zeros(Nparam,size(y,2));
                            for k = 1:Nparam
                                pars0(k,:) = ones([1,size(y,2)],'single') *startPoint(k);
                            end
                    end

                    ET  = duration(0,0,toc(start),'Format','hh:mm:ss');
                    fprintf('Starting points estimated. Elapsed time (hh:mm:ss): %s \n',string(ET));

                    % initiate starting points arrays
                    ni = 5000;
                    NSegment = ceil(Nsample/ni);
                    out = cell(NSegment,1);
                    pars_i = cell(NSegment,1);
                    mask_sample_comb = mask_sample(mask_sample(:)==1);
                    sum_vox = 0;
                    for ii = 1:NSegment
                        if NSegment ~= 1
                            fprintf('Running #Segment = %d/%d \n',ii,NSegment);
                            disp   ('------------------------')
                        end
                        if ii==NSegment
                            sample_rng = ((ii-1)*ni+1):Nsample;
                        else
                            sample_rng = ((ii-1)*ni+1):(ii*ni);
                        end

                        for km = 1:Nparam
                            pars_i{ii}.(model_obj.modelParams{km}) = pars0(km,sample_rng);
                        end
                        yi = y(:,sample_rng);
                        mask_i = mask_sample_comb(sample_rng);
                        w_i = w(:,sample_rng);

                        sum_vox = sum_vox + sum(mask_i(:));
                        fprintf('Running #voxel = %d/%d \n',sum_vox,sum(mask_sample(:)));

                        out{ii} = askadamObj.optimisation(yi, mask_i, w_i, pars_i{ii}, ...
                            fitting, @model_obj.FWD, fitting.lmax);
                    end

                    out_final = out{1}.final;
                    for ii = 2:ceil(Nsample/ni)
                        for pp = 1:Nparam
                            out_final.(modelParams{pp}) = cat(1,out_final.(modelParams{pp}),out{ii}.final.(modelParams{pp}));
                        end
                    end

                    % Fitting error
                    MK_fit = modelFWD(out_final);
                    MK_rmse = sqrt(mean((y-MK_fit).^2));
                    MK_sse = sum((y- MK_fit).^2);
                    ET  = duration(0,0,toc(start_entire),'Format','hh:mm:ss');
                    fprintf('Fitting is Done. Elapsed time for entire fitting(hh:mm:ss): %s \n',string(ET));

                case 'MCMC'
                    %% MCMC fitting
                    disp('Start MCMC fitting')
                    start_entire = tic;
                    model_obj = gpu_KurtModel(t_vec, gt_vec, qw2_int);
                    % Input Data
                    input = [];
                    input.Nt = Nt;
                    input.tr = 1.59;
                    input.t_vec = t_vec;
                    input.gt_vec = gt_vec;
                    input.qw2_int = qw2_int;
                    input.modelParams = model_obj.modelParams;
                    input.intervals = [model_obj.lb, model_obj.ub];

                    % Set up fitting algorithm
                    fitting                     = [];
                    % define model parameter name and fitting boundary
                    fitting.algorithm           = 'GW';
                    fitting.Nwalker             = 50;
                    fitting.modelParams         = model_obj.modelParams;
                    fitting.lb                  = model_obj.lb;    % lower bound
                    fitting.ub                  = model_obj.ub;     % upper bound
                    % Estimation algorithm setting
                    fitting.start               = 'default';
                    fitting.iteration           = 10e5/fitting.Nwalker;
                    fitting.burnin              = 0.1;     % 10% iterations
                    fitting.thinning            = 100;
                    fitting.StepSize            = 2;
                    fitting.diffusion           = 0;
                    fitting.startRadius         = 1e-3;
                    fitting.metric              = {'median','mode'};
                    fitting.Nbin                = 1e3;

                    mcmc_obj    = mcmc_exvivo_DS;
                    modelParams = model_obj.modelParams;
                    out = struct;
                    for i = 1:length(modelParams)
                        xPosterior.(modelParams{i}) = [];
                    end

                    disp('MK fitting')
                    modelFWD    = @ogse_MK_FWD;
                    % modelFWD = @model_obj.FWD;
                    switch ROI
                        case 'WM_ROI'
                            y = median_mk(:,1:nWM);
                        case 'GM_ROI'
                            y = median_mk(:,(nWM+1):end);
                    end

                    % initialization based on a dictionary
                    start = tic;
                    switch fitting.start
                        case 'likelihood'
                            disp('Estimate starting points based on likelihood ...')
                            % Generate ground truth parameter maps
                            Nsamples = 1e6;
                            [x_train, S_train] = model_obj.traindata(modelFWD, Nsamples, lmax, input);
                            pars0 = model_obj.likelihood(y, x_train, S_train, lmax);
                            fprintf('Likelihood method was used with number of sample = %d\n',Nsamples);
                        case 'default'
                            % use fixed points
                            fprintf('Using default starting points for all voxels at [%s]: [%s]\n', cell2str(model_obj.modelParams),replace(num2str(model_obj.startPoint(:).',' %.2f'),' ',','));
                            pars0 = zeros(Nparam,size(y,2));
                            for k = 1:Nparam
                                pars0(k,:) = ones([1,size(y,2)],'single') *model_obj.startPoint(k);
                            end
                    end
                    if fitting.diffusion==0
                        x0 = struct;
                        for i = 1:length(modelParams)-1
                            x0.(modelParams{i}) = pars0(i,:);
                        end
                    else
                        % x0.(modelParams{1}) = single(pars0(1,:));
                        x0.(modelParams{1})  = ones([1,size(y,2)],'single') *model_obj.startPoint(1);
                        for i = 2:length(modelParams)
                            x0.(modelParams{i}) = out_md.(modelParams{i});
                        end
                    end
                    x0.noise = single(ones(1,size(y,2))*model_obj.startPoint(end));
                    ET  = duration(0,0,toc(start),'Format','hh:mm:ss');
                    fprintf('Starting points estimated. Elapsed time (hh:mm:ss): %s \n',string(ET));

                    % decide weights
                    weight_type = 'equal';
                    switch weight_type
                        case 'equal'
                            weights = [];
                        case 'L1'
                            % weight for 'l1' lossFunction is sqrt(Nav/(2l+1))
                            weights = [ones(size(y(1:floor(end/2),:))); 1/sqrt(5)*ones(size(y(floor(end/2)+1:end,:)))];
                        case 'L2'
                            % weight for 'l2' lossFunction is Nav/(2l+1)
                            weights = [ones(size(y(1:floor(end/2),:))); 1/5*ones(size(y(floor(end/2)+1:end,:)))];
                    end

                    % MCMC
                    if strcmpi(fitting.algorithm,'MH')
                        out_t = mcmc_obj.metropolis_hastings(y, x0, weights, fitting, modelFWD ,varargin{:});
                    else
                        out_t = mcmc_obj.goodman_weare(y, x0, weights, fitting, modelFWD, input);
                    end

                    for i = 1:length(modelParams)
                        xPosterior.(modelParams{i}) = [];
                    end
                    for i = 1:length(modelParams)
                        out_t.reshape.(modelParams{i}) = reshape(out_t.(modelParams{i}), ...
                            [size(y,2), prod(size(out_t.(modelParams{i}),2:3))]);
                        out_t.median.(modelParams{i}) = median(out_t.reshape.(modelParams{i}), 2);
                        xPosterior.(modelParams{i}) = cat(1,xPosterior.(modelParams{i}),out_t.reshape.(modelParams{i}));
                    end

                    % generate 3D
                    for i = 1:length(modelParams)
                        out.mean.(modelParams{i}) = mean(xPosterior.(modelParams{i}),2)';
                        out.std.(modelParams{i}) = std(xPosterior.(modelParams{i}),[],2)';
                        out.median.(modelParams{i}) = median(xPosterior.(modelParams{i}),2)';
                        out.mode.(modelParams{i}) = mode(xPosterior.(modelParams{i}),2)';
                    end

                    out_final = out.median;
                    % Fitting error
                    MK_fit = modelFWD(out_final,input);
                    MK_rmse = sqrt(mean((y-MK_fit).^2));
                    MK_sse = sum((y- MK_fit).^2);
                    ET  = duration(0,0,toc(start_entire),'Format','hh:mm:ss');
                    fprintf('Fitting is Done. Elapsed time for entire fitting(hh:mm:ss): %s \n',string(ET));

            end

            %% Save fitting results
            model_label = 'model-KurtFit';
            acq_label = protocol;
            solver_label= sprintf('solver-%s',solver);
            if strcmp(solver,'askAdam')
                loss_label  = fitting.lossFunction;
            else
                loss_label  = fitting.iteration;
            end
            initLR   = '1e-3';
            startmethod = fitting.start;
            if lr_combine==1
                LRcomb = 'combine2';
            else
                LRcomb = 'separate2'; % 'separate' or 'combine'
            end

            output_prefix = strcat(subj_label,'_',ROI,'_',model_label,'_',acq_label,'_',...
                solver_label,'_',loss_label,'_iLR',initLR,'_',startmethod,'_LR',LRcomb,'_exclude');


            save_dir = fullfile(project_dir,sprintf('bids/derivatives/%s/%s',model_label,subj_label));
            if ~isfolder(save_dir)
                mkdir(save_dir)
            end

            save(fullfile(save_dir,sprintf('out_%s.mat',output_prefix)),'out_final')
            fprintf('%s %s Save Complete\n',subj_label,ROI)

            %% Plot fitting result
            if do_plot
                if lr_combine==1
                    switch ROI
                        case 'WM_ROI'
                            RegionLabel = {'WM','CC','PLIC','SCR','SLF'};
                            chosen_roi = [1,31,14,18,25];
                        case 'GM_ROI'
                            RegionLabel = {'Cortex','Frontal','Parietal','Temporal','Occipital'};
                            chosen_roi = [32,69,70,71,72]-nWM;
                    end
                else
                    switch ROI
                        case 'WM_ROI'
                            RegionLabel = {'WM','R-PLIC','R-SCR','R-SLF','CC','L-PLIC','L-SCR','L-SLF'};
                            chosen_roi = [1,20,26,42,52,21,27,43];
                        case 'GM_ROI'
                            RegionLabel = {'Cortex','L-Frontal','L-Parietal','L-Temporal','L-Occipital',...
                                'R-Frontal','R-Parietal','R-Temporal','R-Occipital'};
                            chosen_roi = [1,38,39,40,41,79,80,81,82];
                    end
                end

                % Input parameters for fixed N
                bval=2;
                Freq_Plot = 1:100;
                tau=1000./Freq_Plot;
                N_Plot_fixN1 = 1*ones(size(Freq_Plot));
                N_Plot_fixN2 = 2*ones(size(Freq_Plot));
                N_Plot_fixN3 = 3*ones(size(Freq_Plot));
                Gmax_Plot_fixN1 = (sqrt(bval)*40/0.0107) ./ sqrt( N_Plot_fixN1.*((1/24)*tau.^3 - (2/3)*tr^2*tau + (16/15)*tr^3) + (1/30)*tr^3 );
                Gmax_Plot_fixN2 = (sqrt(bval)*40/0.0107) ./ sqrt( N_Plot_fixN2.*((1/24)*tau.^3 - (2/3)*tr^2*tau + (16/15)*tr^3) + (1/30)*tr^3 );
                Gmax_Plot_fixN3 = (sqrt(bval)*40/0.0107) ./ sqrt( N_Plot_fixN3.*((1/24)*tau.^3 - (2/3)*tr^2*tau + (16/15)*tr^3) + (1/30)*tr^3 );

                % Input parameters for fixed T
                Freq_Plot_fixT1 = [40,60,85];
                tau_fixT1 = 1000./Freq_Plot_fixT1;
                N_Plot_fixT1 = [1,2,3];
                Gmax_Plot_fixT1 = (sqrt(bval)*40/0.0107) ...
                    ./ sqrt( N_Plot_fixT1.*((1/24)*tau_fixT1.^3 - (2/3)*tr^2*tau_fixT1 + (16/15)*tr^3) + (1/30)*tr^3 );
                Tm_Plot_fixT1 = [30.9, 16.5, 10.8]';

                Freq_Plot_fixT2 = [30,65,90];
                tau_fixT2 = 1000./Freq_Plot_fixT2;
                N_Plot_fixT2 = [1,2,3];
                Gmax_Plot_fixT2 = (sqrt(bval)*40/0.0107) ...
                    ./ sqrt( N_Plot_fixT2.*((1/24)*tau_fixT2.^3 - (2/3)*tr^2*tau_fixT2 + (16/15)*tr^3) + (1/30)*tr^3 );
                Tm_Plot_fixT2 = [7.4, 15.2, 10.1]';

                % Find optimal Tm
                Nt = 1e4;
                Tm_Opt1 = zeros(numel(Freq_Plot),3);
                Tm_Opt2 = zeros(numel(Freq_Plot),3);
                f0_list = zeros(numel(Freq_Plot),3);
                for N = 1:3
                    for i = 1:length(Freq_Plot)
                        freq = Freq_Plot(i);
                        tau = 1000/freq;
                        Gmax = (sqrt(bval)*40/0.0107) ./ ...
                            sqrt( N.*((1/24)*tau.^3 - (2/3)*tr^2*tau + (16/15)*tr^3) + (1/30)*tr^3 );
                        % Create cosine trapazoidal waveform
                        [t1,gt1] = ut.costrap_waveform_half(N,freq,Gmax,tr,Nt);
                        [w_tmp,qw2_tmp,bval_tmp] = ut.convQw(t1', gt1');
                        qw2_int_tmp = 1./(2*pi*bval_tmp).*trapz(w_tmp, sqrt(abs(w_tmp)).*qw2_tmp);
                        [~,loc] = findpeaks(qw2_tmp,w_tmp,'SortStr','descend');
                        f0_list(i,N) = 1000*abs(loc(1))/(2*pi);
                        Tm_Opt1(i,N) = 1000*(N)/(f0_list(i,N)) - 1000*N/freq - tr - 1/N; % 1/N was empircally subtracted.
                        Tm_Opt2(i,N) = 1000*(N+1)/(f0_list(i,N)) - 1000*N/freq - tr - 1/N;
                    end
                end
                Tm_Opt1(Tm_Opt1<7.4) = 10000;
                Tm_Opt = min(Tm_Opt1,Tm_Opt2);
                Tm_Opt(Tm_Opt>40) = 7.4;

                K_fit_fixN = zeros(numel(Freq_Plot),3,size(y,2));
                for N = 1:3
                    for i = 1:length(Freq_Plot)
                        freq = Freq_Plot(i);
                        tau = 1000/freq;
                        Gmax = (sqrt(bval)*40/0.0107) ./ ...
                            sqrt( N.*((1/24)*tau.^3 - (2/3)*tr^2*tau + (16/15)*tr^3) + (1/30)*tr^3 );
                        % Create cosine trapazoidal waveform
                        [tvec,gtvec] = ut.costrap_waveform(N,freq,Tm_Opt(i,N),Gmax,tr,Nt); % t[ms]
                        dt = tvec(2)-tvec(1);
                        [w_tmp,qw2_tmp,bval_tmp] = ut.convQw(tvec, gtvec);
                        qw2int = 1./(2*pi*bval_tmp).*trapz(w_tmp, sqrt(abs(w_tmp)).*qw2_tmp');

                        if strcmp(solver,'askAdam')
                            model_obj_tmp = gpu_KurtModel(tvec', gtvec', qw2int);
                            K_fit_fixN(i,N,:) = model_obj_tmp.FWD(out_final,lmax);
                        else
                            input.t_vec = tvec';
                            input.gt_vec = gtvec';
                            input.qw2_int = qw2int;
                            K_fit_fixN(i,N,:) = ogse_MK_FWD(out_final,input);
                        end

                    end
                end

                figureHandle = figure();

                for ri = 1:numel(chosen_roi)
                    MK_Reg = y(:,chosen_roi(ri));
                    tex = out_final.tex(chosen_roi(ri));
                    fn = out_final.fn(chosen_roi(ri));
                    aD = out_final.aD(chosen_roi(ri));
                    aC = out_final.aC(chosen_roi(ri));
                    Kinf = out_final.Kinf(chosen_roi(ri));

                    disp('%%%%%%%%%%%%%%%%%%%%%%%')
                    fprintf('Region: %s\n',RegionLabel{ri});
                    fprintf('t_{ex} = %.4f\n', tex);
                    fprintf('f = %.4f\n', fn);
                    fprintf('aD = %.4f\n', aD);
                    fprintf('aC = %.4f\n', aC);
                    fprintf('Kinf = %.4f\n', Kinf);
                    fprintf('RMSE = %f\n',MK_rmse(chosen_roi(ri)));


                    K_fit_fixN1 = K_fit_fixN(:,1,chosen_roi(ri));
                    K_fit_fixN2 = K_fit_fixN(:,2,chosen_roi(ri));
                    K_fit_fixN3 = K_fit_fixN(:,3,chosen_roi(ri));
                    K_fit_fixT1 = MK_fit(list_fixT,chosen_roi(ri));
                    K_fit_fixT2 = MK_fit(list_fixT,chosen_roi(ri));

                    if lr_combine==0
                        subplot(2,round(length(RegionLabel)/2),ri)
                    else
                        subplot(1,round(length(RegionLabel)),ri)
                    end
                    cmap = colormap('lines');
                    hold on;
                    h_fixT_y = plot(true_freq_list(list_fixT), MK_Reg(list_fixT), 'v', 'Color', cmap(1,:), 'markersize', 18, 'linewidth', 2);
                    h_fixN_y1 = plot(true_freq_list(list_fixN1), MK_Reg(list_fixN1), '.', 'Color', cmap(2,:), 'markersize', 40);
                    h_fixT_fit2 = plot([f0_list(30,1),f0_list(65,2),f0_list(90,3)], K_fit_fixT2, '-', 'Color', cmap(1,:), 'linewidth', 2);
                    h_fixN_fit1 = plot(f0_list(:,1), K_fit_fixN1, '-', 'Color', cmap(2,:), 'linewidth', 2);

                    switch ROI
                        case 'WM_ROI'
                            ylim([0.6, 1]);
                        case 'GM_ROI'
                            ylim([0.4, 0.7]);
                    end
                    xlabel('Frequency [Hz]')
                    ylabel('Mean Kurtosis')
                    title(sprintf('%s',RegionLabel{ri}));
                    set(gca,'FontSize',20)
                end

                fprintf('\nMean RMSE = %f\n', mean(MK_rmse,'omitnan'));

                % Save figures.
                % Get screen size
                screenSize = get(0, 'ScreenSize');  % [left bottom width height]

                % Set figure position:
                % [left, bottom, width, height]
                % Stretch horizontally (full width), fixed vertical height
                fixedHeight = 500;   % you can adjust this
                bottomPos   = 100;   % vertical position from bottom
                if lr_combine==1
                    figurePos = [0, bottomPos, screenSize(3), fixedHeight];
                else
                    figurePos = [0, 0, screenSize(3), screenSize(4)];
                end
                set(figureHandle, 'Position', figurePos);

                % Force update before saving
                drawnow;

                % Save figure as PNG
                figdir = fullfile(project_dir,'bids/derivatives/figures');
                filename = fullfile(figdir,sprintf('fig_%s.png',output_prefix));
                exportgraphics(figureHandle, filename, 'Resolution', 300);  % High-res save
            end
        end
    end
end