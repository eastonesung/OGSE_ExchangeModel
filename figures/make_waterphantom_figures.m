% make_waterphantom_figures.m
% ----------------------------------------------------------------------
% Companion code for:
%   Sung D, et al. Microstructural and exchange imaging with oscillating
%   gradient spin-echo (OGSE) diffusion MRI. Magn Reson Med, 2025.
%   DOI: 10.1002/mrm.70300
%
% Description:
%   Generates the water-phantom validation figures (mean diffusivity and
%   mean kurtosis vs. effective frequency) used to confirm that the
%   fixed-Tm / fixed-N acquisition protocol behaves as expected in an
%   isotropic medium with no exchange.
%
% Author: Dongsuk Sung (dsung2@mgh.harvard.edu)
% Athinoula A. Martinos Center for Biomedical Imaging, MGH/Harvard Medical School
%
% Dependencies (not included in this repo):
%   - gacelle toolbox   https://github.com/kschan0214/gacelle
%   - util_ogse         custom OGSE waveform utility class (see README.md)
%   - MATLAB Image Processing Toolbox (niftiread/niftiinfo)
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

%% Initial setting
subj_label = 'sub-ogse0pw';
phantom_type = 'water';
acq={'30Hz_N1','40Hz_N1','50Hz_N1','60Hz_N1','65Hz_N2','90Hz_N3'};
N_list    = [1,1,1,1,2,3];
freq_list = [30,40,50,60,65,90];                 % [Hz]
Tm_list   = [7.4, 30.9, 24.7, 20.4, 15.2, 10.1]; % Mixing Time [ms]
list_fixN = [1,2,3,4];
list_fixT = [1,5,6];
dwi_proc_dirs = fullfile(project_dir,'bids/derivatives/mrtrix_tensor',subj_label);

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
model_obj = gpu_KurtModel(t_vec, gt_vec, qw2_int);
%% LOAD Mask and b=0
preproc_dir = fullfile(project_dir,'bids/derivatives/preprocessed',subj_label);
synthseg_dir = fullfile(project_dir,'bids/derivatives/synthseg',subj_label);
nii_info = niftiinfo(fullfile(preproc_dir,sprintf('%s_dwi_brain_mask_gncorr_C2.nii.gz',subj_label)));
mask     = niftiread(nii_info)>0;
b0 = niftiread(fullfile(synthseg_dir,sprintf('%s_acq-30Hz_N1_dwi_denoise_degibbs_eddy_gncorr_b0_Tmean.nii.gz',subj_label)));
%% Create mask
switch phantom_type
    case 'aniso'
        x_range = 39:41;
        y_range = 38:40;
        z_range = 2:12;
        mask_r = zeros(size(mask));
        mask_r(x_range,y_range,z_range) = 1;
        
        x_range = 32:47;
        y_range = 45:60;
        z_range = 2:12;
        mask_w = zeros(size(mask));
        mask_w(x_range,y_range,z_range) = 1;
    case 'water'
        x_range = 33:47;
        y_range = 36:50;%31:45;
        z_range = 2:16;
        mask_r = zeros(size(mask));
        mask_w = zeros(size(mask));
        mask_w(x_range,y_range,z_range) = 1;
end
%% Load measured MD & MK
Nfreq = length(N_list);
median_md = zeros(Nfreq,3);
median_mk = zeros(Nfreq,3);
for i = 1:length(acq)
    %%% Load MD & MK %%%
    md = 1000*niftiread(fullfile(dwi_proc_dirs,sprintf('%s_acq-%s_dwi_adc.nii.gz',subj_label,acq{i})));
    mk = niftiread(fullfile(dwi_proc_dirs,sprintf('%s_acq-%s_dwi_mk.nii.gz',subj_label,acq{i})));

    %%% Remove outliers %%%
    lower_lim = mk<0;
    upper_lim = mk>3;    
    mk(lower_lim) = NaN; %0;
    mk(upper_lim) = NaN; %3   

    % Whole brain mask
    md_roi = md.*mask;
    median_md(i,1) = median(nonzeros(md_roi(:)),'omitnan');
    mk_roi = mk.*mask;
    median_mk(i,1) = median(nonzeros(mk_roi(:)),'omitnan');

    % Restricted mask
    md_roi = md.*mask_r;
    median_md(i,2) = median(nonzeros(md_roi(:)),'omitnan');
    mk_roi = mk.*mask_r;
    median_mk(i,2) = median(nonzeros(mk_roi(:)),'omitnan');

    % Water mask
    md_roi = md.*mask_w;
    median_md(i,3) = median(nonzeros(md_roi(:)),'omitnan');
    mk_roi = mk.*mask_w;
    median_mk(i,3) = median(nonzeros(mk_roi(:)),'omitnan');
    
end
%% Draw plots of MD and MK
% MD
figure;
cmap = colormap('lines');
hold on;
h_fixT_y = plot(true_freq_list(list_fixT), median_md(list_fixT,3), 'v-', 'Color', cmap(1,:), 'markersize', 16, 'linewidth', 2);
h_fixN_y1 = plot(true_freq_list(list_fixN), median_md(list_fixN,3), 'o-', 'Color', cmap(2,:), 'markersize', 16,'linewidth', 2);
ylim([1.1,1.5]);

ylabel('Mean Diffusivity');
xlabel('Frequency [Hz]')
set(gca,'FontSize',16)
legend({'Fixed-T','Fixed-N'})
% MK
figure;
cmap = colormap('lines');
hold on;
h_fixT_y = plot(true_freq_list(list_fixT), median_mk(list_fixT,3), 'v-', 'Color', cmap(1,:), 'markersize', 16, 'linewidth', 2);
h_fixN_y1 = plot(true_freq_list(list_fixN), median_mk(list_fixN,3), 'o-', 'Color', cmap(2,:), 'markersize', 16,'linewidth', 2);
ylim([0,0.001]);

ylabel('Mean Kurtosis');
xlabel('Frequency [Hz]')
set(gca,'FontSize',16)
legend({'Fixed-T','Fixed-N'})
%% MK in different ROIs
figure;
cmap = colormap('lines');
list_fixT = [1,5,6];
list_fixN = 1:4;
title_list = {'all','restricted','water-only'};

for ii = 1:3
    subplot(1,3,ii);
    h_fixT_y = plot(true_freq_list(list_fixT), median_mk(list_fixT,ii), 'v-', 'Color', cmap(1,:), 'markersize', 16, 'linewidth', 2);
    hold on;
    h_fixN_y1 = plot(true_freq_list(list_fixN), median_mk(list_fixN,ii), 'o-', 'Color', cmap(2,:), 'markersize', 16,'linewidth', 2);
    ylim([0,0.4]);
    ylabel('Mean Kurtosis');
    xlabel('Frequency [Hz]');
    title(title_list{ii});
    set(gca,'FontSize',16)
    hold on;
end

%% Draw mask
figure
tl = tiledlayout(3,4);
for k = 2:13
    [Br,Lr] = bwboundaries(rot90(mask_r(:,:,k)));
    [Bw,Lw] = bwboundaries(rot90(mask_w(:,:,k)));
    nexttile
    imagesc(rot90(b0(:,:,k))); hold on; colormap gray; clim([0,2000])
    xlim([25,55]); ylim([25,55]);
    if ~isempty(Br)
        boundary1 = Br{1};
        plot(boundary1(:,2), boundary1(:,1), '-','Color',[0.6350 0.0780 0.1840],  'LineWidth', 1)
    end
    if ~isempty(Bw)
        boundary2 = Bw{1};
        plot(boundary2(:,2), boundary2(:,1), '-','Color',[0 0.480 0.6340],  'LineWidth', 1)
    end
    axis('off')
end

%% Draw b=0,1000,2000
figure();
titledlayout(1,3);
imagesc(rot90(b0(:,:,5))); clim([0,1000]);
