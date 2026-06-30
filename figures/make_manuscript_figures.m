% make_manuscript_figures.m
% ----------------------------------------------------------------------
% Companion code for:
%   Sung D, et al. Microstructural and exchange imaging with oscillating
%   gradient spin-echo (OGSE) diffusion MRI. Magn Reson Med, 2025.
%   DOI: 10.1002/mrm.70300
%
% Description:
%   Generates the main manuscript and supplementary figures: gradient
%   waveforms and power spectra (Fig 1A/2), Kvar/h[q]/MK frequency
%   dependence in fixed-T and fixed-N regimes (Fig 1B/6, Fig S5),
%   q(w) power spectra and frequency-matching checks, the mixing-time
%   optimization curves (Fig S7), and the group-level MD/MK summary
%   across subjects (Fig S2).
%
% Author: Dongsuk Sung (dsung2@mgh.harvard.edu)
% Athinoula A. Martinos Center for Biomedical Imaging, MGH/Harvard Medical School
%
% Dependencies (not included in this repo):
%   - gacelle toolbox   https://github.com/kschan0214/gacelle
%   - util_ogse         custom OGSE waveform utility class (see README.md)
%   - MATLAB Image Processing Toolbox (niftiread/niftiinfo)
%
% Note: this script is organized in MATLAB "%%" cell blocks, each
% producing one figure (or figure panel). Run cells individually rather
% than the whole file end-to-end, as several cells expect variables
% defined by earlier cells and some define overlapping variable names
% (e.g. Tm_list, project_dir) intentionally reused across figures.
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

%% Initial set up
figdir = fullfile(project_dir,'figures/Manuscript');
ut = util_ogse;

%% Figures 1A, Figure 2 - Gradient Waveforms and Power Spectrum
close all;
% protocol_list = {'acqs-all','acqs-six','acqs-six2','acqs-ten','acqs-four'};
protocol = 'fixed-Tm';
% acq={'30Hz_N1','40Hz_N1','50Hz_N1','60Hz_N1',''65Hz_N2','90Hz_N3'};
wave_idx = 1;
do_overlap = 0;
do_figsave = 0;
with_legend = 0;
plot_waveform(protocol, wave_idx, figdir, do_overlap, do_figsave, with_legend)
disp('Done plotting Grad Wavefrom and Power Spectrum')

%% Figure 1B, Figure 6 - Kvar, hq, MK change with frequency in fixed-T and fixe-N regimes
% figure(1): Fiugre 6A
% figure(2): Figure 6B
% figure(3): Figure 1B - waveform_type = 'cosine'
% figure(4): Figure S5
close all
waveform_type = 'costrap'; % cosine vs costrap
do_figsave = true;
ut = util_ogse;
regime_list = {'Fixed T','Fixed N'};
screenSize = get(0, 'ScreenSize');  % [left bottom width height]
fixedHeight = 500;   % you can adjust this
fixedWidth   = 1000;   % vertical position from bottom
figurePos = [0, 0, fixedWidth,fixedHeight];
cmap = colormap('lines');

% Initialize model parameters in WM
tex_WM = 140.84; 
f_WM = 0.62;
aD_WM = 2.98;
aC_WM = 6.68;
Kinf_WM = 0.02;

% Initialize model parameters in GM
tex_GM = 13.76;
f_GM = 0.43; 
aD_GM = 1.70;
aC_GM = 10.22;    
Kinf_GM = 0.43;

% Initialize other parameters
bval = 2;
tr = 1.59;
Nt = 1e4;

f_select = {[23.21, 62.53, 88.27],[23.21, 34.03, 42.15, 51.61]};
DT_eff = cell(2,1);
for ri = 1:2
    regime = regime_list{ri};  
    switch regime
        case 'Fixed T'
            N_list = (1:4)';
            freq_list = 30*N_list;
        case 'Fixed N'
            Ni = 1;
            freq_list = (30:120)';
            N_list = Ni*ones(numel(freq_list),1);
    end
    
    p2_list   = 1000./(2*freq_list)-2*tr;
    p1_list   = (p2_list+tr)/2-tr;
    G_list = sqrt(bval)./sqrt(4*(1/40*0.0107).^2.*(( (1/3)*p1_list.^3 + (3/2)*(p1_list.^2.*tr) + (23/12)*(p1_list*tr^2) + (23/30)*(tr^3)) + ...
                                  (2*N_list-1) .*( (1/3)*p1_list.^3 + (3/2)*tr*p1_list.^2 + (23/12)*tr^2*p1_list + (91/120)*tr^3)));
    %%%%%%%%%%%%%%%%%%%%%% Find optimal Tm %%%%%%%%%%%%%%%%%%%%%%%%%%%
    t1 = zeros(length(N_list),Nt);
    gt1 = zeros(length(N_list),Nt);
    w1 = zeros(length(N_list),1);
    k = 1:20;
    Tm_Opt = zeros(numel(N_list),numel(k));
    f1 = zeros(numel(N_list),1);
    for i = 1:length(N_list)
        Nh = N_list(i);
        freq = freq_list(i);
        Gmax = G_list(i);
        % Create cosine trapazoidal waveform
        [t1(i,:),gt1(i,:)] = ut.costrap_waveform_half(Nh,freq,Gmax,tr,Nt);
        [w_tmp,qw2_tmp] = ut.convQw(t1(i,:)', gt1(i,:)');
        [~,loc] = findpeaks(qw2_tmp,w_tmp,'SortStr','descend');
        f1(i) = 1000*abs(loc(1))/(2*pi);
        Tm_Opt(i,:) = 1000*k/(f1(i)) - 1000*Nh/freq - tr - 4/Nh; % 4/Nh is an adjustment term to match with manual findings
    end
    
    Tm_Opt(Tm_Opt<7.4) = 10000;
    Tm_Opt = min(Tm_Opt,[],2);
    Tm_Opt(Tm_Opt>40) = 7.4;
    Tm_list = Tm_Opt;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Tm_list = 1*ones(numel(N_list),1);
    T_list = 2*(1000*N_list./freq_list + tr)+Tm_list;
    
    t = zeros(length(N_list),Nt);
    gt = zeros(length(N_list),Nt);
    qw2 = zeros(length(N_list),1);
    w = zeros(length(N_list),1);
    for i = 1:length(N_list)
        Nh = N_list(i);
        freq = freq_list(i);
        Tm = Tm_list(i);
        Gmax = G_list(i);
        % Create cosine trapazoidal waveform
        [t(i,:),gt(i,:)] = ut.costrap_waveform(Nh,freq,Tm,Gmax,tr,Nt);
        [w_tmp,qw2_tmp] = ut.convQw(t(i,:)', gt(i,:)');
        qw2(i) = 1./(2*pi*bval).*trapz(w_tmp, sqrt(abs(w_tmp)).*qw2_tmp');
        [~,loc] = findpeaks(qw2_tmp,w_tmp,'SortStr','descend');
        w(i) = abs(loc(1));
    end
    
    w1_rng = find((1000*w/(2*pi)>=23) & (1000*w/(2*pi)<52));
    w2_rng = find((1000*w/(2*pi)>=52) & (1000*w/(2*pi)<89));
    w3_rng = find((1000*w/(2*pi)>=23) & (1000*w/(2*pi)<89));


    % Figure 6B. h[q]
    figureHandle1 = figure(1);
    subplot(1,2,ri);
    hq = zeros(size(t,1),2);
    hq_new = zeros(freq_list(end)-freq_list(1)+1,2);
    for ti = 1:2
        if ti==1; tex = tex_WM; else; tex = tex_GM; end
        switch waveform_type
            case 'costrap'
                hq(:,ti) = ut.KMhq(t, gt, 1/tex); 
            case 'cosine'
                hq(:,ti) = ut.KMhqOG(T_list, N_list, 1/tex);
        end
        if ri==1
            w_new = linspace(w(1),w(end),size(hq_new,1));
            w1_rng = find((1000*w_new/(2*pi)>=23) & (1000*w_new/(2*pi)<52));
            w2_rng = find((1000*w_new/(2*pi)>=52) & (1000*w_new/(2*pi)<89));
            w3_rng = find((1000*w_new/(2*pi)>=23) & (1000*w_new/(2*pi)<89));
            hq_new(:,ti) = interp1(w,hq(:,ti),w_new,'spline');
            
            f_list = 1000*w_new/(2*pi);
            closest_points = zeros(size(f_select{ri}));
            for ii = 1:length(f_select{ri})
                [~,idx] = min(abs(f_list-f_select{ri}(ii)));
                closest_points(ii) = idx;
            end
            plot(1000*w_new(w3_rng)/(2*pi),hq_new(w3_rng,ti),'-','LineWidth',3, 'Color', cmap(ti,:));
            hold on, scatter(1000*w_new(closest_points)/(2*pi),...
                hq_new(closest_points,ti),100,'o','MarkerEdgeColor',cmap(ti,:),'LineWidth',2);
            legend({'WM','','Cortex'},'Location','northeast')

        else
            plot(1000*w(w1_rng)/(2*pi),hq(w1_rng,ti),'-','LineWidth',3,'Color',cmap(ti,:)); hold on;
            plot(1000*w(w2_rng)/(2*pi),hq(w2_rng,ti),'--','LineWidth',3,'Color',cmap(ti,:));
            f_list = 1000*w/(2*pi);
            closest_points = zeros(size(f_select{ri}));
            for ii = 1:length(f_select{ri})
                [~,idx] = min(abs(f_list-f_select{ri}(ii)));
                closest_points(ii) = idx;
            end
            hold on, scatter(1000*w(closest_points)/(2*pi),hq(closest_points,ti),...
                100,'o','MarkerEdgeColor',cmap(ti,:),'LineWidth',2);

            legend({'WM','','','Cortex'},'Location','northeast')
        end
        hold on;
    end
    xlabel('Frequency [Hz]')
    ylabel('h[q]')
    xlim([10,100]);
    ylim([0,1.2]);
    xticks([30,60,90])
    set(gca,'FontSize',20)
    title(regime);
    set(figureHandle1, 'Position', figurePos);
    drawnow;


    % Figure 6A. Kvar
    switch waveform_type
        case 'costrap'
            Kvar_WM = ut.KvarCosTrap(qw2, f_WM, aD_WM, aC_WM);   
            Kvar_GM = ut.KvarCosTrap(qw2, f_GM, aD_GM, aC_GM);
        case 'cosine'
            Kvar_WM = ut.KvarCosine(N_list, T_list, f_WM, aD_WM, aC_WM);
            Kvar_GM = ut.KvarCosine(N_list, T_list, f_GM, aD_GM, aC_GM);
    end

    figureHandle2 = figure(2);
    subplot(1,2,ri)
    if ri==1
        Kvar_WM_new = interp1(w,Kvar_WM,w_new,'spline');
        plot(1000*w_new(w3_rng)/(2*pi),Kvar_WM_new(w3_rng),'-','LineWidth',3,'Color',cmap(1,:));
        f_list = 1000*w_new/(2*pi);
        closest_points = zeros(size(f_select{ri}));
        for ii = 1:length(f_select{ri})
            [~,idx] = min(abs(f_list-f_select{ri}(ii)));
            closest_points(ii) = idx;
        end
        hold on, scatter(1000*w_new(closest_points)/(2*pi),...
                Kvar_WM_new(closest_points),100,'o','MarkerEdgeColor',cmap(1,:),'LineWidth',2);
    else
        plot(1000*w(w1_rng)/(2*pi),Kvar_WM(w1_rng),'-','LineWidth',3,'Color',cmap(1,:)); hold on;
        plot(1000*w(w2_rng)/(2*pi),Kvar_WM(w2_rng),'--','LineWidth',3,'Color',cmap(1,:));
        f_list = 1000*w/(2*pi);
        closest_points = zeros(size(f_select{ri}));
        for ii = 1:length(f_select{ri})
            [~,idx] = min(abs(f_list-f_select{ri}(ii)));
            closest_points(ii) = idx;
        end
        hold on, scatter(1000*w(closest_points)/(2*pi),...
                Kvar_WM(closest_points),100,'o','MarkerEdgeColor',cmap(1,:),'LineWidth',2);
    end
    hold on;
    if ri==1
        Kvar_GM_new = interp1(w,Kvar_GM,w_new,'spline');
        plot(1000*w_new(w3_rng)/(2*pi),Kvar_GM_new(w3_rng),'-','LineWidth',3,'Color',cmap(2,:));
        hold on, scatter(1000*w_new(closest_points)/(2*pi),...
                Kvar_GM_new(closest_points),100,'o','MarkerEdgeColor',cmap(2,:),'LineWidth',2);
        legend({'WM','','Cortex'},'Location','northeast')
    else
        plot(1000*w(w1_rng)/(2*pi),Kvar_GM(w1_rng),'-','LineWidth',3,'Color',cmap(2,:)); hold on;
        plot(1000*w(w2_rng)/(2*pi),Kvar_GM(w2_rng),'--','LineWidth',3,'Color',cmap(2,:));
        hold on, scatter(1000*w(closest_points)/(2*pi),...
                Kvar_GM(closest_points),100,'o','MarkerEdgeColor',cmap(2,:),'LineWidth',2);
        legend({'WM','','','Cortex'},'Location','northeast')
    end
    xlabel('Frequency [Hz]')
    ylabel('K_{var}')
    xlim([10,100]);
    ylim([0,1.2]);
    set(gca,'FontSize',20)
    xticks([30,60,90])
    title(regime);
    set(figureHandle2, 'Position', figurePos);
    drawnow;

    % Figure 1B
    if ri==2; w_new = w; end
    
    f_list = 1000*w_new/(2*pi);
    f_select_all = [34.03, 42.15, 51.61, 67.24, 76.55, 88.27];
    closest_points = zeros(length(f_select_all),1);
    for ii = 1:length(f_select_all)
        [~,idx] = min(abs(f_list-f_select_all(ii)));
        closest_points(ii) = idx;
    end
    w3_rng = find((1000*w_new/(2*pi)>=33) & (1000*w_new/(2*pi)<89));
    figureHandle3 = figure(3);
    subplot(1,2,ri)
    if ri==1
        plot(1000*w_new(w3_rng)/(2*pi),Kvar_WM_new(w3_rng)'.*hq_new(w3_rng,1)+Kinf_WM,'-','LineWidth',3,'markersize', 18); hold on;
        plot(1000*w_new(w3_rng)/(2*pi),Kvar_GM_new(w3_rng)'.*hq_new(w3_rng,2)+Kinf_GM,'-','LineWidth',3,'markersize', 18);
        hold on, scatter(1000*w_new(closest_points)/(2*pi),...
            Kvar_WM_new(closest_points)'.*hq_new(closest_points,1)+Kinf_WM,100,'o','MarkerEdgeColor',cmap(1,:),'LineWidth',2);
        hold on, scatter(1000*w_new(closest_points)/(2*pi),...
            Kvar_GM_new(closest_points)'.*hq_new(closest_points,2)+Kinf_GM,100,'o','MarkerEdgeColor',cmap(2,:),'LineWidth',2);

    else
        plot(1000*w_new(w3_rng)/(2*pi),Kvar_WM(w3_rng).*hq(w3_rng,1)+Kinf_WM,'-','LineWidth',3,'markersize', 18); hold on;
        plot(1000*w_new(w3_rng)/(2*pi),Kvar_GM(w3_rng).*hq(w3_rng,2)+Kinf_GM,'-','LineWidth',3,'markersize', 18);
        hold on, scatter(1000*w_new(closest_points)/(2*pi),...
            Kvar_WM(closest_points).*hq(closest_points,1)+Kinf_WM,100,'o','MarkerEdgeColor',cmap(1,:),'LineWidth',2);
        hold on, scatter(1000*w_new(closest_points)/(2*pi),...
            Kvar_GM(closest_points).*hq(closest_points,2)+Kinf_GM,100,'o','MarkerEdgeColor',cmap(2,:),'LineWidth',2);

    end
    xlabel('Frequency [Hz]')
    ylabel('Mean Kurtosis')
    xticks([30,60,90])
    ylim([0.5,0.9]);
    set(gca,'FontSize',20)
    xlim([20,100]);
    % if ri==1
        legend({'Long t_{ex}','Short t_{ex}'},'Location','northeast')
    % else
    %     legend({'Long t_{ex}','Short t_{ex}'},'Location','southeast')
    % end
    title(regime);
    set(figureHandle3, 'Position', figurePos);
    drawnow;

    
    
    % Figure S5. Effective diffusion time
    tau = 1000*N_list./freq_list + tr;
    if ri==1
        DT_eff{ri} = tau.*(1-15./(4*w.^2.*tau.^2));
        DT_eff_new = interp1(w,DT_eff{ri},w_new,'spline');
    else
        DT_eff{ri} = (2*pi*N_list./w).*(1-15./(6*pi^2.*N_list.^2));
    end
    figureHandle4 = figure(4);
    subplot(1,2,ri)
    if ri==1
        plot(1000*w_new(w3_rng)/(2*pi),DT_eff_new(w3_rng),'-','LineWidth',3,'markersize', 18); hold on;
        title('Fixed-T')
    else
        plot(1000*w(w3_rng)/(2*pi),DT_eff{ri}(w3_rng),'-','LineWidth',3,'markersize', 18); hold on;
        title('Fixed-N')
    end
    xlabel('Frequency [Hz]');
    ylabel('Effective DT [ms]');
    xticks([30,60,90]);
    ylim([5,40]);
    set(gca,'FontSize',20)
    xlim([20,100]);
    set(figureHandle4, 'Position', figurePos);
    drawnow;
end

if do_figsave
    filename1 = fullfile(figdir,sprintf('hq_fixedT_fixedN_GM_WM_30-90Hz_fix_%s.png',waveform_type));
    filename2 = fullfile(figdir,sprintf('Kvar_fixedT_fixedN_GM_WM_30-90Hz_fix_%s.png',waveform_type));
    filename3 = fullfile(figdir,sprintf('MK_fixedT_fixedN_short_long_tex_30-90Hz_fix_%s.png',waveform_type));
    filename4 = fullfile(figdir,sprintf('DTeff_fixedT_fixedN_30-90Hz_fix_%s.png',waveform_type));
    
    exportgraphics(figureHandle1, filename1, 'Resolution', 300);  % High-res save
    exportgraphics(figureHandle2, filename2, 'Resolution', 300);  % High-res save
    exportgraphics(figureHandle3, filename3, 'Resolution', 300);  % High-res save
    exportgraphics(figureHandle4, filename4, 'Resolution', 300);  % High-res save
end

%% Figure S6. Figures for Response about fixed-Tm vs optimized-Tm
close all
% Initial setup
is_fixTm = true;
feature={'30Hz-N1','40Hz-N1','50Hz-N1',...
            '60Hz-N1','60Hz-N2','90Hz-N3'};
tr=1.59;  % Ramp Time [ms]
N_list    = [1,1,1,1,2,3];
% Parameters
if ~is_fixTm
    % freq_list = [30,40,50,60,65,90];
    % Gmax_list = [136.85, 213.20, 302.70, 405.06, 326.54, 464.59];
    % Tm_list   = [7.4, 30.9, 24.7, 20.4, 15.2, 10.1]; % Mixing Time [ms]  
    req_list = [30,40,50,60,60,90];
    Gmax_list = [136.85, 213.20, 302.70, 405.06, 286.48, 464.59];
    Tm_list   = [7.4, 30.9, 24.7, 20.4, 16.5, 10.1]; % Mixing Time [ms
else
    freq_list = [30,40,50,60,60,90];
    Gmax_list = [136.85, 213.20, 302.70, 405.06, 286.48, 464.59];
    Tm_list   = 10.1*ones(size(N_list));
end

T_tot = 2*(1000*N_list./freq_list + tr)+Tm_list;
fprintf('Total time = %.2f ms\n',T_tot)
disp('%%%%%%%%%%%%%')

% Cosine Trapazoidal Waveform with Tm
ut = util_ogse;
Nt = 1e4;           % number of time steps
t_all = zeros(Nt,length(N_list));
gt_all = zeros(Nt,length(N_list));

figure();
for i = 1:length(N_list)
    N = N_list(i);
    Tm = Tm_list(i);
    freq = freq_list(i);
    Gmax = Gmax_list(i);

    % Create cosine trapazoidal waveform
    [t,gt] = ut.costrap_waveform(N,freq,Tm,Gmax,tr,Nt);
    t_all(:,i) = t;
    gt_all(:,i) = gt;

end

% Transform the waveform to frequency domain
w_list = zeros(size(N_list));
b_list = zeros(size(N_list));
figure;

cnt = 1;
selection = [1,5,6,2,3,4];
for i =  1:numel(N_list)
    si = selection(i);
    [w,qw2,b_list(si)] = ut.convQw(1e-3*t_all(:,si),gt_all(:,si));
    % [w,qw2,b_list(si)] = ut.convQw(t_all(:,si),gt_all(:,si));
    [pks,loc] = findpeaks(qw2,w,'SortStr','descend');
    fprintf('%s actual frequency = %.2f\n',feature{si},abs(loc(1))/(2*pi));
    w_list(si) = abs(loc(1));
    subplot(2,ceil(length(selection)/2),cnt),plot(w/(2*pi),qw2','LineWidth',5,'Color','k');
    title(feature{si})
    xlabel('f [Hz]')
    ylabel('|Q(f)|^2')
    xlim([-150,150])
    set(gca,'FontSize',20)
    % xlim([-abs(loc(1))*1.5/pi,abs(loc(1))*1.5/pi]);
    cnt = cnt + 1;
end

disp('%%%%%%%%%%%%%%%%%%%%%%%%%')


% Compare angular frequency (w) acquired in two regimes
figure,plot(w_list/(2*pi),'o-')
hold on, plot((2*N_list)./(T_tot/1000),'*-','LineWidth',2);
legend('w_p/2\pi','2\pi*(2*N)/T', 'Location','northwest')
ylabel('Frequency [Hz]');
xticks([1:6]);
xticklabels({'30Hz-N1','40Hz-N1','50Hz-N1','60Hz-N1','60Hz-N2','90Hz-N3'});
set(gca,'FontSize',20)
xlim([0,8]);

figure,plot((2*N_list)./(T_tot/1000),w_list/(2*pi),'o-','LineWidth',2)
x = linspace(0,100,100);
y = x;
hold on, plot(x,y,'LineWidth',2);
ylabel('w_p/2\pi');
xlabel('2\pi*(2*N)/T')
set(gca,'FontSize',20)
legend('data points','y=x line', 'Location','northwest')

%% Figure S7
close all
waveform_type = 'costrap'; % cosine vs costrap
do_figsave = true;
ut = util_ogse;

% Initialize model parameters in WM
tex_WM = 140.84; 
f_WM = 0.62;
aD_WM = 2.98;
aC_WM = 6.68;
Kinf_WM = 0.02;

% Initialize model parameters in GM
tex_GM = 13.76;
f_GM = 0.43; 
aD_GM = 1.70;
aC_GM = 10.22;    
Kinf_GM = 0.43;

% Initialize other parameters
bval = 2;
tr = 1.59;
Nt = 1e4;


qw2 = zeros(length(N_list),1);
w = zeros(length(N_list),1);

Nh = 1;
freq = 30;
Gmax = 136.85;
Tm_list = 7:0.5:42;
hq_GM = zeros(numel(Tm_list),1);
hq_WM = zeros(numel(Tm_list),1);
Kvar_GM = zeros(numel(Tm_list),1);
Kvar_WM = zeros(numel(Tm_list),1);

for i = 1:numel(Tm_list)

    Tm = Tm_list(i);
    fprintf('Tm = %.1f\n',Tm);
    % Create cosine trapazoidal waveform
    [t,gt] = ut.costrap_waveform(Nh,freq,Tm,Gmax,tr,Nt);
    [w_tmp,qw2_tmp] = ut.convQw(t, gt);
    qw2(i) = 1./(2*pi*bval).*trapz(w_tmp, sqrt(abs(w_tmp)).*qw2_tmp');
    [~,loc] = findpeaks(qw2_tmp,w_tmp,'SortStr','descend');


    switch waveform_type
        case 'costrap'
            hq_GM(i) = ut.KMhq(t', gt', 1/tex_GM); 
            hq_WM(i) = ut.KMhq(t', gt', 1/tex_WM); 
            Kvar_WM(i) = ut.KvarCosTrap(qw2(i), f_WM, aD_WM, aC_WM);   
            Kvar_GM(i) = ut.KvarCosTrap(qw2(i), f_GM, aD_GM, aC_GM);
        case 'cosine'
            hq_GM(i) = ut.KMhqOG(T_list, N_list, 1/tex_GM);
            hq_WM(i) = ut.KMhqOG(T_list, N_list, 1/tex_WM);
            Kvar_WM(i) = ut.KvarCosine(N_list, T_list, f_WM, aD_WM, aC_WM);
            Kvar_GM(i) = ut.KvarCosine(N_list, T_list, f_GM, aD_GM, aC_GM);
    end
end

figure,plot(Tm_list,hq_GM.*Kvar_GM,'o-','LineWidth',2);
xlabel('Mixing Time [ms]')
ylabel('Kvar*hq')
set(gca,'FontSize',20)
xlim([5,45])
figure,plot(Tm_list,hq_WM.*Kvar_WM,'o-','LineWidth',2);
xlabel('Mixing Time [ms]')
ylabel('Kvar*hq')
set(gca,'FontSize',20)
xlim([5,45])

%% Figure S2. MK and MD of all subjects
% close all
subject_list = {'sub-ogse005-02','sub-ogse006','sub-ogse007','sub-ogse008','sub-ogse009','sub-ogse010',...
    'sub-ogse011','sub-ogse012','sub-ogse013','sub-ogse014','sub-ogse015'};
Nsub = length(subject_list);
Nfreq = length(N_list);
median_md = zeros(Nfreq,2,Nsub);
median_mk = zeros(Nfreq,2,Nsub);

list_fixN1 = [1,2,3,4];
list_fixT = [1,5,6];


for sub_i = 1:length(subject_list)
    disp(subject_list{sub_i})
    subj_label = subject_list{sub_i};
    dwi_proc_dirs = fullfile(project_dir,'bids/derivatives/mrtrix_tensor',subj_label);
    synthseg_dir = fullfile(project_dir,'bids/derivatives/synthseg',subj_label);

    mask_Cortex = niftiread(fullfile(synthseg_dir,sprintf('%s_cortex_mask_res.nii.gz',subj_label)));
    mask_Cortex = mask_Cortex>0;
    mask_WM = niftiread(fullfile(synthseg_dir,sprintf('%s_wm_mask_res.nii.gz',subj_label)));
    mask_WM = mask_WM>0;

    for i = 1:length(acq)
        %%% Load MD & MK %%%
        md = 1000*niftiread(fullfile(dwi_proc_dirs,sprintf('%s_acq-%s_dwi_adc.nii.gz',subj_label,acq{i})));
        mk = niftiread(fullfile(dwi_proc_dirs,sprintf('%s_acq-%s_dwi_mk.nii.gz',subj_label,acq{i})));
    
        %%% Remove outliers %%%
        lower_lim_mk = mk<0;
        upper_lim_mk = mk>3;    
        mk(lower_lim_mk) = NaN; %0;
        mk(upper_lim_mk) = NaN; %3   

        % lower_lim_md = md<0;
        % upper_lim_md = md>2;    
        md(lower_lim_mk) = NaN; %0;
        md(upper_lim_mk) = NaN; %3   
    
        % WM
        md_roi = md.*mask_WM;
        median_md(i,1,sub_i) = median(nonzeros(md_roi(:)),'omitnan');
        mk_roi = mk.*mask_WM;
        median_mk(i,1,sub_i) = median(nonzeros(mk_roi(:)),'omitnan');
        
        % Cortex
        md_roi = md.*mask_Cortex;
        median_md(i,2,sub_i) = median(nonzeros(md_roi(:)),'omitnan');
        mk_roi = mk.*mask_Cortex;
        median_mk(i,2,sub_i) = median(nonzeros(mk_roi(:)),'omitnan');
    end
end

% draw Plot
RegionLabel = {'WM','Cortex'};
for ri = 1:2
    if ri==1
        MD_Reg = squeeze(median_md(:,1,:));
        MK_Reg = squeeze(median_mk(:,1,:));
    else
        MD_Reg = squeeze(median_md(:,2,:));
        MK_Reg = squeeze(median_mk(:,2,:));
    end

    figure(1)
    subplot(1,2,ri)
    cmap = colormap('lines');
    hold on;
    h_fixT_y = plot(true_freq_list(list_fixT), MD_Reg(list_fixT,:), 'v-', 'Color', cmap(1,:), 'markersize', 18, 'linewidth', 2);
    h_fixN_y1 = plot(true_freq_list(list_fixN1), MD_Reg(list_fixN1,:), '.-', 'Color', cmap(2,:), 'markersize', 40,'linewidth',2);

    xlim([0,100])
    ylim([0.9, 1.3]);

    xlabel('Frequency [Hz]')
    ylabel('Mean Diffusivity')
    title(sprintf('%s',RegionLabel{ri}));
    set(gca,'FontSize',20)
    % legend({'Fixed T(=80ms), data','Fixed N(=1), data'},'Location','southeast')

    figure(2)
    subplot(1,2,ri)
    cmap = colormap('lines');
    hold on;
    h_fixT_y = plot(true_freq_list(list_fixT), MK_Reg(list_fixT,:), 'v-', 'Color', cmap(1,:), 'markersize', 18, 'linewidth', 2);
    h_fixN_y1 = plot(true_freq_list(list_fixN1), MK_Reg(list_fixN1,:), '.-', 'Color', cmap(2,:), 'markersize', 40,'linewidth',2);

    xlim([0,100])
    ylim([0.45, 0.85]);

    xlabel('Frequency [Hz]')
    ylabel('Mean Kurtosis')
    title(sprintf('%s',RegionLabel{ri}));
    set(gca,'FontSize',20)
end