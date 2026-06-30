classdef gpu_KurtModel < handle
% gpu_KurtModel  Forward model for OGSE exchange-sensitized kurtosis
%
% ----------------------------------------------------------------------
% Companion code for:
%   Sung D, et al. Microstructural and exchange imaging with oscillating
%   gradient spin-echo (OGSE) diffusion MRI. Magn Reson Med, 2025.
%   DOI: 10.1002/mrm.70300
%
% Description:
%   Implements the Karger-model-derived forward signal for frequency-
%   dependent (OGSE) mean kurtosis as a function of exchange time (tex),
%   neurite/intra-compartment volume fraction (fn), extracellular
%   tortuosity (aD), a frequency-dependence factor (aC), and a partial-
%   volume term (Kinf). Designed to be fit with the askAdam solver
%   (gradient-descent) or MCMC sampler from the gacelle toolbox.
%
% Author: Dongsuk Sung, MGH (dsung2@mgh.harvard.edu)
% Adapted from the gpuNEXI class by Kwok-Shing Chan, MGH (kchan2@mgh.harvard.edu)
% Date created: Jun 12, 2025 (v0.1.0)
%
% Dependency: gacelle toolbox (askadam, utils classes)
%             https://github.com/kschan0214/gacelle
% ----------------------------------------------------------------------

    properties
        % default model parameters and estimation boundary
        % ra        : exchange rate from neurite to extracellular space [1/ms]
        % fa        : Neurite volume fraction
        % Dinf_para : frequency independent parallel diffusivity [ms/us^2]
        % Dinf_ratio : frequency independent perpendicular diffusivity [ms/us^2]       
        % p2        : non-linear neurite dispersion index

        % this is for askAdam
        modelParams     = {'tex';  'fn'; 'aD'; 'aC'; 'Kinf'};
        ub              = [  200;   0.99;    1;    20;    1]
        lb              = [  1;     0.01;   0.1;  0.01; 0.01];
        startPoint      = [  150;    0.7;   0.5;    5;    0.2];
        % modelParams     = {'tex';  'fn'; 'aD'; 'aC'; 'Kinf'};
        % ub              = [  250;   0.99;    1;    20;    1]
        % lb              = [  1;     0.01;   0.1;  0.01; 0.01];
        % startPoint      = [  150;    0.7;   0.5;    5;    0.2];

        % this is for MCMC
        % modelParams     = {'tex';  'fn'; 'aD'; 'aC'; 'Kinf'; 'noise'};
        % ub              = [  200;   0.99;    1;    20;    1;  0.005]
        % lb              = [  1;     0.01;   0.1;  0.01; 0.01; 0.001];
        % startPoint      = [  150;    0.7;   0.5;    5;    0.2; 0.002];
    end

    properties (GetAccess = public, SetAccess = protected)
        t_vec;
        gt_vec;
        qw2_int;
        Nav
    end
    
    methods

        % constructuor
        function this = gpu_KurtModel(t_vec, gt_vec, qw2_int, varargin)
        % NEXI Exchange rate estimation using NEXI model
        % obj = gpuNEXI(b, Delta, Nav)
        %
        % Input
        % ----------
        % b         : b-value [ms/um2]
        % N         : Number of oscillations
        % Nav       : # gradient direction for each b-shell (optional)
        %
        % Output
        % ----------
        % obj       : object of a fitting class
        %
        % Usage
        % ----------
        % obj                   = NEXI(b, Delta, Nav);
        % [out, fa, Da, De, r]  = obj.fit(S, mask, fitting,);
        % Sfit                  = smt.FWD([fa, Da, De, r]);
        % [x_train, S_train]    = obj.traindata(1e4);
        % pars0                 = smt.likelihood(S, x_train, S_train);
        % [out, fa, Da, De, r]  = smt.fit(S, mask, fitting, pars0);
        %
        % Author:
        %  Kwok-Shing Chan (kchan2@mgh.harvard.edu) 
        %  Hong-Hsi Lee (hlee84@mgh.harvard.edu)
        %  Copyright (c) 2023 Massachusetts General Hospital
        %
        %  Adapted from the code of
        %  Dmitry Novikov (dmitry.novikov@nyulangone.org)
        %  Copyright (c) 2023 New York University
            
            this.t_vec  = t_vec;
            this.gt_vec = gt_vec;
            this.qw2_int = qw2_int;

            if nargin > 3
                % NOTE: varargin{8} looks unusual for a single optional Nav
                % argument (nargin>3 implies only varargin{1} exists at minimum).
                % Left unchanged from the original analysis code -- please verify
                % this branch is exercised/correct before relying on it.
                this.Nav = varargin{8};
            else
                this.Nav =  ones(size(t_vec,1),1) ;
            end
            this.Nav = this.Nav(:) ;
            
        end
        
        % update properties according to lmax
        function this = updateProperty(this, fitting)

            % DIMWI
            if fitting.lmax == 0
                idx = find(ismember(this.modelParams,'p2'));
                this.modelParams(idx)     = [];
                this.lb(idx)              = [];
                this.ub(idx)              = [];
                this.startPoint(idx)      = [];
            end

        end

        % display some info about the input data and model parameters
        function display_data_model_info(this)

            disp('========================');
            disp('NEXI with askAdam solver');
            disp('========================');

            disp('----------------')
            disp('Input data: t, gt, qw2_int');
            disp('----------------')
            % fprintf('b-shells (ms/um2)              : [%s] \n',num2str(this.b.',' %.2f'));
            % fprintf('N of oscillations              : [%s] \n',num2str(this.N.',' %d'));
            % fprintf('Set frequencies (Hz)           : [%s] \n',num2str(this.freq.',' %d'));
            % fprintf('Mixing times (ms)              : [%s] \n',num2str(this.Tm.',' %.1f'));
            disp('----------------');
        end

        %% higher-level data fitting functions
        % Wrapper function of fit to handle image data; automatically segment data and fitting in case the data cannot fit in the GPU in one go
        function  [out] = estimate(this, dwi, mask, mask_CSF, extradata, fitting, pars0)
        % Perform NEXI model parameter estimation based on askAdam
        % Input data are expected in multi-dimensional image
        % 
        % Input
        % -----------
        % dwi       : 4D DWI, [x,y,z,dwi]
        % mask      : 3D signal mask, [x,y,z]
        % extradata : Optional additional data
        %   .bval       : 1D bval in ms/um2, [1,dwi]                
        %   .bvec       : 2D b-table, [3,dwi]                       
        %   .N          : 1D # of oscillations, [1,dwi]             
        %   .freq       : 1D input frequencies in Hz, [1,dwi]     
        %   .Tm         : 1D mixing times in ms, [1,dwi] 
        %   .sigma      : 3D noise map, [x,y,z]                     (Optional, only needed for NEXIrice model)
        % fitting   : fitting algorithm parameters (see fit function)
        % pars0     : (Optional) initial starting points for model parameters
        % 
        % Output
        % -----------
        % out       : output structure contains all estimation results
        % tex        : exchange rate from intra- to extra-neurite compartment
        % fn        : Intraneurite volume fraction
        % aD        : extra-cellular tortuosity (Dinf_para/Dinf_perp)
        % aC        : frequency dependence factor (c/Dinf_perp)
        % Kinf      : Partial volume related term
        % 
            
            % display basic info
            this.display_data_model_info;

            % get all fitting algorithm parameters 
            fitting = this.check_set_default(fitting);

            % get matrix size
            dims = size(dwi,1:3);

            %%%%%%%%%%%%%%%% Step 1: Validate all input data %%%%%%%%%%%%%%%%
            % compute rotationally invariant signal if needed
            [dwi,mask] = this.prepare_dwi_data(dwi,mask,extradata,fitting.lmax);

            % mask sure no nan or inf
            [dwi,mask] = utils.remove_img_naninf(dwi,mask);

            % if no pars input at all (not even empty) then use prior
            if nargin < 7; pars0 = []; end

            % convert datatype to single
            dwi     = single(dwi);
            mask    = mask >0;
            if ~isempty(pars0)
                for km = 1:numel(this.modelParams)
                    pars0.(this.modelParams{km}) = single(pars0.(this.modelParams{km})); 
                end
            end

            %%%%%%%%%%%%%%%% End Step 1 %%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%% Step 2: Validate if GPU has enough memory  %%%%%%%%%%%%%%%%
            % determine if we need to divide the data to fit in GPU
            g = gpuDevice(1); reset(g);
            % memoryFixPerVoxel       = 0.2;   % get this number based on mdl fit
            % memoryDynamicPerVoxel   = 2.4;%memoryFixPerVoxel*size(dwi,4);     % get this number based on mdl fit
            % [NSegment,maxSlice]     = utils.find_optimal_divide(mask,memoryFixPerVoxel,memoryDynamicPerVoxel);
            % parameter estimation
            out = [];
            sum_vox = 0;
            for ks = 1:NSegment
                g = gpuDevice(1); reset(g);
                if NSegment ~= 1
                    fprintf('Running #Segment = %d/%d \n',ks,NSegment);
                    disp   ('------------------------')
                end
    
                % determine slice# given a segment
                if ks ~= NSegment
                    slice = 1+(ks-1)*maxSlice : ks*maxSlice;
                else
                    slice = 1+(ks-1)*maxSlice : dims(3);
                end
                
                % divide the data
                dwi_tmp     = dwi(:,:,slice,:);
                mask_tmp    = mask(:,:,slice);
                mask_csf_tmp = mask_CSF(:,:,slice);
                if ~isempty(pars0); for km = 1:numel(this.modelParams); pars0_tmp.(this.modelParams{km}) = pars0.(this.modelParams{km})(:,:,slice); end
                else;               pars0_tmp = [];                 end
                sum_vox = sum_vox + sum(mask_tmp(:));
                fprintf('Running #voxel = %d/%d \n',sum_vox,sum(mask(:)));
                if sum(mask_tmp(:))==0; continue; end
                
                % run fitting
                [out_tmp]    = this.fit(dwi_tmp,mask_tmp,mask_csf_tmp,fitting,pars0_tmp);

                % restore 'out' structure from segment
                out = utils.restore_segment_structure(out,out_tmp,slice,ks);

            end
            out.mask = mask;
            %%%%%%%%%%%%%%%% End Step 2 %%%%%%%%%%%%%%%%

            % save the estimation results if the output filename is provided
            askadam.save_askadam_output(fitting.outputFilename,out)

        end

        % Data fitting function, can be 2D (voxel-based) or 4D (image-based)
        function [out] = fit(this,dwi,mask,mask_CSF,fitting,pars0)
        %
        % Input
        % -----------
        % dwi       : S0 normalised 4D dwi images, [x,y,slice,diffusion], 4th dimension corresponding to [Sl0_b1,Sl0_b2,Sl2_b1,Sl2_b2, etc.]; the order of bval must match the order in the constructor gpuNEXI
        % mask      : 3D signal mask, [x,y,slice]
        % fitting   : fitting algorithm parameters
        %   .Nepoch             : no. of maximum iterations, default = 4000
        %   .initialLearnRate   : initial gradient step size, defaulr = 0.01
        %   .decayRate          : decay rate of gradient step size; learningRate = initialLearnRate / (1+decayRate*epoch), default = 0.0005
        %   .convergenceValue   : convergence tolerance, based on the slope of last 'convergenceWindow' data points on loss, default = 1e-8
        %   .convergenceWindow  : number of data points to check convergence, default = 20
        %   .tol                : stop criteria on metric value, default = 1e-3
        %   .lambda             : regularisation parameter, default = 0 (no regularisation)
        %   .TVmode             : mode for TV regulariation, '2D'|'3D', default = '2D'
        %   .regmap             : parameter map used for regularisation, 'fa'|'ra'|'Da'|'De', default = 'fa'
        %   .lmax               : Order of rotational invariant, 0|2, default = 0
        %   .lossFunction       : loss for data fidelity term, 'L1'|'L2'|'MSE', default = 'L1'
        %   .display            : online display the fitting process on figure, true|false, defualt = false
        % pars0     : 4D parameter starting points of fitting, [x,y,slice,param], 4th dimension corresponding to fitting  parameters with order [fa,Da,De,ra,p2] (optional)
        % 
        % Output
        % -----------
        % out       : output structure
        %   .final      : final results
        %       .tex        : exchange rate from intra- to extra-neurite compartment
        %       .fn        : Intraneurite volume fraction
        %       .Dinf_para : Frequency independent parallel diffusivity (um2/ms)
        %       .c         : constant to calculate intra-neurite parallel diffusivity
        %       .p2        : dispersion index (if fitting.lax=2)
        %       .loss       : final loss metric
        %   .min        : results with the minimum loss metric across all iterations
        %       .tex        : exchange rate from intra- to extra-neurite compartment
        %       .fn        : Intraneurite volume fraction
        %       .aD        : extra-cellular tortuosity (Dinf_para/Dinf_perp)
        %       .aC        : frequency dependence factor (c/Dinf_perp)
        %       .Kinf      : Partial volume related term
        %       .loss       : loss metric      
        % tex        : Final exchange rate from intra- to extra-neurite compartment
        % fn        : Final intraneurite volume fraction
        % Dinf_para : Final frequency independent parallel diffusivity (um2/ms)
        % Dinf_ratio : Final frequency independent perpendicular diffusivity (um2/ms)
        % c         : Final constant to calculate intra-neurite parallel diffusivity
        % p2        : Final dispersion index (if fitting.lax=2)
        %
        % Description: askAdam Image-based NEXI model fitting
        %
        % Kwok-Shing Chan @ MGH
        % kchan2@mgh.harvard.edu
        % Date created: 8 Dec 2023
        % Date modified: 3 April 2024
        %
        %
            
            % check GPU
            gpool = gpuDevice;
            
            % get image size
            dims = size(dwi,1:3);

            %%%%%%%%%%%%%%%%%%%% 1. Validate and parse input %%%%%%%%%%%%%%%%%%%%
            if nargin < 3 || isempty(mask); mask = ones(dims,'logical'); end % if no mask input then fit everthing
            if nargin < 4; fitting = struct(); end
            % set initial tarting points
            if nargin < 5; pars0 = []; % no initial starting points
            else
                if ~isempty(pars0); for km = 1:numel(this.modelParams); pars0.(this.modelParams{km}) = single(pars0.(this.modelParams{km})); end; end
            end

            % get all fitting algorithm parameters 
            fitting                 = this.check_set_default(fitting);
            % determine fitting parameters
            this                    = this.updateProperty(fitting);
            fitting.modelParams     = this.modelParams;
            % set fitting boundary if no input from user
            if isempty( fitting.ub); fitting.ub = this.ub(1:numel(fitting.modelParams)); end
            if isempty( fitting.lb); fitting.lb = this.lb(1:numel(fitting.modelParams)); end
            
            %%%%%%%%%%%%%%%%%%%% End 1 %%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%% 2. Setting up all necessary data, run askadam and get all output %%%%%%%%%%%%%%%%%%%%
            % 2.1 setup fitting weights
            w = this.compute_optimisation_weights(mask,fitting.lossFunction,fitting.lmax); % This is a customised funtion

            % 2.2 estimate prior if needed
            % if isempty(pars0);  pars0 = this.determine_x0(dwi,mask,fitting); end
            if isempty(pars0);  pars0 = this.determine_x0(dwi,mask,mask_CSF,fitting); end

            % You may add more dispay messages here
            disp('---------------------------');
            disp('Additional model parameters');
            disp('---------------------------');
            disp(['KMRK4 rotational invariant model, lmax = ' num2str(fitting.lmax)]);
            disp('---------------------------');

            % 2.3 askAdam optimisation main
            askadamObj = askadam();
            % initiate starting points arrays
            out = askadamObj.optimisation( dwi, mask, w, pars0, fitting, @this.FWD, fitting.lmax);

            %%%%%%%%%%%%%%%%%%%% End 2 %%%%%%%%%%%%%%%%%%%%

            disp('The estimation is completed.');
            
            % clear GPU
            reset(gpool)
            
        end

        %% Data preparation

        % compute weights for optimisation
        function w = compute_optimisation_weights(this,mask,lossFunction,lmax)
        % 
        % Output
        % ------
        % w         : 1D signal masked wegiths
        %

            dims = size(mask,1:3);
            % lmax dependent weights
            l = 0:2:lmax;
            w = zeros([dims size(this.t_vec,1)*numel(l)],'single');
            % w = zeros(dims,'single');
            for kl = 1:(lmax/2+1)
                for kb = 1:size(this.t_vec,1)
                    w(:,:,:,(kl-1)*size(this.t_vec,1)+kb) = this.Nav(kb) / (2*l(kl)+1);
                end
            end
            % if L1 then take square root
            if strcmpi(lossFunction,'l1')
                w = sqrt(w);
            end
            w = w ./ max(w(:));
        end

        % compute rotationally invariant DWI signal if necessary
        function [dwi, mask] = prepare_dwi_data(this,dwi,mask,extradata,lmax)
            % full DWI data then compute rotaionally invariant signal
            if size(dwi,4)/(lmax/2+1) > size(this.t_vec,1) 
                % compute spherical mean signal
                fprintf('Computing rotationally invariant signal...')

                % if the inout N is one value then create a vector
                if isscalar(extradata.N) % previously if numel(extradata.N) == 1
                    extradata.N = ones(size(extradata.bval)) * extradata.N;
                end
                obj = util_ogse;
                [dwi] = obj.get_Sl_all(dwi,extradata.bval,extradata.bvec,extradata.N,extradata.freq,extradata.Tm,lmax);
                fprintf('done.\n');

            elseif size(dwi,4) < size(this.t_vec,1)
                error('There are more b-shells in the class object than available in the input data. Please check your input data.');
            end

            % make sure no NaN/Inf in signal
            [dwi, mask_naninf] = utils.set_nan_inf_zero(dwi); mask_naninf = max(mask_naninf,[],4);
            mask_valid         = and(and(mask, max(dwi,[],4)<=1.01),~mask_naninf);

            % check signal similarity
            % dwi_2D          = utils.reshape_ND2GD(dwi,mask_valid);
            dwi_2D          = utils.reshape_ND2AD(dwi,mask_valid);
            signalTemplate  = mean(dwi_2D,2,"omitmissing");         % assuming the majority of the signal are from tissues
            signalTemplate  = (signalTemplate - mean(signalTemplate)) ./ std(signalTemplate);
            Rcorr           = zeros(1,size(dwi_2D,2));
            for k = 1:size(dwi_2D,2)
                signalVoxel = dwi_2D(:,k);
                signalVoxel = (signalVoxel - mean(signalVoxel)) ./ std(signalVoxel);
                Rcorr(k)    = corr(signalTemplate,signalVoxel);
            end
            % Rcorr            = utils.reshape_GD2ND(Rcorr,mask_valid);
            Rcorr            = utils.reshape_AD2ND(Rcorr,mask_valid);
            mask_dissimilar  = Rcorr < 0.1;
            mask_valid       = and(mask_valid,~mask_dissimilar);

            if numel(mask_valid(mask_valid)) ~= numel(mask(mask)) 
                disp('Signal mask is updated! Please use the fitted mask in your subsequent analysis.');
                mask = mask_valid;
            end
        end

        %%%%% Prior estimation related functions %%%%%

        % determine how the starting points will be set up
        function x0 = determine_x0(this,y,mask,mask_CSF,fitting) 

            disp('---------------');
            disp('Starting points');
            disp('---------------');

            dims = size(mask,1:3);

            if ischar(fitting.start)
                switch lower(fitting.start)
                    case 'likelihood'
                        % using maximum likelihood method to estimate starting points
                        % x0 = this.estimate_prior(y,mask,[],fitting.lmax);
                        x0 = this.estimate_prior(y,mask,mask_CSF,fitting.Ntrain,fitting.lmax);
    
                    case 'default'
                        % use fixed points
                        fprintf('Using default starting points for all voxels at [%s]: [%s]\n', cell2str(this.modelParams),replace(num2str(this.startPoint(:).',' %.2f'),' ',','));
                        x0 = utils.initialise_x0(dims,this.modelParams,this.startPoint);

                end
            else
                % user defined starting point
                x0 = fitting.start(:);
                fprintf('Using user-defined starting points for all voxels at [%s]: [%s]\n',cell2str(this.modelParams),replace(num2str(x0(:).',' %.2f'),' ',','));
                x0 = utils.initialise_x0(dims,this.modelParams,this.startPoint);

            end
            
            % make sure the input is bounded
            x0 = askadam.set_boundary(x0,fitting.ub,fitting.lb);

            fprintf('Estimation lower bound [%s]: [%s]\n',      cell2str(this.modelParams),replace(num2str(fitting.lb(:).',' %.2f'),' ',','));
            fprintf('Estimation upper bound [%s]: [%s]\n',      cell2str(this.modelParams),replace(num2str(fitting.ub(:).',' %.2f'),'  ',','));
            disp('---------------');
        end

        % using maximum likelihood method to estimate starting points
        function pars0 = estimate_prior(this, dwi, mask, mask_CSF, Nsample, lmax)
        % Estimation starting points for NEXI using likehood method

            start = tic;
            
            disp('Estimate starting points based on likelihood ...')

            % manage pool
            pool            = gcp('nocreate');
            isDeletepool    = false;
            if numel(mask(mask>0)) > 1e4    % only start a pool if many voxel
                if isempty(pool)
                    Nworker = min(max(8,floor(maxNumCompThreads/4)),maxNumCompThreads);
                    pool    = parpool('Processes',Nworker);
                    isDeletepool = true;
                end
            end

            if nargin < 4 || isempty(Nsample)
                Nsample         = 1e4;
                % Nsample         = 1;
            end
            % create training data
            [x_train, S_train] = this.traindata(Nsample,lmax);

            % reshape input data,  put DWI dimension to 1st dim
            dims    = size(dwi);
            dwi     = permute(dwi,[4 1 2 3]);
            dwi     = reshape(dwi,[dims(4), prod(dims(1:3))]);

            % find masked voxels
            ind         = find(mask(:));
            if lmax == 0
                Nparam = 5;
            elseif lmax == 2
                Nparam = 6;
            end

            % pars0_mask  = zeros(Nparam,length(ind));
            pars0_mask = this.likelihood(dwi(:,ind), x_train, S_train,lmax);
            % if ~isempty(pool)
            %     parfor kvol = 1:length(ind)
            %         pars0_mask(:,kvol) = this.likelihood(dwi(:,ind(kvol)), x_train, S_train,lmax);
            %     end
            % else
                % for kvol = 1:length(ind)
                %     pars0_mask(:,kvol) = this.likelihood(dwi(:,ind(kvol)), x_train, S_train,lmax);
                % end
            % end
            pars           = zeros(Nparam,size(dwi,2));
            pars(:,ind)    = pars0_mask;

            % reshape estimation into image
            pars           = permute(reshape(pars,[size(pars,1) dims(1:3)]),[2 3 4 1]);

            % Correction for CSF
            % bval_thres      = max(min(gather(this.b)),1.1);
            % idx             = gather(this.b) <= bval_thres;
            % D0              = real(this.b(idx)\-log(dwi(cat(1,idx,false(size(idx))),:)));
            % D0              = permute(reshape(D0,[size(D0,1) dims(1:3)]),[2 3 4 1]);
            % D0              = max(utils.set_nan_inf_zero(D0),0);
            % mask_CSF        = D0>1.5;
            
            % % ratio to modulate pars0 estimattion
            % pars0_csf = [0.01,0.01,1,1,1,0.01];
            % % pars0_csf = [0.01,1,1,0.01,0.01];
            % for k = 1:size(pars,4)
            %     tmp                 = pars(:,:,:,k);
            %     tmp(mask_CSF==1)    = tmp(mask_CSF==1).*pars0_csf(k);
            %     pars(:,:,:,k)      = tmp;
            % end

            ET  = duration(0,0,toc(start),'Format','hh:mm:ss');
            fprintf('Starting points estimated. Elapsed time (hh:mm:ss): %s \n',string(ET));
            if isDeletepool
                delete(pool);
            end

            for km = 1:size(pars,4)
                pars0.(this.modelParams{km}) = pars(:,:,:,km);
            end

        end

        % create training data for likelihood
        function [x_train, S_train, intervals] = traindata(this, modelFWD, N_samples, lmax, varargin)
            % intervals = [this.lb, this.ub];
            intervals = [ 1 200   ;   % exchange time = (1-fn)/r
                          0.1 1  ;   % fn
                          0.1 1     ;   % 1/aD
                          5 15   ;   % aC
                          0.01 0.5;  % Kinf
                             ];   
            if nargin >4
                input = varargin{1};
            else
                input = [];
            end
            
            numBSample = size(this.t_vec,1);
            numParam   = 5;%size(this.lb,1);

            % % Use Latin Hypercube Sampling (LHS)
            % P = lhsdesign(N_samples, numParam);  % samples in [0, 1]
            % dictionary = zeros(N_samples, numParam);
            % for i = 1:numParam
            %     dictionary(:,i) = this.lb(i) + P(:,i)*(this.ub(i) - this.lb(i));
            % end           

            % batch size can be modified according to available hardware
            batch_size  = 1e3;

            divisors = find(mod(N_samples, 1:N_samples) == 0);
            [~, idx] = min(abs(divisors - batch_size));
            batch_size = divisors(idx);
            % batch_size  = 1;
            reps        = floor(N_samples/batch_size);
            x_train     = zeros(numParam,batch_size,reps);
            S_train     = zeros(numBSample,batch_size,reps);
            parfor k = 1:reps
                % generate random parameter guesses and construct batch for NN signal evaluation
                % pars0 = dictionary(((k-1)*batch_size+1:k*batch_size),:)';
                % pars0 = this.lb(1:numParam) + (this.ub(1:numParam) - this.lb(1:numParam)).*rand(numParam,batch_size);

                % generate random parameter guesses and construct batch for NN signal evaluation
                pars0 = intervals(:,1) + diff(intervals,[],2).*rand(size(intervals,1),batch_size);
                % pars0(1,:) = (1-pars0(2,:))./pars0(1,:);
                pars0 = pars0(1:numParam,:);
                pars = struct;
                for km = 1:size(pars0,1)
                    pars.(this.modelParams{km}) = pars0(km,:);
                end

                % Karger signal evaluation using RK4
                if isempty(input)
                    S = modelFWD(pars);
                else
                    S = modelFWD(pars,input);
                end

                % remaining signals (dot, soma)
                x_train(:,:,k) = pars0;
                S_train(:,:,k) = S;

            end

            if lmax == 2
                intervals(6,:) = [];
            end
        end
        
        % likelihood
        function [pars_best, s_best, sse_best] = likelihood(this, S0, x_train, S_train, lmax)
            wt = kron(this.Nav(:), 1./(2*(0:2:lmax)+1));
            % wt = kron(this.Nav(:), [1,1]);
            wt = wt(:);
            nL = floor(lmax/2);
            S0 = S0(1:size(this.t_vec,1)*(nL+1),:);
            % batch size can be modified according to available hardware
            Nx = size(x_train,1);
            Ns = size(S_train,1);
            [~, Nv] = size(S0);
            pars_best = zeros(Nx,Nv);
            s_best = zeros(Ns,Nv);
            sse_best  = inf(Nv,1);
            % for k = 1:reps
            % pars = x_train(:,:,k);
            % S    = S_train(:,:,k);
            pars = reshape(x_train,Nx,[]);
            S    = reshape(S_train,Ns,[]);
            parfor i = 1:Nv
                S0i = S0(:,i);

                % scale generated signals (fit S0) to input signal
                sse = sum(wt.*(S0i - (S0i'*S)./dot(S,S).*S).^2);
                % mae = mean(wt.*abs(S0i - (S0i'*S)./dot(S,S).*S));

                % store best encountered parameter combination
                % [sse_new,best_index] = min(sse);
                [sse_new,best_index] = min(sse);
                % if sse_new<sse_best(i)
                sse_best(i)    = sse_new;
                pars_best(:,i) = pars(:,best_index);
                s_best(:,i)    = S(:,best_index);
                % end
            end
            % end
        end

        

        %% OGSE diffusivity and kurtosis related functions

        % Forward model to generate kurtosis signal
        function Kc = FWD(this, pars, lmax)
        
            % Model parameters
            % rn          = pars.rn;
            tex          = pars.tex;
            fn           = max(pars.fn, askadam.epsilon); % avoid division by zeros when computing re
            aD           = 1./pars.aD;
            % aD           = pars.aD;
            aC           = pars.aC;
            Kinf         = pars.Kinf;
            % Input parameters
            if isgpuarray(fn)
                t             = gpuArray(single(this.t_vec));
                gt            = gpuArray(single(this.gt_vec));
                qw2int        = dlarray(gpuArray(single(this.qw2_int)));
            else
                t             = single(this.t_vec);
                gt            = single(this.gt_vec);
                qw2int       = single(this.qw2_int);
            end

            % Calculate Kvar
            Kvar = this.KvarCosTrap(qw2int, fn, aD, aC);

            % Calculate h[q]
            rex = 1./tex;
            % rex = rn./(1-fn);
            hq = this.KMhq(t, gt, rex);
            if isgpuarray(rex)
                hq = dlarray(reshape(hq,[size(Kvar,1),size(Kvar,2),size(Kvar,3)]),'SSS');
            else
                hq = reshape(hq,[size(Kvar,1),size(Kvar,2),size(Kvar,3)]);
            end
            
            % Calulate Karger model derived kurtosis, K_KM
            Kc = Kvar.*hq + Kinf;
        end
    end  
    methods(Static)
        
        function K = KvarCosTrap(qw2_int, f, aD, aC)
            Nx = 100;
            x  = zeros([ones(1,ndims(aD)), Nx],'single','gpuArray'); x(:) = linspace(0,0.99,Nx);
            typedl = repmat('S',1,ndims(x));
            if isgpuarray(f)
                xx = dlarray(x,typedl);
                Kvar = dlarray(3*f.*(1-f).*( (aD+aC.*qw2_int).*xx.^2./(1-xx.^2) + (1-f) ).^(-2),typedl);
                K = dlarray(trapz(x(:), Kvar, ndims(x)),typedl);
            else
                Kvar = 3*f.*(1-f).*( (aD+aC.*qw2_int).*x.^2./(1-x.^2) + (1-f) ).^(-2);
                K = trapz(x(:), Kvar, ndims(x));
            end      
        end
        

        function hq = KMhq(t, gt, rex)

            [Nf, Nt] = size(t);
            dt = t(:,2)-t(:,1);
            qt = cumsum(gt, 2).*dt;
            b  = sum(qt.^2, 2).*dt;

            q4t = zeros(Nf, Nt);
            for fi = 1:Nf
                q4 = @(x) 1./b(fi).^2*dt(fi) .* sum( qt(fi,1:end-floor(x/dt(fi))).^2 .* qt(fi,floor(x/dt(fi))+1:end).^2, 2);           
                for i = 1:Nt
                    q4t(fi,i) = q4(t(fi,i));
                end
            end

            % Ensure rex is [1, 1, 5000] for broadcasting
            if isgpuarray(rex) 
                rex = dlarray(reshape(rex, [1, 1, numel(rex)]),'SSS'); % Now rex is [1,1,5000]
            else
                rex = reshape(rex, [1, 1, numel(rex)]);
            end
            % Ensure t is [18, 1000, 1] for broadcasting
            if isgpuarray(rex)
                 t = dlarray(reshape(t, [size(t,1), size(t,2), 1]),'SSS'); % Now t is [18,1000,1]
            else
                t = reshape(t, [size(t,1), size(t,2), 1]); % Now t is [18,1000,1]
            end
            
            % Compute the exponential term with automatic broadcasting
            exp_term = exp(-rex .* t); % Now size(exp_term) = [18, 1000, 5000]
            
            % Ensure q4t is compatible by expanding it to match third dimension
            if isgpuarray(rex)
                q4t_exp = dlarray(repmat(q4t, [1, 1, size(rex,3)]),'SSS'); % Now q4t_exp is [18,1000,5000]
            else
                q4t_exp = repmat(q4t, [1, 1, size(rex,3)]); % Now q4t_exp is [18,1000,5000]
            end
            
            % Compute the summation along the second dimension
            sum_term = sum(exp_term .* q4t_exp, 2); % Summing along dimension 2 -> size is [18, 1, 5000]
            
            % Multiply by 2 * dt, ensuring dt is properly broadcasted
            hq = 2 * dt .* squeeze(sum_term); % Resulting size should be [18, 5000]
        end
       
       
        %% Utilities
        % check and set default fitting algorithm parameters
        function fitting2 = check_set_default(fitting)
            % get basic fitting setting check
            fitting2 = askadam.check_set_default_basic(fitting);

            % get customised fitting setting check
            if ~isfield(fitting,'regmap');          fitting2.regmap     = 'fn';             end
            if ~isfield(fitting,'lmax');            fitting2.lmax       = 0;                end
            if ~isfield(fitting,'start');           fitting2.start      = 'likelihood';     end

            if ~iscell(fitting2.regmap)
                fitting2.regmap = cellstr(fitting2.regmap);
            end

        end

        

    end

end