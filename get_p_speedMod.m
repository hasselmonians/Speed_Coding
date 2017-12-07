function [speedMod] = get_p_speedMod(self, OFF_epochs, ON_epochs, varargin)
% determines if cells in CMBHOME.Session object are signifcantly speed
% modulated based on R. Saves data as speedMod under self.name, if varargin
% name-value pair 'save' is set to 1; exclusion criterium 1a: if spiking
% rate drops below min_SpkRate per individual time epoch, this epoch is
% excluded from analysis. exclusion criterium 1b: if less than 5 epochs are
% present for either the Laser OFF or Laser ON condition, the cell is
% excluded from analysis. exclusion criterium 2: cells with less than
% min_nmbrSpikes spikes in either condition are excluded; exclusion
% criterium 3: if speed bins do not cover the range up to min_speedBin or
% if the speed ranges differ more than twice between conditions, cells are
% excluded.

% default parameters
default_alpha_level = 1;% define significance level for speed tuning
default_min_spkRate = 0.1;% minimum spiking rate per inhbition/non-inhibition time epochs
default_min_nmbrSpikes = 150;% minimum number of spikes in total sum of all inhibition/non-inhibition
% time epochs
default_min_nmbr_epochs = 5;% define minimal number of inhibition and non-inhibition epochs with
% spiking rate higher than min_spkRate
default_min_speedBin = 12; % default minimal maximum speed Bin
% decide if speedMod is saved in root object containing file
default_save = 0; % default: do not save
% select cells to analyze
default_cells_to_analyze = 1:size(self.cells,1); % default: run script over all cells

p = inputParser;
addParameter(p,'alpha_level',default_alpha_level);
addParameter(p,'min_spkRate',default_min_spkRate);
addParameter(p,'min_nmbrSpikes',default_min_nmbrSpikes);
addParameter(p,'min_nmbr_epochs',default_min_nmbr_epochs);
addParameter(p,'min_speedBin',default_min_speedBin);
addParameter(p,'save',default_save);
addParameter(p,'cells_to_analyze',default_cells_to_analyze);
parse(p,varargin{:})
alpha_level = p.Results.alpha_level;
min_spkRate = p.Results.min_spkRate;
min_nmbrSpikes = p.Results.min_nmbrSpikes;
min_nmbr_epochs = p.Results.min_nmbr_epochs;
min_speedBin = p.Results.min_speedBin;
save_speedMod = p.Results.save;
cells_to_analyze = p.Results.cells_to_analyze;

import CMBHOME.Utils.*

% initialize variables
speedMod.Non.p = zeros(size(self.cells,1),1);
speedMod.Inh.p = zeros(size(self.cells,1),1);
speedMod.Non.sig_R = zeros(size(self.cells,1),1);
speedMod.Inh.sig_R = zeros(size(self.cells,1),1);
speedMod.Non.R = zeros(size(self.cells,1),1);
speedMod.Inh.R = zeros(size(self.cells,1),1);
speedMod.Non.b = zeros(size(self.cells,1),1);
speedMod.Inh.b = zeros(size(self.cells,1),1);
speedMod.Non.rsq = zeros(size(self.cells,1),1);
speedMod.Inh.rsq = zeros(size(self.cells,1),1);
speedMod.Non.y_int = zeros(size(self.cells,1),1);
speedMod.Inh.y_int = zeros(size(self.cells,1),1);
speedMod.Non.N_spikes = zeros(size(self.cells,1),1);
speedMod.Inh.N_spikes = zeros(size(self.cells,1),1);
speedMod.IncludeForSpeedModChangeAna = ones(size(self.cells,1),1); % default value is to include for analysis of change in speed tuning
speedMod.Non.epochs = cell(size(self.cells,1),1);
speedMod.Inh.epochs = cell(size(self.cells,1),1);

speedMod.sig_level = alpha_level;

%% order epochs in time. Otherwise root.ts will be empty for all epochs which occur earlier in time than a previous epoch
OFF_epochs = sortrows(OFF_epochs);
ON_epochs = sortrows(ON_epochs);

%% for baseline conditions
self.epoch = OFF_epochs;
% define speed criterium
vel_dim_Non = 2*self.spatial_scale^-1:self.spatial_scale^-1:prctile(CMBHOME.Utils.ContinuizeEpochs(self.vel),95);

%% exclusion criterium 1
for k = cells_to_analyze
    self.cel = [self.cells(k,1),self.cells(k,2)];
    speedMod.Non.epochs{k} = self.epoch;
    counter = 0;
    for j = 1:size(speedMod.Non.epochs{k},1)
        i = j-counter;
        time_epoch = (speedMod.Non.epochs{k}(i,2) - speedMod.Non.epochs{k}(i,1));
        if size(self.cel_ts{j,1},1) < time_epoch*min_spkRate % exclusion criterium 1a: min_spkRate Hz as minimum spiking rate per time bin
            speedMod.Non.epochs{k}(i,:) = [];
            counter = counter + 1;
        end
    end
    if size(speedMod.Non.epochs{k},1) < min_nmbr_epochs % exclusion criterium 1b
        speedMod.IncludeForSpeedModChangeAna(k,1) = 0; % exclude (set value to 0) if exclusion criteria is met
    end
end

for k = cells_to_analyze
    if ~isempty(speedMod.Non.epochs{k})
        self.epoch = speedMod.Non.epochs{k};
        [~,~,N_spikes] = self.VelocityRate_HD([self.cells(k,1),self.cells(k,2)], vel_dim_Non);
        if N_spikes > min_nmbrSpikes % exclusion criterium 2
            [p,sig_R,alpha_level,~,R,rsq,b,y_int] = shuffle_speedBins(self, self.cells(k,:), alpha_level);
            speedMod.Non.p(k,1) = p;
            speedMod.Non.sig_R(k,1) = sig_R;
            speedMod.Non.b(k,1) = b;
            speedMod.Non.rsq(k,1) = rsq;
            speedMod.Non.y_int(k,1) = y_int;
            speedMod.Non.R(k,1) = R;
            speedMod.Non.N_spikes(k,1) = N_spikes;
            clear p sig_R b rsq y_int R N_spikes
        else
            speedMod.Non.p(k,1) = nan;
            speedMod.Non.sig_R(k,1) = nan;
            speedMod.Non.b(k,1) = nan;
            speedMod.Non.rsq(k,1) = nan;
            speedMod.Non.y_int(k,1) = nan;
            speedMod.Non.R(k,1) = nan;
            speedMod.Non.N_spikes(k,1) = N_spikes;
            speedMod.IncludeForSpeedModChangeAna(k,1) = 0;
            clear N_spikes
        end
    else
        speedMod.Non.p(k,1) = nan;
        speedMod.Non.sig_R(k,1) = nan;
        speedMod.Non.b(k,1) = nan;
        speedMod.Non.rsq(k,1) = nan;
        speedMod.Non.y_int(k,1) = nan;
        speedMod.Non.R(k,1) = nan;
        speedMod.Non.N_spikes(k,1) = nan;
        speedMod.IncludeForSpeedModChangeAna(k,1) = 0;
        clear N_spikes
    end
end

%% for inhibition condition
self.epoch = ON_epochs;
% define speed criterium
vel_dim_Inh = 2*self.spatial_scale^-1:self.spatial_scale^-1:prctile(CMBHOME.Utils.ContinuizeEpochs(self.vel),95);

%% exclusion criterium 1
for k = cells_to_analyze
    self.cel = [self.cells(k,1),self.cells(k,2)];
    speedMod.Inh.epochs{k} = self.epoch;
    counter = 0;
    for j = 1:size(speedMod.Inh.epochs{k},1)
        i = j-counter;
        time_epoch = (speedMod.Inh.epochs{k}(i,2) - speedMod.Inh.epochs{k}(i,1));
        if size(self.cel_ts{j,1},1) < time_epoch*min_spkRate % exclusion criterium 1a: min_spkRate Hz as minimum spiking rate per time bin
            speedMod.Inh.epochs{k}(i,:) = [];
            counter = counter + 1;
        end
    end
    if size(speedMod.Inh.epochs{k},1) < min_nmbr_epochs % exclusion criterium 1b
        speedMod.IncludeForSpeedModChangeAna(k,1) = 0; % exclude (set value to 0) if exclusion criteria is met
    end
end

for k = cells_to_analyze
    if ~isempty(speedMod.Inh.epochs{k})
        self.epoch = speedMod.Inh.epochs{k};
        [~,~,N_spikes] = self.VelocityRate_HD([self.cells(k,1),self.cells(k,2)], vel_dim_Inh);
        if N_spikes > min_nmbrSpikes % exclusion criterium 2
            [p,sig_R,alpha_level,~,R,rsq,b,y_int] = shuffle_speedBins(self, self.cells(k,:), alpha_level); % get speed tuning curve (linear regression), and do shuffle to determine significance of firing rate-speed correlation
            speedMod.Inh.p(k,1) = p;
            speedMod.Inh.sig_R(k,1) = sig_R;
            speedMod.Inh.b(k,1) = b;
            speedMod.Inh.rsq(k,1) = rsq;
            speedMod.Inh.y_int(k,1) = y_int;
            speedMod.Inh.R(k,1) = R;
            speedMod.Inh.N_spikes(k,1) = N_spikes;
            clear p sig_R b rsq y_int R N_spikes
        else
            speedMod.Inh.p(k,1) = nan;
            speedMod.Inh.sig_R(k,1) = nan;
            speedMod.Inh.b(k,1) = nan;
            speedMod.Inh.rsq(k,1) = nan;
            speedMod.Inh.y_int(k,1) = nan;
            speedMod.Inh.R(k,1) = nan;
            speedMod.Inh.N_spikes(k,1) = N_spikes;
            speedMod.IncludeForSpeedModChangeAna(k,1) = 0;
            clear N_spikes
        end
    else
        speedMod.Inh.p(k,1) = nan;
        speedMod.Inh.sig_R(k,1) = nan;
        speedMod.Inh.b(k,1) = nan;
        speedMod.Inh.rsq(k,1) = nan;
        speedMod.Inh.y_int(k,1) = nan;
        speedMod.Inh.R(k,1) = nan;
        speedMod.Inh.N_spikes(k,1) = nan;
        speedMod.IncludeForSpeedModChangeAna(k,1) = 0;
        clear N_spikes
    end
end

%% apply exclusion criterium 3: if speed bins do not cover the range up
% to min_speedBin cm/s or if the speed range differs more than twice between
% conditions, cells are excluded.
if vel_dim_Non(end) < min_speedBin*self.spatial_scale^-1 || vel_dim_Inh(end) < min_speedBin*self.spatial_scale^-1 || vel_dim_Non(end)/vel_dim_Inh(end) > 1.5 || vel_dim_Non(end)/vel_dim_Inh(end) < 0.667
    for i = 1:length(speedMod.Non.sig_R)
        speedMod.IncludeForSpeedModChangeAna(i,1) = 0; % exclude (set value to 0) if exclusion criteria is met
    end
end
%% apply exclusion criterium 4: if spiking rates in first and second half of session differs more than twice, cells are excluded
% self.epoch = [-Inf Inf];
% half_time = self.ts(end)/2;
% self.epoch = [-Inf half_time; half_time Inf];
% for i = cells_to_analyze
%     self.cel = self.cells(i,:);
%     test = length(self.cel_ts{1})/length(self.cel_ts{2});
%     if test > 2 || test < 0.5
%         speedMod.IncludeForSpeedModChangeAna(i,1) = 0; % exclude (set value to 0) if exclusion criteria is met
%     end
% end

if save_speedMod
    [~,savename,~] = fileparts(self.name);
    save(savename,'speedMod','-append')
    clear sig_level
end

    function [p,sig,sig_level,SSresid,R,rsq,b,y_int] = shuffle_speedBins(self, cel, sig_level)
        %shuffle_speedBins shuffles the speed bins, fits a linear regression and
        %determines if the cell is significantly speed modulated with significance
        %level sig_level (in per cent) based on the distribution of R values of
        %the shuffled data
        
        if ~exist('sig_level','var') == 1
            sig_level = 1; %default value for sig_level in percent
        end
        
        import CMBHOME.Utils.*
        
        vel_dim = 2*self.spatial_scale^-1:self.spatial_scale^-1:prctile(CMBHOME.Utils.ContinuizeEpochs(self.vel),95);
        
        [F, V, ~] = self.VelocityRate_HD(cel, vel_dim);
        BinHalfWidth = diff(V)/2; %take the middle of each bin
        V(end) = [];
        V = V + BinHalfWidth;
        V = V(~isnan(F));
        F = F(~isnan(F));
        % apply root.spatial_scale to velocity
        V = V*self.spatial_scale;
        
        % fit a regression line to the mean values of each speed bin
        x = [V' ones(length(V),1)];
        Coeffs = x\F';
        b = Coeffs(1,1);
        y_int = Coeffs(2,1);
        [R, ~] = corrcoef(V,F);
        R = R(1,2);
        LinReg = @(x) b*x + y_int; %  Linear fit equation
        Fn = LinReg(V); % values predicted by linear fit
        yresid = F - Fn; % residuals
        SSresid = sum(yresid.^2); %sum of squared residuals
        SStotal = (length(F)-1) * var(F); %total sum of squares (variance)
        rsq = 1 - SSresid/SStotal; % coefficient of determination (equals R^2 for linear regression)
        
        % shuffling
        max_cycles = 1000; %number of permutations
        random = zeros(max_cycles,length(F));
        for kk = 1:max_cycles
            random(kk,:) = randperm(length(F));
        end
        random = unique(random,'rows');
        cycles = size(random,1);
        Coeffs_s = cell(1,cycles);
        b_s = cell(1,cycles);
        R_s = cell(1,cycles);
        P_s = cell(1,cycles);
        y_int_s = cell(1,cycles);
        F_s = cell(1,cycles);
        for ii = 1:cycles
            F_s{ii} = F(random(ii,:));
            Coeffs_s{ii} = x\F_s{ii}';
            b_s{ii} = Coeffs_s(1,1);
            y_int_s{ii} = Coeffs_s{ii}(2,1);
            [R_s{ii}, P_s{ii}] = corrcoef(V,F_s{ii});
            R_s{ii} = R_s{ii}(1,2);
            P_s{ii} = P_s{ii}(1,2);
        end
        
        R_shuffle = cell2mat(R_s);
        
        prc_left = prctile(R_shuffle,sig_level/2); % for finding positively speed modulated cells
        prc_right = prctile(R_shuffle,100-sig_level/2); % for finding negatively speed modulated cells
        if R < prc_left || R > prc_right
            sig = 1;
        else
            sig = 0;
        end
        
        % get exact p-value
        p = get_percentrank(R,R_shuffle,cycles);
    end

    function [p] = get_percentrank(value,distribution,cycles)
        percentiles = tiedrank(distribution)/length(distribution);
        [sorted_percentiles,I] = sort(percentiles);
        distribution = distribution(I); % sort distribution in rank order
        rank = length(distribution(distribution <= value));
        if rank == 0
            p = sorted_percentiles(1)/2;
        elseif rank == cycles
            p = (sorted_percentiles(cycles)+1)/2;
        elseif value == distribution(rank)
            p = sorted_percentiles(rank);
        else
            p = (sorted_percentiles(rank) + sorted_percentiles(rank + 1))/2;
        end
    end
end


