function isSpeedModulated(filepath)

addpath(genpath('/media/craig_kelley/Datastore/Dropbox (hasselmonians)/hdannenb (1)/toolboxes/CMBHOME/'))
import CMBHOME.*

filein = load(filepath);
self = filein.root;
clear filein
self = self.AppendKalmanVel;

pathname = '/media/craig_kelley/Datastore/Dropbox (hasselmonians)/Craig/Speed Modulation/Checking_Speed_Modulation/Caitlins_Cells/';

if size(self.epoch,1) == 1
    alpha_level = 1;
    min_spkRate = .1;
    min_nmbrSpikes = 150;
    min_speedBin = 12;
    [Ncells, ~] = size(self.cells);
    isSpeedMod = nan(Ncells,6);
    isSpeedMod(:,1) = zeros(Ncells,1);

    vel_dim = 2*self.spatial_scale^-1:self.spatial_scale^-1:prctile(CMBHOME.Utils.ContinuizeEpochs(self.vel),95);
    
    speedBinGood = 1;
    if vel_dim(end) < min_speedBin * self.spatial_scale^-1
        speedBinGood = 0;
    end

    spkRateGood = zeros(Ncells,1);
    numSpikesGood = zeros(Ncells,1);
    for i = 1:Ncells
        self.cel = self.cells(i,:);
        time_epoch = self.epoch(1,2) - self.epoch(1,1);
        if length(self.cel_ts{1}) >= time_epoch*min_spkRate
            spkRateGood(i) = 1;
        end
   
        [~, ~, N_spikes] = self.VelocityRate_HD(self.cells(i,:), vel_dim);
        if N_spikes >= min_nmbrSpikes
            numSpikesGood(i) = 1;
        end
        
        [p, sig_R, ~, ~, R, rsq, b, y_int] = shuffle_speedBins(self, self.cells(i,:), alpha_level);
        
        if spkRateGood(i) && speedBinGood && numSpikesGood(i) && sig_R
            isSpeedMod(i,:) = [1, p, R, rsq, b, y_int];
        end
    end
    
    if sum(isSpeedMod(:,1))
        name = strrep(self.name, '\', '/');
        fileParts = strsplit(name, '/');
        filename = fileParts{end};
        newfilename = strcat(filename(1:end-4),'_speedModCells.mat');
        save(strcat(pathname,newfilename),'isSpeedMod','filepath')
        strout = strcat(filepath, ' contains ', num2str(sum(isSpeedMod(:,1))), ' speed cells.\n');
        fprintf(strout)
    else
        fprintf(strcat(filepath, ' has no speed cells\n'))
    end
            
else
    fprintf(strcat(filepath, ' contains multiple epochs\n'))
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