function [F,V,N_spikes,p,SSresid,R,rsq,b,y_int,F_total,V_total,N_spikes_total] = get_shortTimeScale_speedTuning(root,epochs,cel,time_window)
% gets speed tuning curves for subsequent time windows within epochs;
% Inputs: root, root object; epochs, epochs for analysis as matrix(e.g.
% Laser OFF (baseline) epochs in seconds in the form [5 120; 270 385] for
% two epochs); cell, cell as 1-by-1 matrix ([T C]); time_window, time
% window in seconds

import CMBHOME.Utils.*

% define speed range
min_speed = 2; % minimum speed in cm/s
max_speed = prctile(CMBHOME.Utils.ContinuizeEpochs(root.vel),95); % maximum speed defined as the 95% percentile of the speed signal distribution
speed_bin_size = 1; % size of speed bins in cm/s
speed_edges = min_speed*root.spatial_scale^-1:speed_bin_size*root.spatial_scale^-1:max_speed;

% get values for each time window within epochs
windows = cell(size(epochs,1),1);
for w = 1:size(epochs,1)
    counter = 1;
    while epochs(w,1) + counter*time_window < epochs(w,2)
        windows{w}(counter,:) = [epochs(w,1) + (counter-1)*time_window, epochs(w,1) + counter*time_window];
        counter = counter + 1;
    end
end

windows_mat = cell2mat(windows);

% set epochs in root object
root.epoch = windows_mat;

for w = 1:size(windows_mat,1) % for loop over windows
    root.epoch = windows_mat(w,:); % set epoch in root object
    [F{w}, V{w}, N_spikes{w}] = root.VelocityRate_HD(cel, speed_edges);
    BinHalfWidth = diff(V{w})/2; %take the middle of each bin
    V{w}(end) = [];
    V{w} = V{w} + BinHalfWidth;
    V{w} = V{w}(~isnan(F{w})); % take out the speed bins with nan-F values associated with it, because these are non-existing speeds.
    F{w} = F{w}(~isnan(F{w})); % delete nan values also in F to match the speed bins in V
    % apply root.spatial_scale to velocity
    V{w} = V{w}*root.spatial_scale;
    if length(V{w}) > 1 % exclude all time windows with only one speed bin
        % fit a regression line to the mean values of each speed bin
        x{w} = [V{w}' ones(length(V{w}),1)];
        Coeffs{w} = x{w}\F{w}';
        b{w} = Coeffs{w}(1,1);
        y_int{w} = Coeffs{w}(2,1);
        [R{w}, ~] = corrcoef(V{w},F{w});
        R{w} = R{w}(1,2);
        LinReg = @(x) b{w}*x + y_int{w}; %  Linear fit equation
        Fn{w} = LinReg(V{w}); % values predicted by linear fit
        yresid{w} = F{w} - Fn{w}; % residuals
        SSresid{w} = sum(yresid{w}.^2); %sum of squared residuals
        SStotal{w} = (length(F{w})-1) * var(F{w}); %total sum of squares (variance)
        rsq{w} = 1 - SSresid{w}/SStotal{w}; % coefficient of determination (equals R^2 for linear regression)
        % shuffling
        max_cycles = 1000; %number of permutations
        random{w} = zeros(max_cycles,length(F{w}));
        for kk = 1:max_cycles
            random{w}(kk,:) = randperm(length(F{w}));
        end
        random{w} = unique(random{w},'rows');
        cycles{w} = size(random{w},1);
        Coeffs_s{w} = cell(1,cycles{w});
        %     b_s = cell(1,cycles);
        R_s{w} = cell(1,cycles{w});
        %     P_s = cell(1,cycles);
        %     y_int_s = cell(1,cycles);
        F_s{w} = cell(1,cycles{w});
        for ii = 1:cycles{w}
            F_s{w}{ii} = F{w}(random{w}(ii,:));
            Coeffs_s{w}{ii} = x{w}\F_s{w}{ii}';
            %         b_s{ii} = Coeffs_s(1,1);
            %         y_int_s{ii} = Coeffs_s{ii}(2,1);
            [R_s{w}{ii}, P_s{w}{ii}] = corrcoef(V{w},F_s{w}{ii});
            R_s{w}{ii} = R_s{w}{ii}(1,2);
            %         P_s{ii} = P_s{ii}(1,2);
        end
        R_shuffle{w} = cell2mat(R_s{w});
        % get exact p-value
        p{w}= get_percentrank(R{w},R_shuffle{w},cycles{w});
    else
        F{w} = nan;
        V{w} = nan;
        N_spikes{w} = nan;
        p{w} = nan;
        SSresid{w} = nan;
        R{w} = nan;
        rsq{w} = nan;
        b{w} = nan;
        y_int{w} = nan;
    end
end

% for comparison: calculate linear regression for all time windows
% concatenated
% set epochs in root object
root.epoch = windows_mat;
[F_total, V_total, N_spikes_total] = root.VelocityRate_HD(cel, speed_edges);
BinHalfWidth = diff(V_total)/2; %take the middle of each bin
V_total(end) = [];
V_total = V_total + BinHalfWidth;
V_total = V_total(~isnan(F_total)); % take out the speed bins with nan-F values associated with it, because these are non-existing speeds.
F_total = F_total(~isnan(F_total)); % delete nan values also in F to match the speed bins in V
% apply root.spatial_scale to velocity
V_total = V_total*root.spatial_scale;


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


