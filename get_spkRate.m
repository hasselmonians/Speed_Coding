function [spkRate] = get_spkRate(root,varargin)
% plot_spkRate_over_time(root) plots the spiking rate over time for
% a the cell set in the root object.

p = inputParser;
addParameter(p,'plot',1)
addParameter(p,'resolution',0.02) % resolution in seconds
parse(p,varargin{:})
plot_fig = p.Results.plot;
resolution = p.Results.resolution;

% [~,name,~] = fileparts(root.name);
% [~,parent_folder,~] = fileparts(pwd);

for e = 1:length(root.cel_ts) % loop over epochs
    if ~isempty(root.cel_ts{e}) % if no spikes are present, root.cel_ts is empty
        edges = 0:resolution:root.cel_ts{e}(end);
        [N,edges] = histcounts(root.cel_ts{e},edges);
        edges_corrected = edges(1:end-1)+diff(edges)/2;
    else
        N = 0;
        edges_corrected = [];
    end
    % calculate z-score of firing rate
    mean_spkRate = mean(N);
    std_spkRate = std(N);
    zScore = (N-mean_spkRate)./std_spkRate;
    % pull out results
    spkRate{e}.resolution = resolution;
%     spkRate{e}.t = edges_corrected;
    spkRate{e}.t = edges;
    spkRate{e}.spkRate = N/resolution;
%     spkRate{e}.spkRate_smoothed = smooth(spkRate{e}.spkRate,smooth_factor);
    spkRate{e}.z = zScore;
    spkRate{e}.spkRate_avg = mean(N);
    
    % figure
%     if plot_fig
%         figure
%         plot(edges_corrected,N,'k')
%         ylabel('Spiking rate (Hz)')
%         xlabel('Time (s)')
%         title(strcat(parent_folder,',',name,sprintf(', T%dC%d',root.cel(1),root.cel(2))),'Interpreter','none')
%         
%         figure
%         plot(edges_corrected,zScore,'g')
%         ylabel('z-score of firing rate')
%         xlabel('Time (s)')
%         title(strcat(parent_folder,',',name,sprintf(', T%dC%d',root.cel(1),root.cel(2))),'Interpreter','none')
%     end
    clear edges edges_corrected N zScore
end

