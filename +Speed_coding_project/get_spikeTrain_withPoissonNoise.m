function [art_noisy_spkTrain] = get_spikeTrain_withPoissonNoise(root,cel,speedMod)
% creates a spike train with perfect correlation to the speed signal and
% adds noise randomly drawn from a Poisson distribution to it


[~,indx] = ismember(cel,root.cells,'rows'); % find row index of cel in root.cells and speedMod data

% Parameters
slope = speedMod.Non.b(indx); % pull out the slope as shown in the speed vs. firing rate plot in root.Visualize2 from speedMod structure; focus on baseline (NonInhibition) epochs
y_int = speedMod.Non.y_int(indx); % pull out y_int for baseline (NonInhibition) epochs

speed = root.svel;
art_spkTrain = speed.*slope + y_int; % creates the perfectly correlated artificial spike train

% introduce Poisson noise to spike train
art_noisy_spkTrain = zeros(1,length(art_spkTrain)); %initialize vector
for i = 1:length(art_noisy_spkTrain)
    art_noisy_spkTrain(i) = poissrnd(art_spkTrain(i));
end



