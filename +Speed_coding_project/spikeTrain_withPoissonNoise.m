% creates a spike train with perfect correlation to the speed signal and
% adds noise randomly drawn from a Poisson distribution to it

% Parameters
slope = 0.353 % enter the slope as shown in the speed vs. firing rate plot in root.Visualize2

speed = root.svel;
art_spkTrain = speed.*slope; % creates the perfectly correlated artificial spike train

% introduce Poisson noise to spike train
art_noisy_spkTrain = zeros(1,length(art_spkTrain)); %initialize vector
for i = 1:length(art_noisy_spkTrain)
    art_noisy_spkTrain(i) = poissrnd(art_spkTrain(i));
end



