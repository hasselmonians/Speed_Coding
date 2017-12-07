function spkRate = get_spkRate(sampling,spk_times)
% inputs: sampling, sampling rate of spiking rate in seconds; spk_times, spiking times
% (spike train) in seconds

edges = 0:sampling:spk_times(end);
[N,edges] = histcounts(spk_times,edges);
edges_corrected = edges(1:end-1)+diff(edges)/2;
mean_spkRate = mean(N);
std_spkRate = std(N);
zScore = (N-mean_spkRate)./std_spkRate;
% pull out results
spkRate.sampling = sampling;
spkRate.ns = N; % number of spikes per time bin
spkRate.t = edges_corrected;
spkRate.spkRate = N/sampling;
% spkRate.spkRate_smoothed = smooth(spkRate{e}.spkRate,smooth_factor);
spkRate.z = zScore;
spkRate.spkRate_avg = mean(N);
