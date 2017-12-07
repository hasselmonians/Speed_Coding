clear
close all

%% create poisson spike train
% Parameters
fr_mean = 15;
ns = 1000;
isi = Speed_coding_project.poisson_spiketrain_function(fr_mean,ns);

% get spike train in time domain
train_mat = cumsum(isi)./1000; % time in s
train = mat2cell(train_mat,1);

% plot raster plot of spike train
figure
plotSpikeRaster(train,'PlotType','vertline')

%% get spiking rate of spike train
spkRate = Speed_coding_project.get_spkRate(0.02,train_mat);

% plot spiking rate over time
figure
plot(spkRate.t,spkRate.ns,'k')
ylabel('Spiking rate (Hz)')
xlabel('Time (s)')

figure
plot(spkRate.t,spkRate.z,'g')
ylabel('z-score of firing rate')
xlabel('Time (s)')
clear edges edges_corrected N zScore
