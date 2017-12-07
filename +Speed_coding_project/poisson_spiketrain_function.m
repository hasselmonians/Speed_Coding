%% Generation of Poisson spike train with refractoriness (modified from "Fundamentals of Computational Neuroscience", Second edition, by Thomas Trappenberg, Table 3.3)
function isi = poisson_spiketrain_function(fr_mean,ns)
% inputs: fr_mean, mean frequency in Hz; ns, number of total spikes (before
% deleting spikes falling into the refractory period)

fr_mean = fr_mean/1000; % mean firing rate conversion from 1/s to 1/ms.

%% generating poisson spike train
lambda = 1/fr_mean; % inverse firing rate
isi1 = -lambda.*log(rand(ns,1)); % generation of expo. distr. ISIs
%% Delete spikes that are within refractory period
is = 0;
for i = 1:ns;
    if rand > exp(-isi1(i)^2/32);
        is = is + 1;
        isi(is) = isi1(i);
    end
end
%% Plotting histogram and calculating cv
figure
hist(isi,50); % Plot histogram of 50 bins
cv = std(isi)/mean(isi) % coefficent of variation


