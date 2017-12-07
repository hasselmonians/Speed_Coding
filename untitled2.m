close all
clear
% load('/media/craig_kelley/Datastore/Dropbox (hasselmonians)/hdannenb (1)/UnitRecordingData/ChATControl_5/170406_S1_LS.mat')
% load('/media/craig_kelley/Datastore/Dropbox (hasselmonians)/hdannenb (1)/UnitRecordingData/ArchTChAT_14/170814_S1_LS.mat')
% load('/media/craig_kelley/Datastore/Dropbox (hasselmonians)/hdannenb (1)/UnitRecordingData/PVArchT_2/170926_S1_LS.mat')
load('/media/craig_kelley/Datastore/Dropbox (hasselmonians)/hdannenb (1)/UnitRecordingData/PVArchT_2/170925_S1_LS.mat')
root2 = root;
root2.cel = [1,1];

[time, correl] = checkSpeedModulation(root2, .5, .1);

ind = PeriodsAboveThr(correl,.3,5);
first = 1;
epoch = nan(length(ind),2);

for i = 1:length(ind)
    [~, last] = min(abs(root2.ts - time(ind{i}(1))));
    last = last(1) - 1;
    epoch(i,:) = [root2.ts(first), root2.ts(last)];
    [~, first] = min(abs(root2.ts - time(ind{i}(end))));
    first = first(1) + 1;
end

root2.epoch = epoch;

% root.Visualize
% 
% root2.Visualize