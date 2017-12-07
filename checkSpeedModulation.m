function [outTime, correl] = checkSpeedModulation(self, windowWidth, shift)

zSpeed = get_zScore_speed(self);
spkRate = get_spkRate(self);
zSpkRate = spkRate{1}.z;
time = spkRate{1}.t(2:end);

Fs = 1 / (self.ts(2) - self.ts(1));

beginInd = 1;
endInd = beginInd + round(windowWidth*Fs);

count = 1;
while endInd <= length(time)
%     if count == 4;
%         count
%     end
    correl(count) = corr(zSpkRate(beginInd:endInd)', zSpeed(beginInd:endInd), 'type', 'Pearson');
%     correl(count) = corr(zSpkRate(beginInd:endInd)', self.svel(beginInd:endInd), 'type', 'Pearson');
    outTime(count) = mean(time(beginInd:endInd));
    beginInd = beginInd + round(shift*Fs);
    endInd = beginInd + round(windowWidth*Fs);
    count = count + 1;
end

figure()
plot(outTime, correl)
xlabel('Time (s)')
ylabel('Pearson Correlation Coefficient')

end