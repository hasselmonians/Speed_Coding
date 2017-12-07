function [out, speedCenters] = speedPlotInBins(self, windowWidth)

spkRate = get_spkRate(self);
zSpkRate = spkRate{1}.spkRate;
time = spkRate{1}.t(2:end);

speedCenters = 0:2:40;
out = nan(1,length(speedCenters));

Fs = 1 / (self.ts(2) - self.ts(1));

beginInd = 1;
endInd = beginInd + round(windowWidth*Fs);

count = 1;
while endInd <= length(time)
    tempSpeed = self.svel(beginInd:endInd);
    tempFR = zSpkRate(beginInd:endInd);
    
    for i = 1:length(speedCenters)-1
        out(count,i) = nanmean(tempFR(tempSpeed > speedCenters(i) & tempSpeed < speedCenters(i+1)));
    end
    
    out(count,end) = nanmean(tempFR(tempSpeed > speedCenters(end)));
    
    beginInd = endInd + 1;
    endInd = beginInd + round(windowWidth*Fs);
    
    count = count + 1;
    
end

plotOut = nanmean(out);

figure()
plot(speedCenters, plotOut, 'k*')

end