function f = linearSpeedDecoder(self, cells2analyze)

self = self.AppendKalmanVel;
Fs = 1 / (self.ts(2) - self.ts(1));

Inds = PeriodsAboveThr(self.svel, 2, 2);
epoch = zeros(length(Inds),2);
for i = 1:length(Inds)
    epoch(i,1) = self.ts(Inds{i}(1));
    epoch(i,2) = self.ts(Inds{i}(end));
end

[Ncells, ~] = size(cells2analyze);
for i = 1:Ncells
    self.cel = cells2analyze(i,:);
    [~, self] = InstFR(self);
    self.epoch = epoch;
    FR(:,i) = CMBHOME.Utils.ContinuizeEpochs(self.myvar);
end

self.epoch = epoch;
speed = CMBHOME.Utils.ContinuizeEpochs(self.svel);

[time, ~] = size(speed);

if Fs ~= 1
    begin = 1;
    ending = Fs;
    count = 1;
    while ending <= time
        newSpeed(count,1) = nanmean(speed(begin:ending));
        for i = 1:Ncells
            newFR(count,i) = nanmean(FR(begin:ending,i));
        end
        begin = ending + 1;
        ending = begin + Fs - 1;
        count = count + 1;
    end
    speed = newSpeed;
    FR = newFR;
end

finalFR = [FR, ones(length(speed),1)];

f = finalFR \ speed;

predictedSpeed = finalFR * f;

figure()
plot(speed, 'b')
hold on
plot(predictedSpeed, 'r')
legend('Actual Speed', 'Predicted Speed')

end