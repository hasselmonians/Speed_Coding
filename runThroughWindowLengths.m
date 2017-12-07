tic
winder = nan(1,10);
winder(1) = .5;
for i = 2:10
    winder(i) = winder(i-1) * 2;
end

shift = nan(1,length(winder));
for i = 1:length(winder)
    shift(i) = .1 * winder(i);
end

correls_tenP = nan(1,length(winder));
for i = 1:length(winder)
    [~, correl] = checkSpeedModulation(root, winder(i), shift(i));
    correls_tenP(i) = nanmean(correl);
end

shift = repmat(.2,1,length(winder));
correls_const = nan(1,length(winder));
for i = 1:length(winder)
    [~, correl] = checkSpeedModulation(root, winder(i), shift(i));
    correls_const(i) = nanmean(correl);
end

shift = repmat(.5,1,length(winder));
correls_const5 = nan(1,length(winder));
for i = 1:length(winder)
    [~, correl] = checkSpeedModulation(root, winder(i), shift(i));
    correls_const5(i) = nanmean(correl);
end

figure()
subplot(1,3,1)
plot(winder, correls_tenP, 'k*-')
% xlabel('Window Width (s)')
ylabel('Average Correlation Coefficient')
title('Shift = 10% Window Width')

% figure()
subplot(1,3,2)
plot(winder, correls_const, 'k*-')
xlabel('Window Width (s)')
% ylabel('Average Correlation Coefficient')
title('Shift = 200ms')

% figure()
subplot(1,3,3)
plot(winder, correls_const5, 'k*-')
% xlabel('Window Width (s)')
% ylabel('Average Correlation Coefficient')
title('Shift = 500ms')
toc