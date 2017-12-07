% works on root object and speedMod data set

close all

% Parameters
epoch = Non.S1; % epochs
time_window = 10; % time window in s
cel = [1 1]; % specify cell in root object
spkTrain = Speed_coding_project.get_spikeTrain_withPoissonNoise(root,cel,speedMod); % give artificial spike train

% get short time scale speed tuning properties for all time windows within
% epochs
[F,V,~,p,SSresid,R,rsq,b,y_int,F_total,V_total,~] = Speed_coding_project.get_art_spkTrain_shortTimeScale_speedTuning(root,epoch,time_window,spkTrain); % concatenate results for each time window
log_i = zeros(length(F),1); % create logical array to indicate wich cells contain nan arrays
for i = 1:length(F)
    if ~isnan(F{i})
        log_i(i) = 1;
    end
end
log_i = logical(log_i);
F = F(log_i);
V = V(log_i);
% N = N(log_i);
R = R(log_i);
F_matrix = cell2mat(F');
V_matrix = cell2mat(V);
% N_matrix = cell2mat(N');
R_matrix = cell2mat(R');

figure
hist(R_matrix,40)
xlabel('R')

figure
hold on
Fn = cell(1,size(V,2));
for i = 1:length(V)
    %     plot(V{i},F{i},'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 25, 'MarkerEdgeColor', 'k'), xlabel('Running Speed (cm/s)'), ylabel('Firing Frequency'), hold on
    if ~isnan(V{i})
        [~,~,~,Fn{i}] = fit_regression_line(V{i},F{i}');
        line(V{i}, Fn{i})
    end
end
% calculate average across regression lines
V_i = cell2mat(V);
F_i = cell2mat(F');
[A,xbar] = get_histTwoVectorInput(V_i,F_i,[2:1:50]);
plot(xbar,A,'Color','r','LineWidth',3)

figure
title('Calculated from concatenated windows')
[R_total,b_total,p_total,Fn_total] = fit_regression_line(V_total,F_total');
plot(V_total, F_total, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 25, 'MarkerEdgeColor', 'k'), xlabel('Running Speed (cm/s)'), ylabel('Firing Frequency')
% title(strcat(parent_folder,',',name,', T',num2str(cel(1,1)),'C',num2str(cel(1,2))),'Interpreter','none');
hold on
line(V_total, Fn_total, 'Color', 'r', 'LineWidth', 3), hold on
ys = ylim;
xs = xlim;
text(xs(2) - .9*diff(xs), ys(2)-.1*diff(ys), ['R = ' num2str(R_total, 3)], 'FontSize', 12)
text(xs(2) - .9*diff(xs), ys(2)-.2*diff(ys), ['slope = ' num2str(b_total, 3)], 'FontSize', 12)
text(xs(2) - .9*diff(xs), ys(2)-.3*diff(ys), ['p = ' num2str(p_total, 3)], 'FontSize', 12)
axis square


