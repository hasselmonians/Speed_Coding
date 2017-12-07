function [R, root] = InstFR(root, ifPlot)

    if ~exist('ifPlot','var')
        ifPlot = 0;
    end

    %% Instantaneous firing rate
    root.epoch = [-inf inf];
    fr = zeros(size(root.x));
    temp = root.spike(root.cel(1), root.cel(2)).i;
    tt = unique(temp);
    for i = 1:length(tt)
        fr(tt(i)) = sum(temp == tt(i));
    end

    sigma = 33/(1000/(root.fs_video));
    ssize = floor(10*sigma);
    x = linspace(-ssize / 2, ssize / 2, ssize);
    gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
    gaussFilter = gaussFilter / sum (gaussFilter); % normalize
    frFilt = conv (fr, gaussFilter, 'same');
    frFilt = frFilt * root.fs_video;
    root.b_myvar = frFilt;
    
    %% 
    R = corr(root.vel, root.myvar);
    
    if ifPlot == 1
        sp = [0, prctile(root.vel, 95)];
        [y_m, R, P, b, y_int] = CMBHOME.Utils.LinearRegression(root.vel, root.myvar);
        ym = sp * b + y_int;
        plot(sp, ym)
        hold on

    end
    
end
