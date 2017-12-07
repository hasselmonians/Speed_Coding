function [ stats ] = fit_speed( spk_ts,Fs,ts,speed,spk_speed,plotit )
%FIT_SPEED Fits the firing rate versus speed of a neuron
%
% INPUT
%   spk_ts - The spike time stamps
%   Fs - The sampling frequency
%   ts - The sampling timestamps
%   speed - The speed of the animal at each time (cm/second)
% 
% OPTIONAL
%   plotit (false) - Whether to plot
%
% RETURNS
%    stats: Output structure with the following fields:
%       F_binned_linear - F statistic for the linear binned fit
%       F_binned_satexp - F statistic for the saturating exponential binned
%           fit
%       F_binned_satexp_vs_binned_lin - F statistic comparing the
%           saturating exponential binned fit to the linear fit
%       F_lin - F statistic for the linear fit
%       F_satexp - F statistic for the saturating exponential fit
%       F_satexp_v_linear - F statistic for the saturating exponential fit
%           over the linear fit
%       Fr - Mean firing rate
%       LL0 - Log-likelihood of uniform rate
%       LL_lin - Log-likelihood of linear fit
%       LL_satexp - Log-likelihood of saturating exponential fit
%       R2_binned_linear - Coefficient of determination for binned linear
%           fit
%       R2_binned_satexp - Coefficient of determination for binned
%           saturating exponential fit
%       R2_lin - Nagelkerke/Cragg & Uhler's Pseudo R^2 for the linear fit
%       R2_satexp - Nagelkerke/Cragg & Uhler's Pseudo R^2 for the 
%           saturating exponential fit
%       R2_satexp_v_linear - Nagelkerke/Cragg & Uhler's Pseudo R^2 for the 
%           saturating exponential fit versus the linear fit
%       binned_lin_intercept - Linear intercept fit from binned data
%       binned_lin_slope - Linear slope fit from binned data
%       binned_satexp_a - Saturating exponential parameter a from binned
%           data
%       binned_satexp_b - Saturating exponential parameter b from binned
%           data
%       binned_satexp_d - Saturating exponential parameter d from binned
%           data
%       ci_linear - 95% confidence intervals for the linear fit parameters,
%           from the Fisher Information
%       ci_satexp - 95% confidence intervals for the saturating exponential
%           fit parameters, from the Fisher Information
%       hzspeed - The firing rate versus speed binned
%       mle_linear - The parameters for the linear fit ([intercept slope])
%       mle_satexp - The parameters for the saturating exponential fit
%           ([d a b])
%       occupancy - The amount of time spent in each speed bin
%       pF_lin - significance of the linear fit
%       pF_satexp - significance of the saturating exponential fit
%       pF_satexp_v_linear - significance of the saturating exponential fit
%           over the linear fit
%       p_binned_linear - significance of the binned linear fit
%       p_binned_satexp - significance of the binned saturating exponential
%           fit
%       p_binned_satexp_vs_binned_lin - significance of the binned
%           saturating exponential fit over the binned linear fit
%       p_rho_linear - Significance of correlation between binned rates and
%           binned speed
%       rho_linear - Correlation coefficient between the binned rate and
%           the binned speed
%       speed_dim - The speed bins
%   
% Copyright (c) 2015 Trustees of Boston University
% 
% This file is part of mle_rhythmicity.
% 
% This version of mle_rhythmicity is solely for the purposes of review and 
% demonstration of its functionality by the editors and reviewers assigned 
% by Neuron (Elsevier Inc.). Redistribution to others and other uses in 
% source and binary forms, with or without modification, is prohibited. 
% Upon publication, the associated paper will include a different link 
% from which the code will be permanently available, distributed under 
% the open-source BSD license (Available at 
% http://opensource.org/licenses/bsd-license.php) allowing future users to 
% freely distribute and modify the code.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS 
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANOnly frameTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.   
%% Parse input
warning off curvefit:fit:noStartPoint;

if ~exist('plotit','var')
   plotit = false; 
end

% Speed bins
speed_dim = (2:1:quantile(speed,0.95))';

% Velocity at each spike
%spk_speed = interp1(ts(~isnan(speed)),speed(~isnan(speed)),spk_ts,'nearest');

% Number of spikes in each speed bin
count = histc(spk_speed, speed_dim);

% Number of frames spent in each speed bin
occupancy = histc(speed, speed_dim);

% Firing rate in each speed bin
hzspeed = count./(occupancy/Fs);

% Average firing rate
Fr = numel(spk_ts)/(sum(abs(diff(ts)-1/Fs)<1/Fs)/Fs);% Calulates duration of session by time frames with small differences

% Count of spikes in each frame
cnt = histc(spk_ts,ts);

% Linear fit of binned data
temp = polyfit(speed_dim(1:end-1),hzspeed(1:end-1),1);

% Calculate statistics from the linear fit
binned_lin_intercept = temp(2);% Y-intercept
binned_lin_slope = temp(1);% Slope (Hz/(cm/sec))

% Correlations
[rho_linear,p_rho_linear] = corrcoef([speed_dim(1:end-1),hzspeed(1:end-1)]);
% [rho_linear,p_linear] = corrcoef([polyval(temp,speed_dim(1:end-1),2) hzspeed(1:end-1)]);
rho_linear = rho_linear(1,2);
p_rho_linear = p_rho_linear(1,2);

R2_binned_linear =1-sum((polyval(temp,speed_dim(1:end-1),2)-hzspeed(1:end-1)).^2)/sum((hzspeed(1:end-1)-mean(hzspeed(1:end-1))).^2);% Coefficient of determination

% Testing - use the F-test with estimated dispersion
F_binned_linear = ((sum((hzspeed(1:end-1)-mean(hzspeed(1:end-1))).^2)-sum((polyval(temp,speed_dim(1:end-1),2)-hzspeed(1:end-1)).^2))/1)/(sum((polyval(temp,speed_dim(1:end-1),2)-hzspeed(1:end-1)).^2)/(numel(hzspeed)-1-2));
p_binned_linear = 1-fcdf(F_binned_linear,1,numel(hzspeed)-1-2);

% Inhomogenous poisson of unbinned data linear fit 
[mle_linear,ci_linear] = mle(cnt,'pdf',@(~,b0,b1)max(poisspdf(cnt,max([ones(size(speed)) speed]*[b0;b1],realmin)/Fs),realmin),'start',[binned_lin_intercept binned_lin_slope]);%,'options',statset('Display','iter'));

% Calulcate statistics from inhomogenous Poisson fit
LL0 = sum(log(poisspdf(cnt,Fr/Fs)));% Log likelihood of "dumb", uniform fit
LL_lin = sum(log(poisspdf(cnt,[ones(size(speed)) speed]*mle_linear(:)/Fs)));% LL of linear fit
F_lin = 2*(LL_lin-LL0);% F statistic for the linear fit
pF_lin = 1-fcdf(F_lin,1,numel(cnt)-2-1);% Significance of F statistic of the linear fit
% R2_lin = (1-exp(2/numel(cnt)*(LL0-LL_lin)))/(1-exp(2/numel(cnt)*LL0));%  Nagelkerke/Cragg & Uhler's Pseudo R^2
R2_lin=(var(cnt)-mean((cnt-max(([ones(size(speed)) speed]*mle_linear(:))/Fs,realmin)).^2))/(var(cnt)-mean(cnt));

% Saturating exponential fit of binned data
 f = fittype('d - a*exp(-b*x)','independent','x','coefficients',{'d' 'a' 'b'});
 options = fitoptions(f);
 options.lower = [0 0 0];
 options.upper = [1.5*nanmax(hzspeed) inf inf];
 [fitobject, gof] = fit(speed_dim(~isnan(hzspeed)),hzspeed(~isnan(hzspeed)),f,options);
 % Do it 100 times and take the best one
 for i=1:100
%      iFr 
     [fitobject2, gof2] = fit(speed_dim(~isnan(hzspeed)),hzspeed(~isnan(hzspeed)),f,options);
     if gof2.rsquare>gof.rsquare
        fitobject = fitobject2;
        gof = gof2;
     end
 end
 
 % Calculate statistics from the saturating exponential fit of binned data
 binned_satexp_d = fitobject.d;
 binned_satexp_a = fitobject.a;
 binned_satexp_b = fitobject.b;
 
R2_binned_satexp = 1-sum((fitobject(speed_dim(~isnan(hzspeed)))-hzspeed(~isnan(hzspeed))).^2)/sum((hzspeed(1:end-1)-mean(hzspeed(1:end-1))).^2);% Coefficient of determination
F_binned_satexp = ((sum((hzspeed(1:end-1)-mean(hzspeed(1:end-1))).^2)-sum((fitobject(speed_dim(~isnan(hzspeed)))-hzspeed(1:end-1)).^2))/2)/(sum((fitobject(speed_dim(~isnan(hzspeed)))-hzspeed(1:end-1)).^2)/(numel(hzspeed)-1-3));% F statistic
p_binned_satexp = 1-fcdf(F_binned_satexp,2,numel(hzspeed)-1-3);% Significance of f-statistic
F_binned_satexp_vs_binned_lin = ...
    max(((sum((polyval(temp,speed_dim(1:end-1),2)-hzspeed(1:end-1)).^2)-...
    sum((fitobject(speed_dim(~isnan(hzspeed)))-hzspeed(1:end-1)).^2))/1)/...
    (sum((fitobject(speed_dim(~isnan(hzspeed)))-hzspeed(1:end-1)).^2)/(numel(hzspeed)-1-3)),0);% F-statistic comparing linear verus saturating exponential fits
p_binned_satexp_vs_binned_lin = 1-fcdf(F_binned_satexp_vs_binned_lin, 1, numel(cnt)-3);% Significance of comparison

% Inhomogenous poisson of unbinned data saturating exponential fit 
[mle_satexp, ci_satexp] = mle(cnt,'pdf',@(~,d,a,b)max(poisspdf(cnt,max((d-a*exp(-b*speed))/Fs,realmin)),realmin),'start',[binned_satexp_d binned_satexp_a binned_satexp_b],'lowerbound',[0 0 0],'upperbound',[1.5*nanmax(hzspeed) inf inf]);%,'options',statset('Display','iter'));

% Calculate statistics for the saturating exponential fit
LL_satexp = sum(log(max(poisspdf(cnt,max((mle_satexp(1)-mle_satexp(2)*exp(-mle_satexp(3)*speed))/Fs,realmin)),realmin)));% Log-likelihood for the saturating exponential fit
F_satexp = LL_satexp-LL0;% F-statistic compared to dumb, uniform fit
pF_satexp = 1-fcdf(F_satexp, 2, numel(cnt)-3-1);% Significance of F compared to uniform fit
F_satexp_v_linear = max(2*(LL_satexp-LL_lin),0);% F-statistic comparing to the linear fit
pF_satexp_v_linear = 1-fcdf(F_satexp_v_linear, 1, numel(cnt)-3-1);% Significance of comparison to the linear fit
% R2_satexp = (1-exp(2/numel(cnt)*(LL0-LL_satexp)))/(1-exp(2/numel(cnt)*LL0));%  Nagelkerke/Cragg & Uhler's Pseudo R^2
% R2_satexp_v_linear = (1-exp(2/numel(cnt)*(LL_lin-LL_satexp)))/(1-exp(2/numel(cnt)*LL_lin));% Nagelkerke/Cragg & Uhler's Pseudo R^2 comparing to the linear fit
R2_satexp=(var(cnt)-mean((cnt-max((mle_satexp(1)-mle_satexp(2)*exp(-mle_satexp(3)*speed))/Fs,realmin)).^2))/(var(cnt)-mean(cnt));


%% Pack output

% clear ans;
% clear Fs;
% clear cnt;
% clear count;
% clear f;
% clear fitobject;
% clear fitobject2;
% clear gof;
% clear gof2;
% clear i;
% clear options;
% %clear plotit;
% clear speed;
% clear spk_speed;
% clear spk_ts;
% clear temp;
% clear ts;
% clear c;

c = who;
stats = struct;
for i = c'
    eval(['stats.' i{1} '=' i{1} ';']);
end

warning on curvefit:fit:noStartPoint;

%% plotting stuff
if plotit==1
    speed_dim_centers = (speed_dim(1:end-1)+speed_dim(2:end))/2;

    lowerbound_lin = NaN(size(speed_dim_centers));
    upperbound_lin = NaN(size(speed_dim_centers));
    lowerbound_satexp= NaN(size(speed_dim_centers));
    upperbound_satexp= NaN(size(speed_dim_centers));

    for i=1:numel(speed_dim_centers)
        [~,lowerbound_lin(i)]=fmincon(@(phat)[1 speed_dim_centers(i)]*phat(:),mle_linear,[],[],[],[],[],[],@(phat)nonlcon(speed,phat,LL_lin,@(phat)...
            sum(log(poisspdf(cnt,max([ones(size(speed)) speed]*phat(:),realmin)/Fs))),numel(cnt)),optimoptions(@fmincon,'Algorithm','sqp','display','off'));
        [~,upperbound_lin(i)]=fmincon(@(phat)-([1 speed_dim_centers(i)]*phat(:)),mle_linear,[],[],[],[],[],[],@(phat)nonlcon(speed,phat,LL_lin,@(phat)...
            sum(log(poisspdf(cnt,max([ones(size(speed)) speed]*phat(:),realmin)/Fs))),numel(cnt)),optimoptions(@fmincon,'Algorithm','sqp','display','off'));
        [~,lowerbound_satexp(i)]=fmincon(@(phat)phat(1)-phat(2)*exp(-phat(3)*speed_dim_centers(i)),mle_satexp,[],[],[],[],[0 0 0],[1.5*nanmax(hzspeed) inf inf],@(phat)nonlcon(speed,phat,LL_lin,@(phat)...
            sum(log(poisspdf(cnt,max((phat(1)-phat(2)*exp(-phat(3)*speed))/Fs,realmin)))),numel(cnt)),optimoptions(@fmincon,'Algorithm','sqp','display','off'));
        [~,upperbound_satexp(i)]=fmincon(@(phat)-(phat(1)-phat(2)*exp(-phat(3)*speed_dim_centers(i))),mle_satexp,[],[],[],[],[0 0 0],[1.5*nanmax(hzspeed) inf inf],@(phat)nonlcon(speed,phat,LL_lin,@(phat)...
            sum(log(poisspdf(cnt,max((phat(1)-phat(2)*exp(-phat(3)*speed))/Fs,realmin)))),numel(cnt)),optimoptions(@fmincon,'Algorithm','sqp','display','off'));
    end
    upperbound_lin = -upperbound_lin;
    upperbound_satexp = -upperbound_satexp;

    lowerbound_lin_rates = poissinv(0.025,lowerbound_lin.*occupancy(1:end-1)/Fs)./(occupancy(1:end-1)/Fs);
    upperbound_lin_rates = poissinv(0.975,upperbound_lin.*occupancy(1:end-1)/Fs)./(occupancy(1:end-1)/Fs);
    lowerbound_satexp_rates = poissinv(0.025,lowerbound_satexp.*occupancy(1:end-1)/Fs)./(occupancy(1:end-1)/Fs);
    upperbound_satexp_rates = poissinv(0.975,upperbound_satexp.*occupancy(1:end-1)/Fs)./(occupancy(1:end-1)/Fs);

        
    
    % linear
    figure
    hold on
    p_lbnd = poissinv(0.025,hzspeed(1:end-1).*occupancy(1:end-1)/Fs)./(occupancy(1:end-1)/Fs);
    p_upnd = poissinv(0.975,hzspeed(1:end-1).*occupancy(1:end-1)/Fs)./(occupancy(1:end-1)/Fs);
    plot(speed_dim_centers, hzspeed(1:end-1),'ko','MarkerFaceColor','k','MarkerSize',8)
    x = speed_dim_centers(:); x2 = [x; flipud(x)]; 
    y2 = [p_lbnd(:); flipud(p_upnd(:))];
    patch(x2,y2,'k','FaceAlpha',0.4,'EdgeColor','k','EdgeAlpha',0.15,'LineWidth',1);
    x = speed_dim_centers(:); x2 = [x; flipud(x)]; 
    y2 = [lowerbound_lin(:); flipud(upperbound_lin(:))];
    patch(x2,y2,'b','FaceAlpha',0.4,'EdgeColor','b','EdgeAlpha',0.3,'LineWidth',1);
    plot(speed_dim_centers,speed_dim_centers*stats.mle_linear(2) + stats.mle_linear(1),'b','LineWidth',3)
    set(gca,'LineWidth',3)
    
    % Saturating
    figure
    hold on
    p_lbnd = poissinv(0.025,hzspeed(1:end-1).*occupancy(1:end-1)/Fs)./(occupancy(1:end-1)/Fs);
    p_upnd = poissinv(0.975,hzspeed(1:end-1).*occupancy(1:end-1)/Fs)./(occupancy(1:end-1)/Fs);
    plot(speed_dim_centers, hzspeed(1:end-1),'ko','MarkerFaceColor','k','MarkerSize',8)
    x = speed_dim_centers(:); x2 = [x; flipud(x)]; 
    y2 = [p_lbnd(:); flipud(p_upnd(:))];
    patch(x2,y2,'k','FaceAlpha',0.4,'EdgeColor','k','EdgeAlpha',0.15,'LineWidth',1);
    x = speed_dim_centers(:); x2 = [x; flipud(x)]; 
    y2 = [lowerbound_satexp(:); flipud(upperbound_satexp(:))];
    patch(x2,y2,'r','FaceAlpha',0.4,'EdgeColor','r','EdgeAlpha',0.3,'LineWidth',1);
    plot(speed_dim_centers, stats.mle_satexp(1) - stats.mle_satexp(2)*exp(-stats.mle_satexp(3)*speed_dim_centers),'r','LineWidth',3)
    set(gca,'LineWidth',3)
    
    1+1
    
end
%% re-Pack output, because effort

% clear ans;
% clear Fs;
% clear cnt;
% clear count;
% clear f;
% clear fitobject;
% clear fitobject2;
% clear gof;
% clear gof2;
% clear i;
% clear options;
% %clear plotit;
% clear speed;
% clear spk_speed;
% clear spk_ts;
% clear temp;
% clear ts;
% clear c;

c = who;
stats = struct;
for i = c'
    eval(['stats.' i{1} '=' i{1} ';']);
end

warning on curvefit:fit:noStartPoint;


end

function [c,ceq] = nonlcon(speed,phat,LLA,LLF,n)

    ceq = 0;
    c = 2*(LLA-LLF(phat))/numel(phat)-finv(1-0.05,numel(phat),n-numel(phat)-1);

end
