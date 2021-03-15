function [theta,S_dens,Ch,PI,int,s_hist,s_hist_unc, PI_min] = Watson_pragmatic( varargin )

%%% s-hist - histogram smoothed raw data

%%% s_hist_unc - histogram with s phase subtracted

%%% S_dens - probability density function of the s phase

%%% int - linear spacing of the frequencies

%%% Ch - channel index

%%% PI - original input data

%% setup data
PI = varargin{1};

if length(varargin)>1
    opt = varargin{2};
else
    opt =[];
end

n_Ch = round(length(PI)./20);
Ch = 1:n_Ch;

PI_max = max(PI);
PI_min = max([0,min(PI)]);

int = linspace(PI_min,PI_max,n_Ch);

histogram = hist(PI,int);
%% Step 1

for i=1:length(histogram)
    temp(:,i) = histogram(i).*normpdf(Ch,Ch(i),1.5);
end

s_hist = sum(temp,2);

clear temp

%% Step 2

Ch_min = 1;
Ch_max = n_Ch;

%% Step 3

if isempty(opt)
    mu_G1_0 = find(s_hist == max(s_hist));
    sd_G1_0 = find(s_hist > 0.6*max(s_hist),1,'first'); 
    sd_G1_0_width = mu_G1_0 - sd_G1_0;
else
    mu_G1_0 = opt;
    sd_G1_0 = find(s_hist > 0.6*s_hist(opt),1,'first');
    sd_G1_0_width = mu_G1_0 - sd_G1_0;
end

%% Step 4

lb = max([Ch_min,find(Ch > mu_G1_0 - 3.*sd_G1_0_width,1,'first')]); % find 3 std below max
ub = find(Ch > mu_G1_0 + sd_G1_0_width,1,'first'); % find 1 std above max

Ch_reg     = Ch(lb:ub); % channel fit range
s_hist_reg = s_hist(lb:ub);  % histogram data fit range

theta_0 = log([mu_G1_0, sd_G1_0_width]); % log of estimated mean and width 
lb_theta = log([lb, 1e-3]); 
ub_theta = log([ub, (ub-lb)./2]);

% disable output message
opts =  optimset('display','off');

[theta] = lsqnonlin(@(xi)fit_normal(exp(xi),s_hist_reg,Ch_reg),theta_0,lb_theta,ub_theta,opts); % theta fit parameters mean and std deviation
% log of mean and std dev from fit
mu_G1 = round(exp(theta(1))); % actual mean G1 phase
sd_G1_width = round(exp(theta(2))); % actual std G1 phase

clear lb ub s_hist_reg Ch_reg

%% Step 5

lb = find(Ch > 1.5.*mu_G1,1,'first'); % find next max 1.5 std dev from G1 peak

s_hist_reg = s_hist(lb:end);

mu_G2_0 = find(s_hist == max(s_hist_reg)); % approximate mean G2
sd_G2_0 = find(s_hist > 0.6*max(s_hist_reg),1,'last');
sd_G2_0_width = sd_G2_0 - mu_G2_0; % approximate std G2

clear lb s_hist_reg

%% Step 6

lb = find(Ch > mu_G2_0 - sd_G2_0_width,1,'first');
ub = min([Ch_max,find(Ch > mu_G2_0 + 3.*sd_G2_0_width,1,'first')]);

Ch_reg     = Ch(lb:ub);
s_hist_reg = s_hist(lb:ub);

theta_0  = log([mu_G2_0, sd_G2_0_width]);
lb_theta = log([lb, 1e-3]);
ub_theta = log([ub, (ub-lb)./2]);

[theta] = lsqnonlin(@(xi)fit_normal(exp(xi),s_hist_reg,Ch_reg),theta_0,lb_theta,ub_theta,opts);

mu_G2 = round(exp(theta(1)));
sd_G2_width = round(exp(theta(2)));

clear lb ub s_hist_reg Ch_reg

%% Step 7
check = 0;


while check == 0 % loop to check means and std deviations and update if necesscary
    
    
    kG1 = -3; % starting value of kG1
    kG2 = -3; % starting value of kG2

    G1_dens = normpdf(Ch,mu_G1,sd_G1_width); % The fitted Gaussian for G1 phase
    G2_dens = normpdf(Ch,mu_G2,sd_G2_width);
%     crit = 0;
    

    %kG1
    cnt = 0;
    
    figure(300); plot(Ch, G1_dens/max(G1_dens)*s_hist(mu_G1)); hold on
    plot(Ch, s_hist)
    plot(Ch, G2_dens/max(G2_dens)*s_hist(mu_G2))
    
    sfitPar = polyfit(Ch(mu_G1 + 3*sd_G1_width: mu_G2 - 3*sd_G2_width), s_hist(mu_G1 + 3*sd_G1_width: mu_G2 - 3*sd_G2_width), 1);
    plot(Ch, Ch*sfitPar(1) + sfitPar(2))
    sApproxG1 = mu_G1*sfitPar(1) + sfitPar(2);
    sApproxG2 = mu_G2*sfitPar(1) + sfitPar(2);
    
    crit = sApproxG1/((2*G1_dens(mu_G1)/max(G1_dens))*s_hist(mu_G1));
    
    while 0.5*(erf(-kG1) + 1) >= crit % updates the 

        S_dens = 0.5.*(erf(((Ch-mu_G1)./sd_G1_width)-kG1)...
            -erf(((Ch-mu_G2)./sd_G2_width)+kG2));
        
%         kG1 
%         S_dens(mu_G1)
%         (2.*G1_dens(mu_G1))
%         erf(kG1)
        kG1 = kG1+0.05;

        %crit = S_dens(mu_G1)./(2.*G1_dens(mu_G1));
%         crit = sApproxG1/((2*G1_dens(mu_G1)/max(G1_dens))*s_hist(mu_G1));
        
%         figure(10001), plot(cnt,erf(kG1), 'ko', cnt,S_dens(mu_G1)./(2.*G1_dens(mu_G1)), 'b+', ...
%             cnt,S_dens(mu_G1), 'rv' ), hold on
        
%         erf(kG1)
%         S_dens(mu_G1)
%         2.*G1_dens(mu_G1)
cnt = cnt + 1;
    end
%     legend('erf(kG1)', 'S/2G1', 'S')
    %kG2

    crit = sApproxG2/((2*G2_dens(mu_G2)/max(G2_dens))*s_hist(mu_G2));

 cnt = 0;
    while 0.5*(-erf(kG2) + 1) >= crit % updates the 
        
        S_dens = 0.5.*(erf(((Ch-mu_G1)./sd_G1_width)-kG1)...
            -erf(((Ch-mu_G2)./sd_G2_width)+kG2)); % don't need to update every time step
        
        kG2 = kG2 + 0.05;
        cnt = cnt + 1;
    end
    
%     while erf(kG2) <= crit
% 
%         S_dens = 0.5.*(erf(((Ch-mu_G1)./sd_G1_width)-kG1)...
%             -erf(((Ch-mu_G2)./sd_G2_width)+kG2));
% 
%         kG2 = kG2+0.05;
% 
%         crit = S_dens(mu_G1)./(2.*G1_dens(mu_G1));
%         %crit = S_dens(mu_G2)./(2.*G2_dens(mu_G2));
%     end

    %close all
    %figure(1)
    %hold on
    %plot(Ch,S_dens)
    %plot(Ch,G1_dens,'g')
    %plot(Ch,G2_dens,'r')
    %hold off


%% Step 8

    s_hist_unc = s_hist - s_hist.*S_dens';

    %figure(2)
    %plot(Ch,s_hist)
    %hold on
    %plot(Ch,s_hist_unc,'c')

 
%% Step 9

    theta_0  = log([mu_G1, sd_G1_width, mu_G2, min([sd_G1_width,sd_G2_width])]);
    lb_theta = log(1e-3.*ones(4,1));
    ub_theta = log([n_Ch;n_Ch;n_Ch;sd_G1_width]);

    [theta] = lsqnonlin(@(xi)fit_sum_of_normals(exp(xi),s_hist_unc,Ch),theta_0,lb_theta,ub_theta,opts);

    Mu_G1 = round(exp(theta(1)));
    SD_G1_width = round(exp(theta(2)));
    Mu_G2 = round(exp(theta(3)));
    SD_G2_width = round(exp(theta(4)));

    %figure(3)
    %plot(Ch,s_hist_unc)
    %hold on
    %plot(Ch,s_hist_unc(round(Mu_G1))./normpdf(Mu_G1,Mu_G1,SD_G1_width).*normpdf(Ch,Mu_G1,SD_G1_width),'g')
    %plot(Ch,s_hist_unc(round(Mu_G2))./normpdf(Mu_G2,Mu_G2,SD_G2_width).*normpdf(Ch,Mu_G2,SD_G2_width),'r')
    %hold off
    
    if abs((Mu_G1-mu_G1)./n_Ch)<0.025 && abs((SD_G1_width-sd_G1_width)./n_Ch)<0.025 && abs((Mu_G2-mu_G2)./n_Ch)<0.025 && abs((SD_G2_width-sd_G2_width)./n_Ch)<0.025
        check = 1;
    else
        mu_G1 = Mu_G1;
        sd_G1_width = SD_G1_width;
        mu_G2 = Mu_G2;
        sd_G2_width = SD_G2_width;
    end
    
    theta(1) = Mu_G1;
    theta(2) = SD_G1_width;
    theta(3) = Mu_G2;
    theta(4) = SD_G2_width;
    
end


end

function [sos] = fit_normal(theta,hist,Ch)

n_pdf  = normpdf(theta(1),theta(1),theta(2));
n_hist = max(hist); 

sos = (normpdf(Ch,theta(1),theta(2))'./n_pdf - hist./n_hist);

end

function [sos] = fit_sum_of_normals(theta,hist,Ch)

mu1 = theta(1);
sd1 = theta(2);
mu2 = theta(3);
sd2 = theta(4);

pdf1  = normpdf(Ch,mu1,sd1);
pdf2  = normpdf(Ch,mu2,sd2);

a = hist(round(mu1))./normpdf(mu1,mu1,sd1);
b = hist(round(mu2))./normpdf(mu2,mu2,sd2); 

sos = a.*pdf1+b.*pdf2-hist';

end
