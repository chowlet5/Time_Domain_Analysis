function[max_est,min_est] = peak(data,dur_ratio)

if nargin==0 || isempty(data)
    error('Time series input expected');
end

plot_on = 0; % turns plotting on (1) [for diagnostics] or off (0)

num_CDF = 1000;                                % number of points in CDF for peaks
min_CDF = 0.0005;                              % minimum value of CDF for peaks
max_CDF = 0.9995;                              % maximum value of CDF for peaks
CDF_pk = linspace(min_CDF,max_CDF, num_CDF);   % linearly spaced CDF values for peaks

sdata = size(data);
max_est = zeros(sdata(1),1);
min_est = zeros(sdata(1),1);

for i=1:sdata(1)
    X = data(i,:);
    n = length(X);
    avg_X = mean(X);
    sorted_X = sort(X);                                   % sort in ascending order:
    
    CDF_X = (1:n)/(n+1);                                  % Empirical Cumulative Distribution Function
    snv = isnc(CDF_X);                                    % standard normal variate which is [(x-mu)/sigma] is obtained 
    avg_snv = mean(snv);
    % linear regression:
    sigma = (sum(snv(:).*sorted_X(:))-n*avg_snv*avg_X)/(sum(snv.^2)-n*avg_snv^2);
    mu = avg_X - sigma*avg_snv;
    X_fit = mu + sigma*snv;
              
    % Probability Plot Correlation Coefficient:
    norm_PPCC = sigma*std(snv)/std(sorted_X);

    if plot_on
        figure(3)
        plot(snv,sorted_X,'.',snv,X_fit,'-');
        xlabel('Standard normal variate');
        ylabel('Value of time series');
        title({'Normal distribution fit',...
            ['Probability Plot Correlation Coefficient: ' num2str(norm_PPCC)]});
        legend({'Data','Best-fit Normal distribution'},'Location','NorthWest')
        set(gcf,'Name','Fit of normal distribution(press any key to continue...)',...
            'NumberTitle','off');
        pause;
        figure(4);
        plot(sorted_X,CDF_X,'k.',X_fit,CDF_X,'r-');
        xlabel('Value of time series');
        ylabel('Cumulative Probability');
        legend({'Empirical CDF: all data','Normal distribution fit'},...
            'Location','SouthOutside');
        set(gcf,'Name','Cumulative Distribution Function: empirical and fitted (press any key to continue...)',...
            'NumberTitle','off');
        pause;
    end
%--------------------------------------------------------------------------
%Estimate the mean zero upcorssing rate of a process y(t) with standard
%normal probability distribution using the classical Rice(1954) results as
%follow:
%--------------------------------------------------------------------------
%Estimate the interval of integration in frequency domain. This is done
%by matching variances obtained from the timehistory data and integration
%of the spectra in the frequency domain
%Variance from time history
stdX = std(X);   % actual std from the time history data
varX = stdX.*stdX;
%Variance from frequency domain 
df = 65536;
fs=5.12;
[S_X,f]=pwelch(X,[],[],df,fs);
sf=(df/2)+1; 
si=1;
var_X=(fs/df)*sum(S_X(si:sf));         % total q from area of Sq
temp=var_X;
si=2;
var_X=(fs/df)*sum(S_X(si:sf));
while (abs(var_X-varX)<abs(temp-varX)&& var_X>varX)
    temp=var_X;
    si=si+1;        
    var_X=(fs/df)*sum(S_X(si:sf)); 
end
var_X=temp;
%Integration limits in frequeny domain
si=si-1;
sf=sf;
%--------------------------------------------------------------------------
%Mean upcrossing rate for process y(t)
numer=trapz(f(si:sf),f(si:sf).*f(si:sf).*S_X(si:sf));
denom=trapz(f(si:sf),S_X(si:sf));
nu_y=sqrt(numer/denom);   
%maximum peaks of process y(t) corresponding to specified cumulative probabilities
y_pk = sqrt(2.0*log(-dur_ratio*nu_y*3600./ log(CDF_pk))); 
%Mapping peak values from a non-Gaussian process X(t) to a Gaussian process y(t)
X_max = y_pk*sigma + mu;
X_min = -y_pk*sigma + mu;
%Compute mean of the Peaks for process X(t)
pdf_pk = -y_pk .* CDF_pk .* log(CDF_pk);
max_est(i) = trapz(y_pk,pdf_pk.*X_max);
min_est(i) = trapz(y_pk,pdf_pk.*X_min);
    
 end