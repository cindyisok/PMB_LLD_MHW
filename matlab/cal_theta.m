function [mclim,m90] = cal_theta(temp,cli_start,cli_end,varargin)

% Description
%
% [mclim,m90] = cal_theta(temp,cli_start,cli_end) returns the climatology
% (lon x lat x 366) and percentile (lon x lat x 366) of the fixed 31-year 
% baseline. 

% Input Arguments
%
% temp - the detrended SST (lon x lat x 31year) from the fixed 31-yr baseline. 
%
% cli_start - datenum(start_yr,1,1) the first day of climatology period.
% cli_end - datenum(end_yr,12,31) the last day of climatology period.
%
% Threshold - Default is 0.9. Threshold percentile to detect MHW events.
%
% windowHalfWidth - Default is 5. Width of sliding window to calculate
% climatology and threshold. 
%
% smoothPercentileWidth - Default is 31. Width of moving mean window to 
% smooth spatial climatology and threshold.

paramNames = {'Threshold','WindowHalfWidth','smoothPercentileWidth'};
defaults   = {0.9,5,31};

[vThreshold,vWindowHalfWidth,vsmoothPercentileWidth]...
                                       = internal.stats.parseArgs(paramNames, defaults, varargin{:});
                                   
%%  WindowHalfWidth supplementation

temp_clim = cat(3,NaN(size(temp,1),size(temp,2),vWindowHalfWidth), ...
    temp,NaN(size(temp,1),size(temp,2),vWindowHalfWidth));
                                  
%% Calculating climatology and thresholds

date_true=datevec(cli_start-vWindowHalfWidth:cli_end+vWindowHalfWidth); 
date_true=date_true(:,1:3);

date_false = date_true;
date_false(:,1) = 2012; 

fake_doy = day(datetime(date_false),'dayofyear'); 
ind = 1:length(date_false); 

mclim=NaN(size(temp,1),size(temp,2),366);
m90=NaN(size(temp,1),size(temp,2),366);

for i=1:366
    if i == 60
        
    else
        ind_fake=ind;
        ind_fake(fake_doy==i & ~ismember(datenum(date_true),cli_start:cli_end))=nan;
       
        m90(:,:,i) = quantile(temp_clim(:,:,any(ind_fake'>=(ind_fake(fake_doy == i)-vWindowHalfWidth) & ind_fake' <= (ind_fake(fake_doy ==i)+vWindowHalfWidth),2)),vThreshold,3);
        mclim(:,:,i) = mean(temp_clim(:,:,any(ind_fake'>=(ind_fake(fake_doy == i)-vWindowHalfWidth) & ind_fake' <= (ind_fake(fake_doy ==i)+vWindowHalfWidth),2)),3,'omitnan');

    end
end

% -------------------------------- Dealing with Feb29 ----------------------------------
m90(:,:,60) = mean(m90(:,:,[59 61]),3,'omitnan');
mclim(:,:,60) = mean(mclim(:,:,[59 61]),3,'omitnan');

% -------------------------------- Smoothing ---------------------------------------------
m90long = smoothdata(cat(3,m90,m90,m90),3,'movmean',vsmoothPercentileWidth); 
m90 = m90long(:,:,367:367+365);
mclimlong = smoothdata(cat(3,mclim,mclim,mclim),3,'movmean',vsmoothPercentileWidth);
mclim = mclimlong(:,:,367:367+365);

end