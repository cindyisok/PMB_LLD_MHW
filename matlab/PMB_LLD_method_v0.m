% This is for PMB-LLD method
% Ver.0
% update date: Jun 3, 2022

start_yr = 1982;
end_yr = 2014;
mov_base = 31;
length_yr = end_yr - start_yr + 1 - (mov_base-1) - mov_base; % whole length minus the first and last 15 yr and the baseline

sst_all = [];

for m = 1:length_yr
    
    sst = [];
    
    for k = m:m+mov_base-1
        %%%%%%%%%%%%%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load(['sst_' num2str(k+start_yr-1) '.mat']) % the data now should be sst_1.mat, ..., sst_161.mat
        eval(['sst=cat(3,sst,sst_' num2str(k+start_yr-1) ');']) % concatenate in the time dimension
    
        %%%%%%%%%%%%%%%%%%%%%%%%%% detrending %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i = 1:size(sst,1) % lon
            for j = 1:size(sst,2) % lat
                a = detrend(squeeze(sst(i,j,:)),1); % detrend.m in matlab minus the linear function, not the linear trend ...
                                                    % so a is around 0 
                b = squeeze(sst(i,j,:)); % the raw sst
                sst(i,j,:) = a + b(1) - a(1); % b(1)-a(1) is the intercept in y-axis
            end
        end
        disp(k)
    end
    
    % extract the centered detrend year 
    sst_all = cat(3,sst_all,sst(:,:, ...
                                    datenum( start_yr + (mov_base-1)/2 + mov_base + m-1, 1, 1 ) - datenum( start_yr,1,1 ) + 1: ...
                                    datenum( start_yr + (mov_base-1)/2 + mov_base + m-1,12,31 ) - datenum( start_yr,1,1 ) + 1  ));
                                    % e.g. m = 1: 1997.1.1 - 1997.12.31
    disp(m)
end

% calculate theta and climatology using year 1-31

% detect MHWs using year 47-146