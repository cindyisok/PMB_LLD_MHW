% This is PMB-LLD method for users
% Ver.1
% Author: Shengpeng Wang and Di Sun
% update date: Jun 10, 2022

% =================================================================================================  
%    This is the call function of PMB-LLD Method for users. It applies to
%    grided datasets.
%  
%    Input part is for users to input the start_yr, end_yr of the whole timeseries, as well 
%    as the moving baseline period.
%    
%    Theta and climatology detrending part is to detrend linear trend of the 31-yr fixed SST
%    for theta and climatology calculation.
%    
%    SST detrending part is to detrend linear trend of moving 31-yr SST for MHW detecting.
%    
%    MHW detecting part firstly calculates the theta and climatology of the fixed 31-yr and
%    then detect MHWs using detrended SST.
% ================================================================================================

%%
% --------------------------------------- Input -----------------------------------------
%  Please input the start year, end year of the whole timeseries ...
%  (including both theta calculating and MHW detetcing) and the moving baseline period

start_yr  = 1901;
end_yr    = 2100;
mov_base  = 31; % moving baseline period
mov_base2 = (mov_base-1)/2; % half of the moving baseline period


% ----------------------------- Theta and climatology detrending -----------------------------
% This part is for a fixed 31-year baseline period for the detrending of seasonally varying ...
% theta and climatology

% Load data of the first 31yr
sst_theta = []; % SST used to calculate theta and climatology

for k = start_yr:start_yr + mov_base - 1
    load(['/home/jingzhao/sundi/data/SST_NA_LR/SST_year_' num2str(k) '.mat']) % the data now should be sst_1900.mat, ..., sst_2100.mat
    
    % if the leap year is 365 day (e.g in CESM model)
    if ( mod(k,4) == 0 && mod(k,100) ~= 0 ) || mod(k,400) == 0
        sst_leap_day = mean(SST0(:,:,[59 60]),3,'omitnan');
        SST0 = cat(3,SST0(:,:,1:59),sst_leap_day,SST0(:,:,60:end));
    end
    
    sst_theta=cat(3,sst_theta,SST0); % concatenate in the time dimension
end

% This is for extracting the seasonal cycle
[mclim,~] = cal_theta(sst_theta,datenum(start_yr,1,1),datenum(start_yr+mov_base-1,12,31));
mclim365 = mclim(:,:,[1:59 61:end]);
for k = start_yr:start_yr+mov_base-1
    if ( mod(k,4) == 0 && mod(k,100) ~= 0 ) || mod(k,400) == 0
         sst_theta(:,:,datenum(k,1,1)-datenum(start_yr,1,1)+1:datenum(k,12,31)-datenum(start_yr,1,1)+1) = ...
            sst_theta(:,:,datenum(k,1,1)-datenum(start_yr,1,1)+1:datenum(k,12,31)-datenum(start_yr,1,1)+1) - ...
            mclim;
    else
        sst_theta(:,:,datenum(k,1,1)-datenum(start_yr,1,1)+1:datenum(k,12,31)-datenum(start_yr,1,1)+1) = ...
            sst_theta(:,:,datenum(k,1,1)-datenum(start_yr,1,1)+1:datenum(k,12,31)-datenum(start_yr,1,1)+1) - ...
            mclim365;
    end
end

% Linear detrending
for i = 1:size(sst_theta,1) % lon
    for j = 1:size(sst_theta,2) % lat
        sst_theta(i,j,:) = detrend(squeeze(sst_theta(i,j,:)),1); 
        % detrend.m in matlab minus the linear function, not the linear trend
    end
end

% This is to add the seasonal variation
% for k = start_yr:start_yr+mov_base-1
%     if ( mod(k,4) == 0 && mod(k,100) ~= 0 ) || mod(k,400) == 0
%          sst_theta(:,:,datenum(k,1,1)-datenum(start_yr,1,1)+1:datenum(k,12,31)-datenum(start_yr,1,1)+1) = ...
%             sst_theta(:,:,datenum(k,1,1)-datenum(start_yr,1,1)+1:datenum(k,12,31)-datenum(start_yr,1,1)+1) + ...
%             mclim;
%     else
%         sst_theta(:,:,datenum(k,1,1)-datenum(start_yr,1,1)+1:datenum(k,12,31)-datenum(start_yr,1,1)+1) = ...
%             sst_theta(:,:,datenum(k,1,1)-datenum(start_yr,1,1)+1:datenum(k,12,31)-datenum(start_yr,1,1)+1) + ...
%             mclim365;
%     end
% end

% ------------------------------------- SST detrending --------------------------------------
% This part is for a moving 31-year baseline period for the detrending of SST used for ...
% MHW detetcing

sst_all = []; % detrended SST used to detect marine heatwaves

for m = start_yr + mov_base2 + mov_base: end_yr - mov_base2 
    % start year excluding the first 15 yr and baseline for theta
    
    % Load data of moving 31yr
    sst = []; % SST for detrending
    for k = m - mov_base2 : m + mov_base2 % m is the centered year, 1947, ..., 2085
        load(['G:/SST_NA_LR/SST_year_' num2str(k) '.mat']) % the data now should be sst_1900.mat, ..., sst_2100.mat

        % if the leap year is 365 day (e.g in CESM model)
        if ( mod(k,4) == 0 && mod(k,100) ~= 0 ) || mod(k,400) == 0
            sst_leap_day = mean(SST0(:,:,[59 60]),3,'omitnan');
            SST0 = cat(3,SST0(:,:,1:59),sst_leap_day,SST0(:,:,60:end));
        end

        sst=cat(3,sst,SST0); % concatenate in the time dimension
    end
    
    [mclim,~] = cal_theta(sst,datenum(m-mov_base2,1,1),datenum(m+mov_base2,12,31));
    mclim365 = mclim(:,:,[1:59 61:end]);
    
% This is for extracting the seasonal cycle
    for k = m - mov_base2 : m + mov_base2
        if ( mod(k,4) == 0 && mod(k,100) ~= 0 ) || mod(k,400) == 0
             sst(:,:,datenum(k,1,1)-datenum(m-mov_base2,1,1)+1:datenum(k,12,31)-datenum(m-mov_base2,1,1)+1) = ...
                sst(:,:,datenum(k,1,1)-datenum(m-mov_base2,1,1)+1:datenum(k,12,31)-datenum(m-mov_base2,1,1)+1) - ...
                mclim;
        else
            sst(:,:,datenum(k,1,1)-datenum(m-mov_base2,1,1)+1:datenum(k,12,31)-datenum(m-mov_base2,1,1)+1) = ...
                sst(:,:,datenum(k,1,1)-datenum(m-mov_base2,1,1)+1:datenum(k,12,31)-datenum(m-mov_base2,1,1)+1) - ...
                mclim365;
        end
    end
    
    % Linear detrending
    for i = 1:size(sst,1)
        for j = 1:size(sst,2)
            sst(i,j,:) = detrend(squeeze(sst(i,j,:)),1);
        end
    end
    
% This is to add the seasonal variation
%     for k = m - mov_base2 : m + mov_base2
%         if ( mod(k,4) == 0 && mod(k,100) ~= 0 ) || mod(k,400) == 0
%              sst(:,:,datenum(k,1,1)-datenum(m-mov_base2,1,1)+1:datenum(k,12,31)-datenum(m-mov_base2,1,1)+1) = ...
%                 sst(:,:,datenum(k,1,1)-datenum(m-mov_base2,1,1)+1:datenum(k,12,31)-datenum(m-mov_base2,1,1)+1) + ...
%                 mclim;
%         else
%             sst(:,:,datenum(k,1,1)-datenum(m-mov_base2,1,1)+1:datenum(k,12,31)-datenum(m-mov_base2,1,1)+1) = ...
%                 sst(:,:,datenum(k,1,1)-datenum(m-mov_base2,1,1)+1:datenum(k,12,31)-datenum(m-mov_base2,1,1)+1) + ...
%                 mclim365;
%         end
%     end
    
    % Extract the centered detrended year 
    sst_all = cat(3,sst_all,sst(:,:,datenum( m, 1, 1 ) - datenum( m-mov_base2,1,1 ) + 1: ...
                                    datenum( m,12,31 ) - datenum( m-mov_base2,1,1 ) + 1  ));
                                    % eg. m = 1947, m-mov_base2 = 1932
    % for testing
    disp(m)
    
end
%%
% --------------------------------------- MHW detecting --------------------------------------- 
% Calculate theta and climatology using year 1901-1931 (31yr)
[mclim,m90] = cal_theta(sst_theta,datenum(start_yr,1,1),datenum(start_yr+mov_base-1,12,31));
save('mclim.mat','mclim','-v7.3')
save('m90.mat','m90','-v7.3')

% Detect MHWs using year 1947-2085
% Transfer mclim and m90 into detect_mhw.m
[MHW,mhw_ts] = detect_mhw(sst_all,mclim,m90,datenum(start_yr+mov_base2+mov_base,1,1),datenum(end_yr-mov_base2,12,31));
save('MHW.mat','MHW','-v7.3')
save('mhw_ts.mat','mhw_ts','-v7.3')
% ============================================= Ending =============================================
%%
% MHW metrics
% metric_used={'Frequency','MeanInt','MaxInt','CumInt','Duration','Days'};
% for i=1:6
%     eval(['[mean_' metric_used{i} ',annual_' metric_used{i} ',trend_' metric_used{i} ',p_' metric_used{i} ']=mean_and_trend(MHW,mhw_ts,1947,' '''' 'Metric' '''' ',' 'metric_used{i}' ');'])
% end
%% Plot
% trends
% trends of annual MHW days
figure('pos',[100 80 800 300],'color',[1 1 1]);
h1 = axes('position',[-0.13 0.23 0.8 0.55]);
m_proj('Equidistant Cylindrical','lon',[120 240],'lat',[20 70]);
m_coast('patch',[.1 .1 .1],'edgecolor','k'); 
m_grid('xtick',120:30:240,'ytick',20:16:70,'fontsize',10)
hold on

trend_Days(p_Days>0.05) = nan;
m_pcolor(lon,lat,trend_Days*100)
hold on
colormap(h1,jet)
ch = colorbar;
set(get(ch,'title'),'string','days/century','fontsize',9);
caxis([-60 60])
title('\it{k_E} \rm{from 1947-2085 : PMB-LLD}','fontsize',15)

% trends of cum intensity
h2 = axes('position',[0.36 0.23 0.8 0.55]);
m_proj('Equidistant Cylindrical','lon',[120 240],'lat',[20 70]);
m_coast('patch',[.1 .1 .1],'edgecolor','k'); 
m_grid('xtick',120:30:240,'ytick',20:16:70,'fontsize',10)
hold on

trend_CumInt(p_CumInt>0.05) = nan;
m_pcolor(lon,lat,trend_CumInt*100)
hold on
colormap(h2,jet)
ch = colorbar;
set(get(ch,'title'),'string','\circC¡¤days/century','fontsize',9);
caxis([-150 150])
title('\it{k_I} \rm{from 1947-2085 : PMB-LLD}','fontsize',15)
% %% mean metrics
% figure('pos',[100 80 800 400],'color',[1 1 1]);
% h1 = subplot(2,2,1);
% m_proj('Equidistant Cylindrical','lon',[120 240],'lat',[20 70]);
% m_coast('patch',[.1 .1 .1],'edgecolor','k'); 
% m_grid('xtick',120:30:240,'ytick',20:16:70,'fontsize',10)
% hold on
% m_pcolor(lon,lat,mean_Frequency)
% hold on
% colormap(jet)
% ch = colorbar;
% u = get(ch,'Position');
% set(ch,'position',u,'ticks',1:0.5:2,'fontsize',10);
% title('Mean frequency','fontsize',12)
% caxis([1 2])
% set(get(ch,'title'),'string','annual every count','fontsize',10);
% 
% h2 = subplot(2,2,2);
% m_proj('Equidistant Cylindrical','lon',[120 240],'lat',[20 max(max(lat))]);
% m_coast('patch',[.1 .1 .1],'edgecolor','k'); 
% m_grid('xtick',120:30:240,'ytick',20:16:70,'fontsize',10)
% hold on
% m_pcolor(lon,lat,mean_MeanInt)
% hold on
% colormap(h2,jet)
% ch = colorbar;
% u = get(ch,'Position');
% set(ch,'position',u,'ticks',0.5:0.5:3,'fontsize',10);
% title('Mean intensity','fontsize',12)
% caxis([0.5 3])
% set(get(ch,'title'),'string','\circC','fontsize',10);
% 
% h3 = subplot(2,2,3);
% m_proj('Equidistant Cylindrical','lon',[120 240],'lat',[20 max(max(lat))]);
% m_coast('patch',[.1 .1 .1],'edgecolor','k'); 
% m_grid('xtick',120:30:240,'ytick',20:16:70,'fontsize',10)
% hold on
% m_pcolor(lon,lat,mean_CumInt)
% hold on
% colormap(h3,jet)
% ch = colorbar;
% u = get(ch,'Position');
% set(ch,'position',u,'ticks',0:20:100,'fontsize',10);
% title('Cumulative intensity','fontsize',12)
% caxis([0 100])
% set(get(ch,'title'),'string','\circC¡¤day','fontsize',10);
% 
% h4 = subplot(2,2,4);
% m_proj('Equidistant Cylindrical','lon',[120 240],'lat',[20 max(max(lat))]);
% m_coast('patch',[.1 .1 .1],'edgecolor','k'); 
% m_grid('xtick',120:30:240,'ytick',20:16:70,'fontsize',10)
% hold on
% 
% m_pcolor(lon,lat,mean_Duration)
% hold on
% colormap(h4,jet)
% ch = colorbar;
% u = get(ch,'Position');
% set(ch,'position',u,'ticks',5:10:45,'fontsize',10);
% title('Mean duration','fontsize',12)
% caxis([5 45])
% set(get(ch,'title'),'string','days','fontsize',10);