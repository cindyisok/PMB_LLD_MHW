function [MHW,mhw_ts] = detect_mhw(temp,mclim,m90,mhw_start,mhw_end,varargin)

% Description
%
% [MHW,mhw_ts] = detect_mhw(temp,mclim,m90,mhw_start,mhw_end) returns 
% all detected MHW events for the (lon x lat x t) matrix temp from the
% detrended SST by using the moving baseline and also mhw_ts for MHW
% timeseries.

% Input Arguments
%
% temp - the detrended SST (lon x lat x 31year) from the moving 31-yr baseline. 
%
% mclim - cimatology from the fixed 31-yr baseline.
% 
% m90 - threshold from the fixed 31-yr baseline.
%
% mhw_start - datenum(start_yr,1,1) the first day of MHW detecting period.
% mhw_end - datenum(end_yr,12,31) the last day of MHW detecting period.
%
% minDuration - Default is 5. Minimum duration to accept a detection of MHW
% event.
%
% maxGap - Default is 2. Maximum gap accepting joining of MHW events.

paramNames = {'minDuration','maxGap'};
defaults   = {5,2};

[vminDuration,vmaxGap] = internal.stats.parseArgs(paramNames, defaults, varargin{:});

%%
[x_size,y_size]=deal(size(m90,1),size(m90,2));

mbigadd=temp;
date_mhw=datevec(mhw_start:mhw_end);
date_mhw(:,1)=2000;
indextocal = day(datetime(date_mhw),'dayofyear');

ts=str2double(string(datestr(mhw_start:mhw_end,'YYYYmmdd')));

mhw_ts=zeros(x_size,y_size,length(ts));

MHW=[];

%% Detecting MHW in each grid
       
for i=1:x_size
    for j=1:y_size

        mhw_ts(i,j,isnan(squeeze(mbigadd(i,j,:))))=nan;

        if sum(isnan(squeeze(mbigadd(i,j,:))))~=size(mbigadd,3)

            maysum=zeros(size(mbigadd,3),1);

            maysum(squeeze(mbigadd(i,j,:))>=squeeze(m90(i,j,indextocal)))=1;

            trigger=0;
            potential_event=[];

            for n=1:size(maysum,1)
                if trigger==0 && maysum(n)==1
                    start_here=n;
                    trigger=1;
                elseif trigger==1 && maysum(n)==0
                    end_here=n-1;
                    trigger=0;
                    potential_event=[potential_event;[start_here end_here]];
                elseif n==size(maysum,1) && trigger==1 && maysum(n)==1
                    trigger=0;
                    end_here=n;
                    potential_event=[potential_event;[start_here end_here]];
                end
            end

            if ~isempty(potential_event)

                potential_event=potential_event((potential_event(:,2)-potential_event(:,1)+1)>=vminDuration,:);

                if ~isempty(potential_event)

                    gaps=(potential_event(2:end,1) - potential_event(1:(end-1),2) - 1);

                    while min(gaps)<=vmaxGap
                        potential_event(find(gaps<=vmaxGap),2)=potential_event(find(gaps<=vmaxGap)+1,2);
                        loc_should_del=(find(gaps<=vmaxGap)+1);
                        loc_should_del=loc_should_del(~ismember(loc_should_del,find(gaps<=vmaxGap)));
                        potential_event(loc_should_del,:)=[];
                        gaps=(potential_event(2:end,1) - potential_event(1:(end-1),2) - 1);
                    end

                    mhwstart=NaN(size(potential_event,1),1);
                    mhwend=NaN(size(potential_event,1),1);
                    mduration=NaN(size(potential_event,1),1);
                    mhwint_max=NaN(size(potential_event,1),1);
                    mhwint_mean=NaN(size(potential_event,1),1);
                    mhwint_var=NaN(size(potential_event,1),1);
                    mhwint_cum=NaN(size(potential_event,1),1);

                    for le=1:size(potential_event,1)
                        event_here=potential_event(le,:);
                        endtime=ts(event_here(2));
                        starttime=ts(event_here(1));
                        mcl=squeeze(mclim(i,j,indextocal(event_here(1):event_here(2))));
                        mrow=squeeze(mbigadd(i,j,event_here(1):event_here(2)));
                        manom=mrow-mcl;
                        mhw_ts(i,j,event_here(1):event_here(2))=manom;

                        [maxanom,~]=nanmax(squeeze(manom));

                        mhwint_max(le)=...
                            maxanom;
                        mhwint_mean(le)=...
                            mean(manom);
                        mhwint_var(le)=...
                            std(manom);
                        mhwint_cum(le)=...
                            sum(manom);
                        mhwstart(le)=starttime;
                        mhwend(le)=endtime;
                        mduration(le)=event_here(2)-event_here(1)+1;
                    end
                    MHW=[MHW;[mhwstart mhwend mduration mhwint_max mhwint_mean mhwint_var mhwint_cum repmat(i,size(mhwstart,1),1) repmat(j,size(mhwstart,1),1)]];
                end
            end
        end
        
    end
end
        
MHW=table(MHW(:,1),MHW(:,2),MHW(:,3),MHW(:,4),MHW(:,5),MHW(:,6),MHW(:,7),MHW(:,8),MHW(:,9),...
    'variablenames',{'mhw_onset','mhw_end','mhw_dur','int_max','int_mean','int_var','int_cum','xloc','yloc'});