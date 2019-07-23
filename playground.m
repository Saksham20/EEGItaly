% 
% for as=2:10
%     eegcleantrial(as).data=cleaneeg(:,(onoffsettimes(1,as)-500):onoffsettimes(2,as));
%     
% end
% 
% chlist=[18:26];
% chlist=[1:120];
% chlist=20;
% count=0;
% colo=hsv(length(chlist));figure('Name','trial plots');
% for ch=chlist
%     count=count+1;
%     plot(eegcleantrial(5).data(ch,:),'Color',colo(count,:));hold on;
% end




% %% --------- extracting inter trial time across all wues 
% intertrltm=[];starttm=[];
% for i=1:6
%     nosess=length(kin_data(i).sessno);
%     for j=1:nosess 
%         onset=kin_data(i).sessno(j).data.Onset_event;
%         offset=kin_data(i).sessno(j).data.Offset_event;
%         diffs=onset(2:end)-offset(1:end-1);
%         intertrltm=[intertrltm;[repmat([i,j],length(diffs),1),diffs(:)]];
%         starttm=[starttm;[[i,j],kin_data(i).sessno(j).data.Onset_event(1,1)]];
%     end
% end
% 
% % checking starttime is how many seconds after the end of the
% % TENS spikes so as to get a pre onset time reading. 
% startvals=xlsread('C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\TimeOffsetData.xlsx','B28:E33');
% startvals=startvals'; 
% startvals=startvals(:);
% startvals(startvals==0)=[];
% diffsrtval=-startvals+(starttm(:,3))*1000;
% find(diffsrtval<1750)
% 
% matout=[];
% for i=1:6
%     nosess=length(kin_data(i).sessno);
%     for j=1:nosess 
%         matout=[matout;i,j,any(double(new_onsetimes(i).sessno(j).times(1:end-1))-double(events_onoffsetkin_data(i).sessno(j).data(1,:)))];    
%     end
% end




% % checking endtime (TENS start) is how many seconds after the end of the
% % last trial so as to get some time from which a baseline can be extracted  

% endvals=xlsread('C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\TimeOffsetData.xlsx','H28:K33');
% endvals=endvals'; 
% endvals=endvals(:);
% endvals(endvals==0)=[];
% % diffsrtval=-startvals+(starttm(:,3))*1000;
% % find(diffsrtval<1750)
% 
% matout1=[];c=0;
% for i=1:6
%     nosess=length(kin_data(i).sessno);
%     for j=1:nosess
%         c=c+1;
%         tend=int32(1000*kin_data(i).sessno(j).data.Offset_event(end));
%         matout1=[matout1;i,j,(endvals(c)-tend)];    
%     end
% end

% %% 2) band specific ------------
% alpha=[8 12];beta=[13 30];
% bands=[alpha; beta];
% 
% for wueno=1:1
%     %--- taking the mean of spectrograms of the 3 movement epocs
%     sz=size(spectdata_SUBEPOC(wueno).data);
%     sz(1,3)=sz(1,3)-2;
%     NData=zeros(sz);
%     NData(:,:,2,:,:)=mean(spectdata_SUBEPOC(wueno).data(:,:,2:end,:,:),3);
%     NData(:,:,1,:,:)=(spectdata_SUBEPOC(wueno).data(:,:,1,:,:));
%     pvalbandtemp=zeros(2,size(NData,2));
%     sigValbandtemp=zeros(size(pvalbandtemp));
%     freqs=size(spectdata_SUBEPOC(wueno).data,2);
%     freqs=[1:freqs].*500/freqs;
%     bandid=arrayfun(@(x) dsearchn(freqs',x),bands);
%    
%     for chans=1:length(goodelectrodes)
%         for freq=1:size(bands,2)
%             dat1=mean(NData(chans,bandid(freq,1):bandid(freq,2),1,:,:),2);
%             dat1=reshape(dat1(chans,:,1,:,:),size(NData,4),size(NData,5));
%             dat1=dat1(:);
%             dat2=mean(NData(chans,bandid(freq,1):bandid(freq,2),2,:,:),2);
%             dat2=reshape(dat2(chans,:,2,:,:),size(NData,4),size(NData,5));
%             dat2=dat2(:);          
%             pvalbandtemp(chans,freq)=anova1([dat1,dat2],[1 2],'off');% 2 columns. 1st is alpha, 2nd beta. 
%             sigValbandtemp(chans,freq)=(pvalbandtemp(chans,freq)<0.05);
%         end
%     end
%     pval_band(wueno).data=pvalbandtemp;
%     sigVal_band(wueno).data=sigValbandtemp;
% end
% 
% 

%% 
% for i=1:6
%     spectsz=size(spectdata_SUBEPOC(i).data);
%     temp=reshape(spectdata_SUBEPOC(i).data,spectsz(1),prod(spectsz(2:end)));
%     maxval=zeros(spectsz(1),prod(spectsz(2:end)));
%     currdatadB=zeros(size(temp));
%     
%     for j=1:prod(spectsz(3:end))
%         currdata=temp(:, ((j-1)*spectsz(2)+1):((j)*spectsz(2)-1));
%         maxval(:,j)=max(currdata,[],2);
%         currdatadB(:,((j-1)*size(currdata,2)+1):((j)*size(currdata,2)))=...
%             10*log10(bsxfun(@rdivide,currdata,maxval(:,j)));
%     end
%     spectdata_SUBEPOC_dB(i).data=reshape(currdatadB,spectsz);
%     
% end

%% diffstraval modification to have sess and trial info in the array
% counter=0;
% diffstrvalappend=zeros(size(diffsrtval,2),2);
% for i=1:6
%     leng=length(TimeDelay(i).Tdiff);      
%     for j=1:leng
%         counter=counter+1;
%         diffstrvalappend(counter,:)=[i,j];
%     end    
% end
% diffsrtval=[diffstrvalappend, diffsrtval];       
        
        
%% 2) band specific ------------
alpha=[8 12];beta=[13 30];
bands=[alpha; beta];
spectdata_SUBEPOC1=spectdata_SUBEPOC_norm;
for wueno=currwue
    %--- taking the mean of spectrograms of the 3 movement epocs
    
    sz=size(spectdata_SUBEPOC1(wueno).data);
    sz(1,3)=sz(1,3)-2;% since one is rest and other is movement 
    NData=zeros(sz);
    NData(:,:,2,:,:)=nanmean(spectdata_SUBEPOC1(wueno).data(:,:,[2 3 4],:,:),3);% avg of last 3 epocs of movement
    NData(:,:,1,:,:)=(spectdata_SUBEPOC1(wueno).data(:,:,1,:,:));
    pvalbandtemp=zeros(sz(1),size(bands,2));pvaltemp_anova=pvalbandtemp;
    sigValbandtemp=zeros(size(pvalbandtemp));
    sigValtemp_anova=sigValbandtemp;
    moreORless=zeros(size(pvalbandtemp));
    freqs=size(spectdata_SUBEPOC1(wueno).data,2);
    freqs=[1:freqs].*500/freqs;
    bandid=arrayfun(@(x) dsearchn(freqs',x),bands);
     
    for chans=1:sz(1)
        for freq=1:size(bands,2)
            dat1=nanmean(NData(chans,bandid(freq,1):bandid(freq,2),1,:,:),2);
            dat1=reshape(dat1,size(NData,4),size(NData,5));
            dat1=dat1(:);dat1(isnan(dat1))=[]; 
            dat2=nanmean(NData(chans,bandid(freq,1):bandid(freq,2),2,:,:),2);
            dat2=reshape(dat2,size(NData,4),size(NData,5));
            dat2=dat2(:);dat2(isnan(dat2))=[];  
            %[~,pvalbandtemp(chans,freq)]=ttest(dat1,dat2);
            pvalbandtemp(chans,freq)=permutation_test(dat1,dat2);
            pvaltemp_anova(chans,freq)=anova1([dat1,dat2],[1 2],'off');
            sigValbandtemp(chans,freq)=(pvalbandtemp(chans,freq)<0.05);
            sigValtemp_anova(chans,freq)=(pvaltemp_anova(chans,freq)<0.05);
%             pvalbandtemp(chans,freq)=anova1([dat1,dat2],[1 2],'off');% 2 columns. 1st is alpha, 2nd beta. 
%             sigValbandtemp(chans,freq)=(pvalbandtemp(chans,freq)<0.05);
            moreORless(chans,freq)=(nanmean(dat1)>nanmean(dat2));% baseline should have more power when significant 
        end
    end
    pval_band_norm_anova(wueno).data=[pvaltemp_anova goodelectrodes' moreORless];
    sigVal_band_norm_anova(wueno).data=[sigValtemp_anova goodelectrodes' moreORless];
    pval_band_norm(wueno).data=[pvalbandtemp goodelectrodes' moreORless];
    sigVal_band_norm(wueno).data=[sigValbandtemp goodelectrodes' moreORless];
end
clearvars NData pvalbandtemp sigValbandtemp freqs dat1 dat2 pvalbandtemp sigValbandtemp moreORless bandid



  
        