%% borrowed
function [pval sigVal]=signifreq(spectdata_SUBEPOC1);
%% Finding sig difference between 4 conditions for each frequencies across wues. 
%                 | All 4 different | Pre/post & Movmt 
% ----------------------------------------------------
% 1)Indi freq     |        x        |      _/
% 2)Band specific |        x        |      _/ 

goodelectrodes=[1:126];
currwue=1:6;
freqcutoffout=100;freqcutoffin=1;

% 1)indi freq ----------
for wueno=currwue
    %--- taking the mean of spectrograms of the 3 movement epocs
    freqlen=size(spectdata_SUBEPOC1(wueno).data,2);
    freqlen=[1:freqlen].*500/freqlen;
    goodidx=find(freqlen<freqcutoffout);
    
    sz=size(spectdata_SUBEPOC1(wueno).data);
    sz(1,3)=sz(1,3)-2;
    NData=zeros(sz);
    NData(:,:,2,:,:)=nanmean(spectdata_SUBEPOC1(wueno).data(:,:,2:end,:,:),3);
    NData(:,:,1,:,:)=(spectdata_SUBEPOC1(wueno).data(:,:,1,:,:));
    pvaltemp=zeros(length(goodelectrodes),length(goodidx));
    pvaltemp_anova=pvaltemp;
    sigValtemp=zeros(size(pvaltemp));
    sigValtemp_anova=sigValtemp;
    for chans=1:length(goodelectrodes)
        for freq=1:length(goodidx)
            dat1=reshape(NData(chans,freq,1,:,:),size(NData,4),size(NData,5));
            dat1=dat1(:);dat1(isnan(dat1))=[]; 
            dat2=reshape(NData(chans,freq,2,:,:),size(NData,4),size(NData,5));
            dat2=dat2(:);dat2(isnan(dat2))=[];  
            pvaltemp(chans,freq)=permutation_test(dat1,dat2);
            pvaltemp_anova(chans,freq)=anova1([dat1,dat2],[1 2],'off');
            sigValtemp(chans,freq)=(pvaltemp(chans,freq)<0.05);
            sigValtemp_anova(chans,freq)=(pvaltemp_anova(chans,freq)<0.05);
        end
    end
    
    pval(wueno).data=pvaltemp;
    sigVal(wueno).data=sigValtemp;
    
end
end