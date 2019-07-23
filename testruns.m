wue=1;
ses=1;
trlno=2;
close all;
trlmean=AllEEGfiltered_trial_mean(wue).sessno(ses).trial(trlno).eegdata;
trlstd=std(AllEEGfiltered_reshaped_avgref(wue).sessno(ses).trial(trlno).eegdata);
labellist={}; 
%chnos=[11:14, 17:26]; col=hsv(length(chnos));count=0;meanchanval=zeros(length(chnos),2);
chnos=1:126;col=hsv(length(chnos));count=0;meanchanval=zeros(length(chnos),2);

for h=chnos
    count=count+1;
    labellist=[labellist,{['ch ' (num2str(h))]}];
    meanchanval(count,:)=[mean(AllEEGfiltered_reshaped_avgref(wue).sessno(ses).trial(trlno).eegdata(h,:))...
        std(AllEEGfiltered_reshaped_avgref(wue).sessno(ses).trial(trlno).eegdata(h,:))];
    AllEEGfiltered_trial_zscored(wue).sessno(ses).trial(trlno).eegdata(h,:)=...
        AllEEGfiltered_reshaped_avgref(wue).sessno(ses).trial(trlno).eegdata(h,:)./trlstd;

end; 
%plot(trlmean,'k');hold off;
%legend(labellist);

figure('Name','spectogram');
spectrogram(AllEEGfiltered_reshaped_prepoc_avgref(wue).sessno(ses).trial(trlno).eegdata(20,:),128,120,128,1e3);

figure('Name','zscored');count=0;
for i=chnos
    count=count+1;
    plot(AllEEGfiltered_trial_zscored(wue).sessno(ses).trial(trlno).eegdata(count,:),'Color',col(count,:));hold on;
end
figure('Name','avgplot');count=0;
for i=chnos
    count=count+1;
    %plot(( AllEEGfiltered_reshaped_avgref(wue).sessno(ses).trial(trlno).eegdata(count,:)),'Color',col(count,:));hold on;
    plot(( AllEEGData_cut_practical(wue).sessno(ses).eegData(count,:)),'Color',col(count,:));hold on;
end    

figure('Name','meanchanvals')
plot(meanchanval(:,1),'r*');hold on;plot(meanchanval(:,2),'b+');legend({'mean','std'});