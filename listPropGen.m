function list_properties = listPropGen(sessno,trlnum,eeg_chans)

%sessno is the full sessno matrix with all trialnos and the data in it. 

% 1 Median diff value
measure=1;
list_properties(:,measure) = median(diff(sessno.trlno(trlnum).data(eeg_chans,:),[],2),2);
measure = measure + 1;

% 2 Variance of the channels
list_properties(:,measure) = var(sessno.trlno(trlnum).data(eeg_chans,:),[],2);
list_properties(isnan(list_properties(:,measure)),measure)=0;
measure = measure + 1;

% 3 Max difference of each channel
list_properties(:,measure)=(max(sessno.trlno(trlnum).data(eeg_chans,:),[],2)-min(sessno.trlno(trlnum).data(eeg_chans,:),[],2));
measure = measure + 1;

% 4 Deviation from channel mean
sz=sessno.trlno;
sz=length(sz);%gives the number of trials
meantot=zeros(length(eeg_chans),1);totlen=0;
for k=1:sz
    meantot=meantot+sum(sessno.trlno(k).data(eeg_chans,:),2);
    totlen=totlen+size(sessno.trlno(k).data,2);
end
meantot=meantot/totlen;

list_properties(:,measure)=abs(mean(sessno.trlno(trlnum).data(eeg_chans,:),2)-meantot);
measure = measure + 1;

for u = 1:size(list_properties,2)
	list_properties(:,u) = list_properties(:,u) - median(list_properties(:,u));
end


end