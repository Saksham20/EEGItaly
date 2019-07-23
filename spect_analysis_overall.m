
% this code finds the overall spectral peaks for each of the patients. 
% it uses pwelch
%% extraction and averaging: 
%currentdata=AllEEGData_MARAcleaned;
currentdata=AllEEGData_complete;

if isfield(currentdata(1).sessno,'eegData')
    for i=1:6
        for j=1:length(currentdata(i).sessno)
                currentdata(i).sessno(j).data=currentdata(i).sessno(j).eegData;
        end
    end
end

nowues=length(currentdata);

for i=currwue
    
    nosess=length(currentdata(i).sessno);totlen=[];
    for j=1:nosess
        totlen=[totlen,size(currentdata(i).sessno(j).data,2)];
    end
    [minlen,idx]=min(totlen);
    nchan=size(currentdata(i).sessno(1).data,1);
    pwelchlen=(2^ceil((log(floor(minlen/4.5))/log(2))))/2+1; % pwelch calculation: check pwelch doc. 
    spectog=zeros(nchan,pwelchlen,nosess); 
    
    for j=1:nosess
        spectog(:,:,j)=(pwelch((currentdata(i).sessno(j).data(:,1:minlen))'))';
    end
    
    avgspectog=mean(spectog,3);%mean for all sessions
    globalmean=(mean(avgspectog(:,[1:floor(0.14*pwelchlen)]),2));% avg pwr in 0-70 hz 
    globalmean=repmat(globalmean,1,floor(0.14*pwelchlen));
    temp=(avgspectog(:,[1:floor(0.14*pwelchlen)])-globalmean)';
    avgspectogsm=zeros(size(temp));
    for g=1:nchan
        avgspectogsm(:,g)=smooth(temp(:,g),70,'moving');
    end
    
    spectog_wue(i).rawdata=spectog;
    spectog_wue(i).avgdata=avgspectog(:,[1:floor(0.14*pwelchlen)])-globalmean;
    spectog_wue(i).smavgdata=avgspectogsm';
    spectog_wue(i).freqs=[1:0.14*pwelchlen].*(500/pwelchlen);
    spectog_wue(i).globlmean=globalmean;
end
spectog_wue_filtered_sasica=spectog_wue;



