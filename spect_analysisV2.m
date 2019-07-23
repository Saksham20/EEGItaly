
% this code finds the epoc specific time independent spectrograms using
% pmtm()
%% 
nwues=length(AllEEGData_complete_reshaped);

%Analyzing complete spectrogram(not temporal) for each of the sub epocs of every trial: 
%which frequencies get modulated max during the sub epocs will be analyzed.
%based on these freqs then will the spectral purturbation be analyzed. 
events_onoffsetkin_data(2).sessno(2).data(5,11)=events_onoffsetkin_data(2).sessno(2).data(5,11)-1; % for some reason only the 8th traial in 2 wue and 2nd sesion is one shorter in Mara_reshaped 

for i=currwue
    nsess=length(AllEEGData_complete_reshaped(i).sessno);
    max_subepoc1=0;   
    for j=1:nsess
        
        ntrials=length(AllEEGData_complete_reshaped(i).sessno(j).trial);          
        max_subepoc=diff(events_onoffsetkin_data(i).sessno(j).data);
        max_subepoc=max(max(max_subepoc));
        if max_subepoc>max_subepoc1
           max_subepoc1=max_subepoc;
        end        
    end
    numchans=size(AllEEGData_complete_reshaped(i).sessno(1).trial(1).data,1);
    specOutmat=NaN*ones(numchans,max_subepoc1,4,ntrials,nsess);
    datOutmat=zeros(size(specOutmat));
    
    for j=1:nsess
        ntrials=length(AllEEGData_complete_reshaped(i).sessno(j).trial); 
        
%         max_subepoc=diff(events_onoffsetkin_data(i).sessno(j).data);
%         max_subepoc=max(max(max_subepoc));
%         specOutmat=zeros(126,max_subepoc,5,ntrials);
        
        for k=1:ntrials
            trialen=size(AllEEGData_complete_reshaped(i).sessno(j).trial(k).data,2);
   
            for m=1:4; %looping over sub epocs
                startime=events_onoffsetkin_data(i).sessno(j).data(m,k);
                endtime=events_onoffsetkin_data(i).sessno(j).data(m+1,k);
                
                for l=1:numchans
                   dat=AllEEGData_complete(i).sessno(j).data(l,startime:endtime);
                   specOut=pmtm(dat);
                   querypnts=linspace(1,length(dat),max_subepoc1); %interpolating to same length as max subepoc
                   interpdSpec=interp1(specOut,querypnts,'PCHIP');
                   %interpddat=interp1(dat,querypnts,'PCHIP');
                   specOutmat(l,:,m,k,j)=interpdSpec; %specOut: chan x freq x epoc x trial x session  
                   %datOutmat(l,:,m,k,j)=interpddat;
                end
            end
            
        end  
        
    end
    spectdata_SUBEPOC(i).data=specOutmat;
    %alldata_SUBEPOC(i).data=datOutmat;
end

%% borrowed
%% Finding sig difference between 4 conditions for each frequencies across wues. 
%                 | All 4 different | Pre/post & Movmt 
% ----------------------------------------------------
% 1)Indi freq     |        x        |      _/
% 2)Band specific |        x        |      _/ 

% 1)indi freq ----------
spectdata_SUBEPOC1=spectdata_SUBEPOC;
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
    sigValtemp=zeros(size(pvaltemp));
    for chans=1:length(goodelectrodes)
        for freq=1:length(goodidx)
            dat1=reshape(NData(chans,freq,1,:,:),size(NData,4),size(NData,5));
            dat1=dat1(:);dat1(isnan(dat1))=[]; 
            dat2=reshape(NData(chans,freq,2,:,:),size(NData,4),size(NData,5));
            dat2=dat2(:);dat2(isnan(dat2))=[];          
            pvaltemp(chans,freq)=anova1([dat1,dat2],[1 2],'off');
            sigValtemp(chans,freq)=(pvaltemp(chans,freq)<0.05);
        end
    end
    
    pval(wueno).data=pvaltemp;
    sigVal(wueno).data=sigValtemp;
    
end
clearvars dat1 dat2 pvaltemp NData spectdata_SUBEPOC1

%% NOrmalizing the spectdata subepoc for every trial:

spectdata_SUBEPOC_meansub=spectdata_SUBEPOC;
for wueno=currwue
    sz=size(spectdata_SUBEPOC(wueno).data);
    loopsz=prod(sz(3:5));
    for i=1:sz(3)
        for j=1:sz(4)
            for k=1:sz(5)
                mn=nanmean(spectdata_SUBEPOC(wueno).data(:,:,i,j,k),2);
                stde=nanstd(spectdata_SUBEPOC(wueno).data(:,:,i,j,k),0,2);
                mn=repmat(mn,1,size(spectdata_SUBEPOC(wueno).data,2));
                stde=repmat(stde,1,size(spectdata_SUBEPOC(wueno).data,2));
                spectdata_SUBEPOC_meansub(wueno).data(:,:,i,j,k)=...
                    (spectdata_SUBEPOC(wueno).data(:,:,i,j,k)-mn);%./stde;
            end
        end
    end  
end
clearvars mn stde loopsz sz

%% converting the spectdata_subepoc to db : 
for i=currwue
    spectsz=size(spectdata_SUBEPOC(i).data);
    spectsz=reshape(spectdata_SUBEPOC(i).data,spectsz(1),prod(spectsz(2:end)));
    maxval=zeros(spectsz(1),prod(spectsz(3:end)));
    currdatadB=zeros(size(temp));
    
    for j=1:prod(spectsz(3:end))
        currdata=temp(:, ((j-1)*spectsz(2)+1):((j)*spectsz(2)-1));
        maxval(:,j)=max(currdata,[],2);
        currdatadB(:,((j-1)*size(currdata,2)+1):((j)*size(currdata,2)))=...
            10*log10(bsxfun(@rdivide,currdata,maxval(:,j)));
    end
    spectdata_SUBEPOC_dB(i).data=reshape(currdatadB,spectsz);
    
end
clearvars spectsz spectsz maxval currdatadB currdata

%% Channel selection based on a specified freq band variation: 

% band=[9;11];% taking the beta band : 12:28 hZ
% 
% for i=currwue %loop over all wues: 
%         freqz=size(spectdata_SUBEPOC(i).data,2);
%         freqz=[1:freqz]*500/freqz; 
%         freqz_idx=dsearchn(freqz',band);
%         meanVariation=[];
%         thresh=0.50; % some threshold for percentage variation
%         for l=1:numchans;
%             
%             datnow=mean(mean(spectdata_SUBEPOC(i).data(l,freqz_idx(1):freqz_idx(2),:,:,:),4),5);% 1 x freqz x 5 matrix 
%             
%             
%             bandActivity=mean(datnow(:,:,2:4,:,:),3); % mean of reaching, stop and back 
%             baseActivity=mean(datnow(:,:,[1 5],:,:),3); % mean of the pre and post activity as baseline
%             
%             Variation=(bandActivity-baseActivity)./baseActivity;
%             meanVariation(l,1)=mean(Variation);
%             
%             
%         end
%         ActiveChans(i).chans= find(meanVariation>thresh);
%         
%         
% end




