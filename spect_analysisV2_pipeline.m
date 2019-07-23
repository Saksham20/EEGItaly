
% this code finds the epoc specific time independent spectrograms using
% pmtm()
%% 

nwues=length(AllEEGData_complete_reshaped);
goodelectrodes=[1:126];
currwue=1:6;
freqcutoffout=100;freqcutoffin=1;

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
    numchans=size(AllEEGData_complete(3).sessno(1).data,1);
    specOutmat=NaN*ones(numchans,max_subepoc1,4,ntrials,nsess);
    datOutmat=zeros(size(specOutmat));
    
    for j=1:nsess
        ntrials=length(AllEEGData_complete_reshaped(i).sessno(j).trial);         
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
                   specOutmat(l,:,m,k,j)=interpdSpec; %specOut: chan x freq x epoc x trial x session  
                end
            end
            
        end  
        
    end
    spectdata_SUBEPOC(i).data=specOutmat;

end


[pval sigVal]=signifreq(spectdata_SUBEPOC);


%% NOrmalizing the spectdata subepoc for every trial:

spectdata_SUBEPOC_norm=spectdata_SUBEPOC;
for wueno=currwue
    sz=size(spectdata_SUBEPOC(wueno).data);
    loopsz=prod(sz(3:5));
    freqlen=size(spectdata_SUBEPOC(wueno).data,2);
    freqlen=[1:freqlen].*500/freqlen;
    goodidx=find(freqlen<freqcutoffout);
    for i=1:sz(3)
        for j=1:sz(4)
            for k=1:sz(5)
                mn=nanmean(spectdata_SUBEPOC(wueno).data(:,goodidx,i,j,k),2);
                stde=nanstd(spectdata_SUBEPOC(wueno).data(:,goodidx,i,j,k),0,2);
                mn=repmat(mn,1,length(goodidx));
                stde=repmat(stde,1,length(goodidx));
                spectdata_SUBEPOC_norm(wueno).data(:,goodidx,i,j,k)=...
                    (spectdata_SUBEPOC(wueno).data(:,goodidx,i,j,k)-mn)./stde;
            end
        end
    end  
end
clearvars mn stde loopsz sz
[pval_norm sigVal_norm]=signifreq(spectdata_SUBEPOC_norm);


%% conversion to dB scale for every trial: 
spectdata_SUBEPOC_dB=spectdata_SUBEPOC;
for wueno=currwue
    sz=size(spectdata_SUBEPOC(wueno).data);
    loopsz=prod(sz(3:5));
    freqlen=size(spectdata_SUBEPOC(wueno).data,2);
    freqlen=[1:freqlen].*500/freqlen;
    goodidx=find(freqlen<freqcutoffout);
    for i=1:sz(3)
        for j=1:sz(4)
            for k=1:sz(5)
                mx=max(spectdata_SUBEPOC(wueno).data(:,goodidx,i,j,k),[],2);
                
                mx=repmat(mx,1,length(goodidx));
               
                spectdata_SUBEPOC_dB(wueno).data(:,goodidx,i,j,k)=...
                    10*log10((spectdata_SUBEPOC(wueno).data(:,goodidx,i,j,k))./mx);
            end
        end
    end  
end
clearvars mx loopsz sz
[pval_dB sigVal_dB]=signifreq(spectdata_SUBEPOC_dB);


%% plotting: 
imgstoreloc1='C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\scripts\sasicaout\sasica\wicaOn\pipeout\';
imgstoreloc_norm='C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\scripts\sasicaout\sasica\wicaOn\pipeout\norm\';
imgstoreloc_dB='C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\scripts\sasicaout\sasica\wicaOn\pipeout\dB\';

imgst={imgstoreloc1, imgstoreloc_norm,imgstoreloc_dB };

for i=1:3
imgstoreloc=imgst{1,i};

%---------
if i==1
    spectdata_SUBEPOC1=spectdata_SUBEPOC;sigval=sigVal;tag='';
elseif i==2
    spectdata_SUBEPOC1=spectdata_SUBEPOC_norm;sigval=sigVal_norm;tag='normalized,';
elseif i==3
    spectdata_SUBEPOC1=spectdata_SUBEPOC_dB;sigval=sigVal_dB;tag='dB,';
end


for wueno=currwue
    freqlen=size(spectdata_SUBEPOC1(wueno).data,2);
        freqlen=[1:freqlen].*500/freqlen;
        goodidx=find(freqlen<freqcutoffout);
        if ~mkdir([imgstoreloc sprintf('wueno%d',wueno)]);
            mkdir([imgstoreloc sprintf('wueno%d',wueno)]);
        end
        
    for chnonew=1:(size(spectdata_SUBEPOC1(wueno).data,1))
        figure();
        col=hsv(4);
        
        %goodelectrodes=[11:36 53:58 76:77 122:123 125:126];
        %totsesstrials=(size(spectdata_SUBEPOC(wueno).data,4))*(size(spectdata_SUBEPOC(wueno).data,5));
        for t=1:4; 
            plot(freqlen(goodidx),nanmean(nanmean(spectdata_SUBEPOC1(wueno).data(chnonew,(goodidx),t,:,:),4),5),...
                'Color',col(t,:)); hold on; 
        end; 
        ax=gca;yvl=ax.YTick;yvl=yvl(2);
        sigdat=sigval(wueno).data(chnonew,:);
        plot(freqlen(goodidx),yvl*sigdat(goodidx),...
                '-*');
        title(sprintf('wue no. %d, chno. %d (%s avg across trial/sess)',wueno,goodelectrodes(chnonew),tag));
        legend({'pre','reach','grasp','reverse'});grid on
        xlabel(['Hz']);ylabel(['PSD avg across sessions']);
        print([imgstoreloc sprintf('wueno%d',wueno) filesep...
            sprintf('chno_%d ',goodelectrodes(chnonew))],'-djpeg');close all;
    end
end

end






