%plot the time freq for each : 
imgstoreloc='C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\scripts\sasicaout\sasica\wicaOn\pipeout\timefreq\';
if ~isdir(imgstoreloc)
    mkdir(imgstoreloc);
end

freqband=[13 30]; %beta band 
nFrex = ((diff(freqband)*2)-1);
frex = linspace(freqband(1),freqband(2),nFrex);

for wueno=1:currwue
    nsess=length(AllEEGData_complete_reshaped(wueno).sessno);
    
    min_posttrlen1=Inf;min_pretrlen1=Inf;   
    
    for j=1:nsess
        
        ntrials=length(AllEEGData_complete_reshaped(wueno).sessno(j).trial);
        
            
        posttrlen=(events_onoffsetkin_data(wueno).sessno(sessno).data(end,:))-...
            (events_onoffsetkin_data(wueno).sessno(sessno).data(2,:));
        min_posttrlen=min(posttrlen);        
        if min_posttrlen<min_posttrlen1
           min_posttrlen1=min_posttrlen;
        end 

        pretrlen=(events_onoffsetkin_data(wueno).sessno(sessno).data(2,:))-...
            (events_onoffsetkin_data(wueno).sessno(sessno).data(1,:));
        min_pretrlen=min(pretrlen);
        min_pretrlen=min_pretrlen+1;

        if min_pretrlen<min_pretrlen1
           min_pretrlen1=min_pretrlen;
        end 
            
       
    end
    
    numchans=size(AllEEGData_complete_reshaped(1).sessno(1).trial(1).data,1);
    timeFreq=NaN*ones(numchans,(min_posttrlen1+min_pretrlen1),nFrex,ntrials,nsess);
    
    for i=1:numchans;
        freqtime1=[];count=0;
        for sessno=1:nsess
            ntrials=length(AllEEGData_complete_reshaped(wueno).sessno(sessno).trial);
            trl_movmt_start=diff(events_onoffsetkin_data(wueno).sessno(sessno).data(1:2,:));
            for ntrl=1:ntrials
                count=count+1;
                trl_movmt_start_pt=(trl_movmt_start(ntrl))+2;
                
                data1=AllEEGData_complete_reshaped(wueno).sessno(sessno).trial(ntrl).data(i,:); 
                [ freqtime] = wavelet_morlet( data1,frex,[wueno,sessno,ntrl],events_onoffsetkin_data );
                
%                 querypnts=linspace(1,length(data1),min_trlen1); 
%                 freqtime_interp=zeros(nFrex,min_trlen1);
%                 for u=1:nFrex                      
%                 freqtime_interp(u,:)=interp1(freqtime(u,:),querypnts,'PCHIP');
%                 end
                
                % truncation: 
                leng=size(freqtime,2);
                rightlen=leng-trl_movmt_start_pt+1;
                rightlen_truncate=rightlen-min_posttrlen1;
                leftlen=leng-rightlen;
                leftlen_truncate=leftlen-min_pretrlen1;
                
                freqtime=freqtime(:,(leftlen_truncate+1):(end-rightlen_truncate)); 
                
                
                
                trl_movmt_start_id=dsearchn(int32(querypnts'),trl_movmt_start_pt);
                
                timeFreq(i,:,:,ntrl,sessno)=rot90(rot90(rot90(freqtime)));  
                
%                 print([imgstoreloc sprintf('wueno%d',wueno) filesep...
%                         sprintf('chno_%d ',goodelectrodes(chnonew))],'-djpeg');close all;               
%                 contourf(querypnts,frex,freqtime_interp,40,'linecolor','none'); colorbar;
%                 yran=get(gca,'YLim');
%                 line([trl_movmt_start_pt trl_movmt_start_pt],yran,'Color','k','LineStyle','--');
%                 pause;
            end

        end        
        plot(nanmean(nanmean(nanmean(timeFreq(i,:,:,:,:),3),4),5));title(sprintf('chn no. %d',i));
        yran=get(gca,'YLim');
        line([(min_pretrlen1+1) (min_pretrlen1+1)],yran,'Color','k','LineStyle','--');
        pause;
    end
    
    timeFreq_all(wueno).data=timeFreq;
end