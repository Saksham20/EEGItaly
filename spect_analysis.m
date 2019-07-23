
nwues=length(AllEEGData_complete_reshaped);
twin=500; ovwin=400; nft=500;

for i=1:nwues
    nsess=length(AllEEGData_complete_reshaped(i).sessno);
    for j=1:nsess
        ntrials=length(AllEEGData_complete_reshaped(i).sessno(j).trial);
        [~,maxlentrl_idx]=max(kin_data(i).sessno(j).data.Trial_time(:));
        maxlentrl=size(AllEEGData_complete_reshaped(i).sessno(j).trial(maxlentrl_idx).data,2);
        maxlentime=floor((maxlentrl-ovwin)/(twin-ovwin));
        
        speND=zeros((nft/2+1)+1,... %plus one to have the last row as sequence of tspe 
            maxlentime,...
            126,ntrials);
        
        
        for k=1:ntrials
            trialen=size(AllEEGData_complete_reshaped(i).sessno(j).trial(k).data,2);
            refpntsmax=(trialen-ovwin)/(twin-ovwin);
            
            for l=1:126
               [spe,wspe,tspe,pspe]=spectrogram(AllEEGData_complete_reshaped(i).sessno(j).trial(k).data(l,:),...
                   twin,ovwin,nft,1000,'yaxis','power'); 
               
             
               
               querypnts=linspace(1,length(tspe),maxlentime);
               intempspe=interp1(pspe',querypnts,'PCHIP');
               qrtimepnts=linspace(1,tspe(end),maxlentime);            
               intemptspe=interp1(tspe(:),tspe(:),qrtimepnts);
               speND(:,:,l,k)=[intempspe';intemptspe];
            end
        end
        
        
        spectdata(i).sessno(j).data=speND;
        
    end
end


