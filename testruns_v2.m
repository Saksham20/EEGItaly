wue=1;
ses=2;
close all;
goodelectrodes=[11:36 53:58 76:77 122:123 125:126];

chnos=1:length(goodelectrodes);col=hsv(length(chnos));


figure('Name','rawdata');count=0;
for i=chnos
    count=count+1;
    plot(( AllEEGData_cut_practical(wue).sessno(ses).eegData(count,:)),'Color',col(count,:));hold on;
end  

figure('Name','filtered');count=0;
for i=chnos
    count=count+1;
    plot(( AllEEGData_cut_practical_filtered(wue).sessno(ses).eegData(count,:)),'Color',col(count,:));hold on;
end  

nointerps=1;
eegdatin=eegdata_struct_gen(AllEEGData_cut_practical(wue).sessno(ses).eegData,goodelectrodes);chrem=[];idxtot=[];
for i1=1:nointerps
    
    [interpeeg,idx] = pop_rejchan(eegdatin, 'elec',[1:eegdatin.nbchan],'threshold',[-3 3],'norm','on','measure','spec','freqrange',[1 70]);
    %interpeeg = eeg_interp(eegdatin, idx, 'spherical');
    chrem=[chrem;[i1*ones(length(idx),1) idx(:)]];    
    idxtot=[idxtot;idx(:)];
             
    figure('Name',sprintf('interpolated %d',i1));count=0;
    for i=1:interpeeg.nbchan
        count=count+1;
        plot(( interpeeg.data(count,:)),'Color',col(count,:));hold on;
    end  
    eegdatin=interpeeg;
    
end

interpeeg = eeg_interp(eegdatin, unique(idxtot), 'spherical');
figure('Name','interpolated_last');count=0;
for i=chnos
        count=count+1;
        plot(( interpeeg.data(count,:)),'Color',col(count,:));hold on;
end 