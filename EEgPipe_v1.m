wue_no=6;
wuelist={'wue02','wue03','wue06','wue09','wue10','wue11'};
kinloc='C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\KIN events';
load('C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\MAT_files\notrials.mat');
load('C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\MAT_files\TimeDelay.mat');

bpFilt = designfilt('bandpassfir','FilterOrder',20,'CutoffFrequency1',1,'CutoffFrequency2',70,'SampleRate',1000);
taskdetails=cell(wue_no,3);
inx=[];
for i=1:wue_no
    leng=length(TimeDelay(i).Tdiff); %gives the number of sessions per patient 
    taskdetails(i,:)={wuelist{i},leng,notrials(i).sessno};
    for j=1:leng
        if any(TimeDelay(i).Tdiff)
            if i==4 && j==2
                [EEG, com] = pop_loadbv(['C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\EEGdata\' wuelist{i} '\EEG\'], ['grasp_' num2str(j+1) '.vhdr']);
            else
                [EEG, com] = pop_loadbv(['C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\EEGdata\' wuelist{i} '\EEG\'], ['grasp_' num2str(j) '.vhdr']);
            end
            % only truncation/addition of mean vector when time is negative
            if TimeDelay(i).Tdiff(j)>0
                AllEEGData_cut(i).sessno(j).eegData=EEG.data(:,floor(TimeDelay(i).Tdiff(j)):end);
                AllEEGData(i).sessno(j).eegData=[EEG.data];
            else
                inx=[inx; [i,j]];
                AllEEGData_cut(i).sessno(j).eegData=[EEG.data];
                AllEEGData(i).sessno(j).eegData=[EEG.data];
                %AllEEGData(i).sessno(j).eegData=[(repmat((mean(EEG.data,2)),1,-1*(TimeDelay(i).Tdiff(j)))) EEG.data];
            end
            % 1) filtering each 
            for chno=1:126
                AllEEGfiltered(i).sessno(j).eegData(chno,:)=filter(bpFilt,AllEEGData_cut(i).sessno(j).eegData(chno,:));
            end
        end 
    end
end 


% for that session where the eeg data starts later than kinematic, padding
%with channel vise mean values at the beginning: to maintain the mean when calculated  
for inxno=1:size(inx,1)
    AllEEGfiltered(inx(inxno,1)).sessno(inx(inxno,2)).eegData=[(repmat((mean(AllEEGfiltered(inx(inxno,1)).sessno(inx(inxno,2)).eegData,2)),1,floor(-1*(TimeDelay(inx(inxno,1)).Tdiff(inx(inxno,2))))))...
        AllEEGfiltered(inx(inxno,1)).sessno(inx(inxno,2)).eegData]; 
end
                                    

%% 2.1)reshaping filtered matrix to reflect each trial.
 %2.2) also creating a mega epoched 1 trial/session from onset of trial 1 to offset of end trial  
load('C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\MAT_files\onsetoffset.mat');
for i1=1:wue_no
    for j1=1:taskdetails{i1,2}
            for trialno=1:notrials(i1).sessno(j1).nos
                starttime1=(floor(1000*(onsetoffset(i1).sessno(j1).nos(1,trialno))));%1000* : conversion to ms
                starttime1_preonset=(floor(1000*(onsetoffset(i1).sessno(j1).nos(1,trialno))))-500;% -500 above to get 500 ms(500points) of pre-onset data
                endtime1=floor(1000*(onsetoffset(i1).sessno(j1).nos(2,trialno))); 
                AllEEGfiltered_reshaped(i1).sessno(j1).trial(trialno).eegdata=AllEEGfiltered(i1).sessno(j1).eegData(:,starttime1:endtime1);
                AllEEGfiltered_reshaped_prepoc(i1).sessno(j1).trial(trialno).eegdata=AllEEGfiltered(i1).sessno(j1).eegData(:,starttime1_preonset:endtime1);
            end
            starttime2=floor(1000*(onsetoffset(i1).sessno(j1).nos(1,1)))-500;%getting additional time before and
            endtime2=floor(1000*(onsetoffset(i1).sessno(j1).nos(2,end)))+500;%after the main epoc. Arbitrary
            AllEEGfiltered_reshaped_megaepoc(i1).sessno(j1).trial(trialno).eegdata=AllEEGfiltered(i1).sessno(j1).eegData(:,starttime2:endtime2);
            
    end
end
%% 3.1) removing bad channels: 
rej_count=0;rej_val=[];
for i1=1:wue_no
    for j1=1:taskdetails{i1,2}
        if i1==4 && j1==2
            [EEG, com] = pop_loadbv(['C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\EEGdata\' wuelist{i1} '\EEG\'], ['grasp_' num2str(j1+1) '.vhdr']);
            else
            [EEG, com] = pop_loadbv(['C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\EEGdata\' wuelist{i1} '\EEG\'], ['grasp_' num2str(j1) '.vhdr']);
        end
        for trialno=1:notrials(i1).sessno(j1).nos
            chnum=zeros(126,1);
            AllEEGfiltered_trial_mean(i1).sessno(j1).trial(trialno).eegdata=...
                mean(AllEEGfiltered_reshaped(i1).sessno(j1).trial(trialno).eegdata,1);
            avgtrialval=AllEEGfiltered_trial_mean(i1).sessno(j1).trial(trialno).eegdata;
            strdev=std(AllEEGfiltered_reshaped(i1).sessno(j1).trial(trialno).eegdata);
            for chno=1:126               
                temp=AllEEGfiltered_reshaped(i1).sessno(j1).trial(trialno).eegdata(chno,:);
                temp=bsxfun(@minus,temp,(avgtrialval));% from the global mean of all channels for that trial
                temp=bsxfun(@rdivide,temp,(strdev));
                percentage=(((sum(abs(temp)>4))/length(temp))*100);% 4x s.d.

                if percentage>2.5% 2.5% or greater : value based on Roy et.al.(2018)
                    chnum(chno,1)=1;
                    rej_count=rej_count+1;
                    rej_val=[rej_val;[i1 j1 trialno chno]];
                end
            end
            Channelrej(i1).sessno(j1).trlno(trialno).chrej=chnum;
            
%% 3.2) running channel interpolation: 
            % NOTE: the wue/sessions that have channels to be interpolated are not those that have -ve timediff values!           
            if ~isempty(rej_val)
            badchans=rej_val(ismember(rej_val(:,1:3),[i1 j1 trialno],'rows'),4);
            if ~isempty(badchans)              
                EEGTEMP=EEG;
                EEGTEMP.data=AllEEGfiltered_reshaped(i1).sessno(j1).trial(trialno).eegdata;
                EEGTEMP.pnts=size(AllEEGfiltered_reshaped(i1).sessno(j1).trial(trialno).eegdata,2);
                EEGOUT = eeg_interp(EEGTEMP, badchans, 'spherical');
                AllEEGfiltered_reshaped(i1).sessno(j1).trial(trialno).eegdata=EEGOUT.data;
                
                % same thing for prepoc
                EEGTEMP.data=AllEEGfiltered_reshaped_prepoc(i1).sessno(j1).trial(trialno).eegdata;
                EEGTEMP.pnts=size(AllEEGfiltered_reshaped_prepoc(i1).sessno(j1).trial(trialno).eegdata,2);
                EEGOUT = eeg_interp(EEGTEMP, badchans, 'spherical');
                AllEEGfiltered_reshaped_prepoc(i1).sessno(j1).trial(trialno).eegdata=EEGOUT.data;
            end
            end
%% 3.3) running channel average referencing: 
            % for -ve tdiff: it makes sense to have the padded values as
            % mean values of each channel. This will not create bias when
            % channel referencing if the values were set to zero. (these
            % channels would be then contributing less to the all channel avg.
            %reference
            cur_trl=AllEEGfiltered_reshaped(i1).sessno(j1).trial(trialno).eegdata;
            cur_trl=bsxfun(@minus, cur_trl, sum(cur_trl,1)/(126+1)); %divide by nchan+1: check Makoto's preprocessing pipeline
            AllEEGfiltered_reshaped_avgref(i1).sessno(j1).trial(trialno).eegdata=cur_trl;
            
            %prepoc implementation: 
            cur_trl=AllEEGfiltered_reshaped_prepoc(i1).sessno(j1).trial(trialno).eegdata;
            cur_trl=bsxfun(@minus, cur_trl, sum(cur_trl,1)/(126+1)); %divide by nchan+1: check Makoto's preprocessing pipeline
            AllEEGfiltered_reshaped_prepoc_avgref(i1).sessno(j1).trial(trialno).eegdata=cur_trl;
            
%% 3.4) cleanline 
            %{
            EEGclin=EEG;
            EEGclin.data=AllEEGfiltered_reshaped_avgref(i1).sessno(j1).trial(trialno).eegdata;
            EEGclin.pnts=size(AllEEGfiltered_reshaped_avgref(i1).sessno(j1).trial(trialno).eegdata,2);            
            EEG = pop_cleanline(EEGclin, 'Bandwidth',2,'ChanCompIndices',[1:2] ,...
                'SignalType','Channels','ComputeSpectralPower',false,'LineFrequencies',[50] ,...
                'NormalizeSpectrum',false,'LineAlpha',0.01,'PaddingFactor',2,'PlotFigures',false,...
                'ScanForLines',true,'SmoothingFactor',100,'VerbosityLevel',1,'SlidingWinLength',...
                EEGclin.pnts/EEGclin.srate,'SlidingWinStep',EEGclin.pnts/EEGclin.srate);
            
            % on prepoc
            EEGclin_prepoc=EEGTEMP;
            EEGclin_prepoc.data=AllEEGfiltered_reshaped_prepoc_avgref(i1).sessno(j1).trial(trialno).eegdata;
            EEGclin_prepoc.pnts=size(AllEEGfiltered_reshaped_prepoc_avgref(i1).sessno(j1).trial(trialno).eegdata,2);  
            EEG = pop_cleanline(EEGclin_prepoc, 'Bandwidth',2,'ChanCompIndices',[1:EEGclin_prepoc.nbchan] ,...
                'SignalType','Channels','ComputeSpectralPower',false,'LineFrequencies',[50] ,...
                'NormalizeSpectrum',false,'LineAlpha',0.01,'PaddingFactor',2,'PlotFigures',false,...
                'ScanForLines',true,'SmoothingFactor',100,'VerbosityLevel',1,'SlidingWinLength',...
                EEGclin_prepoc.pnts/EEGclin_prepoc.srate,'SlidingWinStep',EEGclin_prepoc.pnts/EEGclin_prepoc.srate);
           %}
        end
    end
end

%% 4) runing ICA:
% PCA is turned on since the data points are not enough for reliable ica. 
% doing this on epoc concatenated vector

            
            
            
            