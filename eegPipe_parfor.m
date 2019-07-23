
wue_no=6;
wuelist={'wue02','wue03','wue06','wue09','wue10','wue11'};
kinloc='C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\KIN events';
load('C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\MAT_files\notrials.mat');
load('C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\MAT_files\TimeDelay.mat');
load('C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\MAT_files\onsetoffset.mat');
load('C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\MAT_files\diffsrtval.mat');

%goodelectrodes=[11:36 47:50 53:58];
goodelectrodes=[1:126];
taskdetails=cell(wue_no,3);
freqcutoffout=100;freqcutoffin=1;

inx=[];rej_val=[];chnos=1:length(goodelectrodes);
currwue=1:6;

o.epoch_interp_options.rejection_options.measure = [1 1 1 1];
o.epoch_interp_options.rejection_options.z = [3 3 3 3];
            
for i=currwue
    leng=length(TimeDelay(i).Tdiff); %gives the number of sessions per patient 
    taskdetails(i,:)={wuelist{i},leng,notrials(i).};
    for j=1:leng  
            fprintf('starting operation for wue %d session %d \n\n\n',i,j);            
            %% 1) Data extraction------------ 
            fprintf('data extraction\n');tic;
            if i==4 && j==2
                [EEG, ~] = pop_loadbv(['C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\EEGdata\' wuelist{i} '\EEG\'], ['grasp_' num2str(j+1) '.vhdr']);
            else
                [EEG, ~] = pop_loadbv(['C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\EEGdata\' wuelist{i} '\EEG\'], ['grasp_' num2str(j) '.vhdr']);
            end         
            if TimeDelay(i).Tdiff(j)>0
                AllEEGData_cut(i).sessno(j).data=EEG.data(goodelectrodes,floor(TimeDelay(i).Tdiff(j)):end);
                AllEEGData(i).sessno(j).data=[EEG.data(goodelectrodes,:)];
            else
                inx=[inx; [i,j]];             
                AllEEGData(i).sessno(j).data=[EEG.data(goodelectrodes,:)];
                AllEEGData_cut(i).sessno(j).data=[zeros(length(goodelectrodes),int32(-1*(TimeDelay(i).Tdiff(j)))) EEG.data(goodelectrodes,:)];
            end
            % clear EEG structures to save RAM: 
            clear ('vars','EEG');
            %starting the eegdata from the onset time of the first trial (can be padded with zeros later for compatibility): 
            diffsrtval_mod=diffsrtval;
            diffsrtval_mod(:,3)=diffsrtval(:,3)-50;% 50ms is setting a relief from the TENS stop and start of pre-onset time
            diffsrtval_mod((diffsrtval(:,3)>2000),3)=2000;%setting upper limit for the starttime offset. dont want it to be very large since the mean of all the inter trial times is 1750ms           
            tstart=int32((onsetoffset(i).sessno(j).nos(1,1))*1000)...
                -diffsrtval_mod(ismember(diffsrtval_mod(:,1:2),[i,j],'rows'),3);% diffsrtval is a vector of the differences between the end of the tens spikes and start of trial times.  
            tend=int32((onsetoffset(i).sessno(j).nos(2,end))*1000); % no need for end offset 
            AllEEGData_cut_practical(i).sessno(j).data=AllEEGData_cut(i).sessno(j).data(:,tstart:tend);
            fprintf('done in %d s\n\n',toc);
            
            %pre assignments: 
%             AllEEGData_cut_practical_filtered(i).sessno(j).data=0;
%             AllEEGData_cut_practical_filtered_chrej(i).sessno(j).data=0;
%             AllEEGData_wICAcleaned(i).sessno(j).data=0;
%             AllEEGData_SASICAcleaned(i).sessno(j).data=0;
%             AllEEGData_SASICAcleaned_reshaped_trlinterpd(i).sessno(j).data=0;
%             AllEEGData_complete(i).sessno(j).data=0;
%             
%             for trlnum=1:size(onsetoffset(i).sessno(j).nos,2)
%                 AllEEGData_SASICAcleaned_reshaped(i).sessno(j).trlno(trlnum).data=0;
%                 AllEEGData_complete_reshaped(i).sessno(j).trlno(trlnum).data=0;
%                 
%             end
    end
end

%% variable declarations: 
eegnames={'AllEEGData_cut_practical_filtered','AllEEGData_cut_practical_filtered_chrej'...
                ,'AllEEGData_wICAcleaned','AllEEGData_SASICAcleaned','AllEEGData_complete'};
reshaped_eegnames={'AllEEGData_SASICAcleaned_reshaped'...
    ,'AllEEGData_SASICAcleaned_reshaped_trlinterpd','AllEEGData_complete_reshaped'};

othernames={'wicaArt','var_afterWICA','ICs_keptsasica','var_afterSASICA','chrej',...
    'new_onsetimes'};
icadata={'icaact','icawinv','icasphere','icaweights','icachansind'};
allvars=[eegnames reshaped_eegnames othernames icadata];


%% 
for i=currwue
    leng=length(TimeDelay(i).Tdiff); %gives the number of sessions per patient 
    taskdetails(i,:)={wuelist{i},leng,notrials(i).sessno};
    sesscel=cell(leng,(numel(allvars)));
        parfor j=1:leng  
            fprintf('starting operation for wue %d session %d \n\n\n',i,j); 
            %% 2) filtering----------------- 1-100HZ
            fprintf('filtering\n');tic;
            eegtemp=eegdata_struct_gen(AllEEGData_cut_practical(i).sessno(j).data,goodelectrodes);
            [eegoutfilt, com, b] = pop_eegfiltnew(eegtemp,freqcutoffin ,freqcutoffout);
            AllEEGData_cut_practical_filtered=eegoutfilt.data;
            fprintf('done in %d s\n\n',toc);
            
            %% 3) Cleanline------------------
            fprintf('CleanLine OP\n');tic;
            eegoutcleanline = pop_cleanline(eegoutfilt, 'Bandwidth',2,'ChanCompIndices',[1:eegoutfilt.nbchan] ,...
                'SignalType','Channels','ComputeSpectralPower',false,'LineFrequencies',[50] ,...
                'NormalizeSpectrum',false,'LineAlpha',0.01,'PaddingFactor',2,'PlotFigures',false,...
                'ScanForLines',true,'SmoothingFactor',100,'VerbosityLevel',1,'SlidingWinLength',...
                eegoutfilt.pnts/eegoutfilt.srate,'SlidingWinStep',eegoutfilt.pnts/eegoutfilt.srate);
            
            fprintf('done in %d s\n\n',toc);
            
            %% 4) Bad channel rejection---------------
            %4.1) HAAPE inplementation: 
            fprintf('Bad channel rejection\n');tic;
            no_iters=2; 
            eegforiter=eegoutcleanline;
            
            chlist_orig={eegforiter.chanlocs.labels};
            for iter=1:no_iters              
                [EEG_temp,idx] = pop_rejchan(eegforiter, 'elec',[1:eegforiter.nbchan],'threshold',[-3 3],'norm','on','measure','spec','freqrange',[1 70]);
                eegforiter=EEG_temp;
                fprintf('.');
            end
            chlist_new={eegforiter.chanlocs.labels};
            chan_acept=chnos(ismember(chlist_orig,chlist_new));
            chan_rej=chnos(~ismember(chlist_orig,chlist_new));
            rej_val=chan_rej(:);
            
            AllEEGData_cut_practical_filtered_chrej=eegforiter.data;
            fprintf('done in %d s\n\n',toc);
            % 4.2) PREP implementation: 
            
            
            
            %% 6) w-ICA-------------------
            fprintf('wICA operation:\n');tic;
            [wIC, A, W, IC] = wICA(eegforiter,'runica', 1, 0,1000,5);% ##vary the threshold multiplier(3rd arg)
            wicaArt=A*wIC;
            AllEEGData_wICAcleaned=... %clean data gen
                eegforiter.data-wicaArt;
            % 6.2) % of variance kept after wica: 
            [~,vaf1]=compvar(eegforiter.data,wIC,A,[1:eegforiter.nbchan]);
            var_afterWICA=vaf1;
            fprintf('done in %d s\n\n',toc);
            
            %% 7) ICA---------------- % perform PCA before ICA to improve ? 
            fprintf('normal ICA:\n');tic;
            eegforICA=eegforiter;
            
            
            % change here for wica off or on
            eegforICA.data=AllEEGData_wICAcleaned
            
            eegforICA = pop_runica(eegforICA, 'extended',1,'interupt','on');
            icaact=eegforICA.icaact; % useful to store to plot the ICs later for manual artifact selection
            icawinv=eegforICA.icawinv;
            icasphere=eegforICA.icasphere;
            icaweights=eegforICA.icaweights;
            icachansind=eegforICA.icachansind;
            fprintf('done in %d s\n\n',toc);
            

            %% 8.2) SASICA------------------------ 
            fprintf('performing SASICA\n');tic;
            [eegSasica, ~] = eeg_SASICA(eegforICA,1,1,0,0,1,0,0,0,0,1);% check 'def = SASICA('getdefs');' to know the seq of 0,1
            
            %save figures
            figHandles = get(groot, 'Children');
            fign={figHandles(:).Name};
            savfig=regexp(fign,'Reject');
            savfig=find(cellfun(@(x) ~isempty(x), savfig));
            savloca='C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\scripts\sasicaout\';
            currfol=[savloca sprintf('wue%d',i)];
            mkdir(currfol);
            for ir=1:length(savfig)
                figure(ir);
                print([currfol filesep sprintf('sessnono_%d_%d',j,ir)],'-djpeg');
            end
            close all;
            
            %ICs accepted info:
            ICs_keptsasica=eegSasica.reject.SASICA; 
            
            %signal cleaning:
            fldnames=fieldnames(eegSasica.reject.SASICA);templog=zeros(1,size(eegSasica.icaact,1));
            for fld=1:numel(fldnames)
                if islogical(eegSasica.reject.SASICA.(fldnames{fld}));
                   templog=[templog | eegSasica.reject.SASICA.(fldnames{fld})];
                end
            end
            tempids=~(logical(templog));
            
            AllEEGData_SASICAcleaned= eegSasica.icawinv(:,tempids)*eegSasica.icaact(tempids,:);
            ICA_actsas=eegSasica.icaact;ICA_winvsas=eegSasica.icawinv;
            
            %variance kept 
            [~, var_afterSASICA] =compvar(eegSasica.data, ICA_actsas, ICA_winvsas, find(tempids));
            fprintf('done in %d s\n\n',toc);
            
            
            %% 9) Reshaping for trials and then pad for consistency with TENS/KINdata------------------ 
            fprintf('reshaping operation\n');tic;
            temponoffsettimes=(1000.*(onsetoffset(i).sessno(j).nos))-...
                ((1000.*onsetoffset(i).sessno(j).nos(1,1)))+...
                (diffsrtval_mod(ismember(diffsrtval_mod(:,1:2),[i,j],'rows'),3))+1; %since the offset is defined by diffstrval
            temponoffsettimes=int32(temponoffsettimes);
            AllEEGData_SASICAcleaned_reshaped=struct();
            for trlnum=1:size(temponoffsettimes,2)                      
                if trlnum==1
                    startm=1;           
                else startm=(temponoffsettimes(2,trlnum-1))+1;
                end
                endtm=(temponoffsettimes(2,trlnum));
                AllEEGData_SASICAcleaned_reshaped(trlnum).data=...
                    AllEEGData_SASICAcleaned(:,startm:endtm);
                fprintf('.');
            end
            fprintf('done in %d s\n',toc);
            
            %% 10)Channel rejection for each individual trial: 
            % following the FASTER criteria for rejecting channels if they
            % fail even in one of the z values. 
            fprintf('Channel rejection for individual trial operation\n');tic;
            
            eeg_chans=[1:size(AllEEGData_SASICAcleaned_reshaped,1)];
            numberoftrials=size(temponoffsettimes,2);
            eegdata_final_concatd=[];new_onoffset=1;szcount=1;
            chrej=struct();
            for trlnum=1:numberoftrials %looping once again over all the trials
                
                list_properties = listPropGen(AllEEGData_SASICAcleaned_reshaped,trlnum,eeg_chans);
                chrej(trlnum).data=eeg_chans(logical(min_z(list_properties,o.epoch_interp_options.rejection_options)));
                
                eegforinterp=eegdata_struct_gen_new(AllEEGData_SASICAcleaned_reshaped(trlnum).data,eegSasica);
                interpeeg = eeg_interp(eegforinterp, chrej(trlnum).data, 'spherical');
                
                AllEEGData_SASICAcleaned_reshaped_trlinterpd=interpeeg.data;
                eegdata_final_concatd=[eegdata_final_concatd, AllEEGData_SASICAcleaned_reshaped_trlinterpd];
                szcount=szcount+size(AllEEGData_SASICAcleaned_reshaped_trlinterpd,2);
                new_onoffset=[new_onoffset, szcount];% the onoffset times will change since its just appending the snipped data
            end
            new_onsetimes=new_onoffset;
            fprintf('done in %d s\n',toc);
            
            %% 11) Interpolation of the mega epoc rejected channels: 
            fprintf('Interpolation of the mega epoc rejected channels\n');tic;
            
            %rebuilding the main eegdata structure with new conct data
            ch_rejnos=rej_val(ismember((rej_val(:,1:2)),[i,j],'rows'),3); % since the data will not be graphically continuous...           
            eegformegaInterp=eegdata_struct_gen(eegdata_final_concatd,goodelectrodes);% will interp use any neighbouring values also for every point?            
            rebuild=zeros(length(goodelectrodes),size(eegdata_final_concatd,2));
            rebuild(setdiff([1:length(goodelectrodes)],ch_rejnos),:)=eegdata_final_concatd;
            eegformegaInterp.data=rebuild;
            
            interpeeg = eeg_interp(eegformegaInterp, ch_rejnos, 'spherical');
            fprintf('done in %d s\n',toc);

            %% 12) Referencing: 
             % check Happe
            fprintf('Final referencing of concatened channels and reshaping again\n');tic;
            refrd_eeg_final = pop_reref(interpeeg, []);
            AllEEGData_complete=refrd_eeg_final.data;% alleegdatacomplete is made by appending the reshaped mara. so middle sections are missing. 
            % final reshaping
            AllEEGData_complete_reshaped=struct();
            for trlnum=1:(size(new_onoffset,2)-1)
                AllEEGData_complete_reshaped.trial(trlnum)=...
                    AllEEGData_complete(:,new_onoffset(trlnum):(new_onoffset(trlnum+1)-1));
            end
            fprintf('done in %d s\n',toc);
            
            %% 13) Reassignment into sliced variable:            
                               
            for g=1:numel(sesscel)
                eval(['sesscel{j,' sprintf('%d',g) '=allvars{1,' sprintf('%d',g) '};']);
            end
        end
        EEGALL(i).data=sesscel;
end 
MATinfo='select channels with sasica with wica';
save('C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\MAT_files\PipeV2allout_v4_sasica_selectchans_WICAon');  


