
% frequency plotting for each of the wues: 
% should run spect_analysisV2.m before. 
imgstoreloc='C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\scripts\sasicaout\sasica\selectchans\wicaOn\pipeout\norm\';
mkdir(imgstoreloc);
%goodelectrodes=[11:36 53:58 76:77 122:123 125:126];

spectog_wue_current=spectog_wue_filtered_sasica;

for t=currwue
    %% temperature plot for all chans-------
    freqlab=(spectog_wue_current(t).freqs<100);
    imge1=imagesc(spectog_wue_current(t).smavgdata(:,freqlab));
    
    colorbar; xlabel(['Hz']);ylabel(['channels']);
    curax=gca;
    curax.XTickLabelRotation=90;
    curax.XTick=linspace(0,sum(freqlab),20);
    xticklabels = round(linspace(0,70,length(curax.XTick)));
    set(gca, 'XTickLabel', xticklabels);
    
    print([imgstoreloc sprintf('wue_%d_all',t)],'-djpeg');close all;
    favfreqs=find(spectog_wue_current(t).freqs<70);
    
    %% all chans graph plot-------
    nbchans=length(goodelectrodes);
    count=0;col=hsv(nbchans);
    figure();
    for allchan=1:nbchans;
            count=count+1;
        plot(spectog_wue_current(t).freqs(favfreqs),(spectog_wue_current(t).smavgdata(allchan,favfreqs)),...
            'Color',col(count,:),'LineStyle','--');
        hold on;
    end
    plot(spectog_wue_current(t).freqs(favfreqs),nanmean(spectog_wue_current(t).smavgdata(:,favfreqs)),...
        'Color',[0 0 0],'LineStyle','-','LineWidth',2);grid on;
    xlabel(['Hz']);ylabel(['PSD avg across sessions']);
    print([imgstoreloc sprintf('wue_%d_avgallchans',t)],'-djpeg');close all;
    
%     %% SMC-------
%     count=0;SMC=[11:36 47:50 53:58];col=hsv(length(SMC));
%     leg={};
%     for favchans=SMC
%         count=count+1;
%         plot(spectog_wue_current(t).freqs(favfreqs),(spectog_wue_current(t).smavgdata(favchans,favfreqs)),...
%             'Color',col(count,:),'LineStyle','--');
%         hold on;leg=[leg, sprintf('chno. %d',favchans)];
%     end
%     plot(spectog_wue_current(t).freqs(favfreqs),mean(spectog_wue_current(t).smavgdata(SMC,favfreqs)),...
%             'Color',[0 0 0],'LineStyle','-','LineWidth',2);legend(leg);grid on;
%         xlabel(['Hz']);ylabel(['PSD avg across sessions']);
%     print([imgstoreloc sprintf('wue_%d_avgSMCchans',t)],'-djpeg','-r600');close all;
end

%% NOrmalizing the spectdata subepoc for every trial:
% run spect_analysisv2.m before
spectdata_SUBEPOC_select=spectdata_SUBEPOC;
spectdata_SUBEPOC1=spectdata_SUBEPOC_select;
spectdata_SUBEPOC_select_norm=spectdata_SUBEPOC1;
for wueno=currwue
    sz=size(spectdata_SUBEPOC1(wueno).data);
    loopsz=prod(sz(3:5));
    freqlen=size(spectdata_SUBEPOC1(wueno).data,2);
    freqlen=[1:freqlen].*500/freqlen;
    goodidx=find(freqlen<freqcutoffout);
    for i=1:sz(3)
        for j=1:sz(4)
            for k=1:sz(5)
                mn=nanmean(spectdata_SUBEPOC1(wueno).data(:,goodidx,i,j,k),2);
                stde=nanstd(spectdata_SUBEPOC1(wueno).data(:,goodidx,i,j,k),0,2);
                mn=repmat(mn,1,length(goodidx));
                stde=repmat(stde,1,length(goodidx));
                spectdata_SUBEPOC_select_norm(wueno).data(:,goodidx,i,j,k)=...
                    (spectdata_SUBEPOC1(wueno).data(:,goodidx,i,j,k)-mn)./stde;
            end
        end
    end  
end
%% Finding sig difference between 4 conditions for each frequencies across wues. 
%                 | All 4 different | Pre/post & Movmt 
% ----------------------------------------------------
% 1)Indi freq     |        x        |      _/
% 2)Band specific |        x        |      _/ 

% 1)indi freq ----------
spectdata_SUBEPOC1=spectdata_SUBEPOC_select_dB;
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
            dat1=dat1(:);
            dat2=reshape(NData(chans,freq,2,:,:),size(NData,4),size(NData,5));
            dat2=dat2(:);          
            pvaltemp(chans,freq)=anova1([dat1,dat2],[1 2],'off');
            sigValtemp(chans,freq)=(pvaltemp(chans,freq)<0.05);
        end
    end
    pval(wueno).data=pvaltemp;
    sigVal(wueno).data=sigValtemp;
end

%% 2) band specific ------------
alpha=[8 12];beta=[13 30];
bands=[alpha; beta];
spectdata_SUBEPOC1=spectdata_SUBEPOC_meansub;
for wueno=currwue
    %--- taking the mean of spectrograms of the 3 movement epocs
    
    sz=size(spectdata_SUBEPOC1(wueno).data);
    sz(1,3)=sz(1,3)-2;% since one is rest and other is movement 
    NData=zeros(sz);
    NData(:,:,2,:,:)=nanmean(spectdata_SUBEPOC1(wueno).data(:,:,[2 3 4],:,:),3);% avg of last 3 epocs of movement
    NData(:,:,1,:,:)=(spectdata_SUBEPOC1(wueno).data(:,:,1,:,:));
    pvalbandtemp=zeros(sz(1),size(bands,2));pvaltemp_anova=pvalbandtemp;
    sigValbandtemp=zeros(size(pvalbandtemp));
    sigValtemp_anova=sigValbandtemp;
    moreORless=zeros(size(pvalbandtemp));
    freqs=size(spectdata_SUBEPOC1(wueno).data,2);
    freqs=[1:freqs].*500/freqs;
    bandid=arrayfun(@(x) dsearchn(freqs',x),bands);
    permu_no_iters=1000; 
    surrogate_pxl_map=zeros(sz(1),size(bands,2),permu_no_iters);
    for chans=1:sz(1)
        for freq=1:size(bands,2)
            dat1=nanmean(NData(chans,bandid(freq,1):bandid(freq,2),1,:,:),2);%rest stage
            dat1=reshape(dat1,size(NData,4),size(NData,5));
            dat1=dat1(:);dat1(isnan(dat1))=[]; 
            dat2=nanmean(NData(chans,bandid(freq,1):bandid(freq,2),2,:,:),2);%moving stage
            dat2=reshape(dat2,size(NData,4),size(NData,5));
            dat2=dat2(:);dat2(isnan(dat2))=[];  
            %[~,pvalbandtemp(chans,freq)]=ttest(dat1,dat2);
            rawvalues(chans,freq)=mean(dat1-dat2);
            [pvalbandtemp(chans,freq),surrogate_pxl_map(chans,freq,:)]=permutation_test(dat1,dat2,permu_no_iters);
            
            pvaltemp_anova(chans,freq)=anova1([dat1,dat2],[1 2],'off');
            sigValbandtemp(chans,freq)=(pvalbandtemp(chans,freq)<0.05);
            sigValtemp_anova(chans,freq)=(pvaltemp_anova(chans,freq)<0.05);
            
            moreORless(chans,freq)=(nanmean(dat1)>nanmean(dat2));% baseline should have more power when significant 
        end
    end
    
    sigValtemp_pixel=zeros(size(rawvalues));sigValtemp_cluster=sigValtemp_pixel;
    for freq=1:size(bands,2) %looping for every freq band since Im assuming alpha, beta, gamma are not related and each are modulated independantly by the task. this will change in case of all freq looping 
        %pixel wise threshold---------------:     
        max_pxls=zeros(permu_no_iters,1);
        for h=1:permu_no_iters
            max_pxls(h,1)=max(max(surrogate_pxl_map(:,freq,h)));
        end
        uprthresh_rawvalues=prctile(max_pxls,95);
        sigValtemp_pixel((rawvalues(:,freq)>uprthresh_rawvalues),freq)=1;


        %cluster wise thresholds--------------: 
        max_clusters=zeros(permu_no_iters,1);
        for h=1:permu_no_iters
            temp=surrogate_pxl_map(:,freq,h);
           
            zscored=zeros(size(temp));zscored_sig=zscored;
            zscored=(temp-mean(temp(:)))/(std(temp(:)));
            zscored_sig(zscored>norminv(1-0.05))=1;
            CC = bwconncomp(zscored_sig);
            temp1=temp(:);
            sum_teststat = cellfun(@(x) sum(temp1(x)),CC.PixelIdxList);
            if ~isempty(sum_teststat)
                max_clusters(h,1) = max(sum_teststat);
            else 
                max_clusters(h,1) = 0;
            end
        end
        uprthresh_rawvalues_cluster=prctile(max_clusters,95);

        temp=rawvalues(:,freq);
        zscored=zeros(size(temp));zscored_sig=zscored;
        zscored=(temp-mean(temp(:)))/(std(temp(:)));
        zscored_sig(zscored>norminv(1-0.05))=1;
        CC = bwconncomp(zscored_sig);
        temp1=temp(:);
        sum_teststat = cellfun(@(x) sum(temp1(x)),CC.PixelIdxList);
        goodclusters=zeros(size(sum_teststat));
        goodclusters(sum_teststat>uprthresh_rawvalues_cluster)=1;
        goodcluster_pixels=CC.PixelIdxList(logical(goodclusters));%cell containing all good pixels 
        sigValtemp_cluster1=zeros(size(temp));pxlsout=[];
        for te=1:numel(goodcluster_pixels)
            temp=goodcluster_pixels{te};
            pxlsout=[pxlsout;temp(:)];
        end
        sigValtemp_cluster1=sigValtemp_cluster1(:);
        sigValtemp_cluster1(pxlsout)=1;
        sigValtemp_cluster1=reshape(sigValtemp_cluster1,size(rawvalues(:,freq)));
        sigValtemp_cluster(:,freq)=sigValtemp_cluster1;
    end
    clearvars uprthresh_rawvalues sigValtemp_cluster1
    
    %result collection in structures----------
    pval_band_norm_anova(wueno).data=[pvaltemp_anova goodelectrodes' moreORless];
    sigVal_band_norm_anova(wueno).data=[sigValtemp_anova goodelectrodes' moreORless];
    pval_band_norm(wueno).data=[pvalbandtemp goodelectrodes' moreORless];
    sigVal_band_norm(wueno).data=[sigValbandtemp goodelectrodes' moreORless];  
    sigVal_band_pixel(wueno).data=[sigValtemp_pixel goodelectrodes' moreORless];  
    sigVal_band_cluster(wueno).data=[sigValtemp_cluster goodelectrodes' moreORless];  
end

% implement bonferroni correction to get new sigval
pvalbon=0.03;
sigVal_band_norm1_bon=sigVal_band_norm;
for wueno=currwue
    temp1=pval_band_norm(wueno).data(:,1:2);
    temp2=zeros(size(sigVal_band_norm(wueno).data(:,1:2)));
    temp2(temp1<pvalbon)=1;
    sigVal_band_norm1_bon(wueno).data(:,1:2)=temp2;
end
fprintf('done\n');
clearvars NData pvalbandtemp sigValbandtemp freqs dat1 dat2 pvalbandtemp sigValbandtemp moreORless bandid temp
%%
% wues with sig alpha and beta variation:
% only those where baseline power decreses during task
sigVal_band1=sigVal_band_norm;
sigVal_band2=sigVal_band_norm_anova;
sigVal_band3=sigVal_band_norm1_bon;%this is after bonferroni correction
sigVal_band4=sigVal_band_pixel;
sigVal_band5=sigVal_band_cluster;

for wueno=currwue
    idx1=(sigVal_band1(wueno).data(:,1)).*double(sigVal_band1(wueno).data(:,4));
    wue(wueno).alphachans=sigVal_band1(wueno).data(logical(idx1),3);
    idx2=(sigVal_band1(wueno).data(:,2)).*double(sigVal_band1(wueno).data(:,5));
    wue(wueno).betachans=sigVal_band1(wueno).data(logical(idx2),3);
    
    %with bonferroni
    idx11=(sigVal_band3(wueno).data(:,1)).*double(sigVal_band3(wueno).data(:,4));
    wue_bon(wueno).alphachans=sigVal_band3(wueno).data(logical(idx11),3);
    idx12=(sigVal_band3(wueno).data(:,2)).*double(sigVal_band3(wueno).data(:,5));
    wue_bon(wueno).betachans=sigVal_band3(wueno).data(logical(idx12),3);
    %with pixel wise thresholding
    idx41=(sigVal_band4(wueno).data(:,1)).*double(sigVal_band4(wueno).data(:,4));
    wue_pixel(wueno).alphachans=sigVal_band4(wueno).data(logical(idx41),3);
    idx42=(sigVal_band4(wueno).data(:,2)).*double(sigVal_band4(wueno).data(:,5));
    wue_pixel(wueno).betachans=sigVal_band4(wueno).data(logical(idx42),3);
    %with cluster wise thresholding
    idx51=(sigVal_band5(wueno).data(:,1)).*double(sigVal_band5(wueno).data(:,4));
    wue_cluster(wueno).alphachans=sigVal_band5(wueno).data(logical(idx51),3);
    idx52=(sigVal_band5(wueno).data(:,2)).*double(sigVal_band5(wueno).data(:,5));
    wue_cluster(wueno).betachans=sigVal_band5(wueno).data(logical(idx52),3);
 
    idx1_anova=(sigVal_band2(wueno).data(:,1)).*double(sigVal_band2(wueno).data(:,4));
    wue_anova(wueno).alphachans=sigVal_band2(wueno).data(logical(idx1_anova),3);
    idx2_anova=(sigVal_band2(wueno).data(:,2)).*double(sigVal_band2(wueno).data(:,5));
    wue_anova(wueno).betachans=sigVal_band2(wueno).data(logical(idx2_anova),3);
end

fprintf('done\n');


%% -------plotting freq responses for each of the sub epocs within trial per
%wue (avg across trials and sessions) 
% mkdir('C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\scripts\spectograms\EPOCs\wICAondB\');
imgstoreloc='C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\scripts\sasicaout\sasica\wicaOn\pipeout\norm\';

spectdata_SUBEPOC1=spectdata_SUBEPOC_norm;
for wueno=2
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
%         ax=gca;yvl=ax.YTick;yvl=yvl(2);
%         sigdat=sigVal(wueno).data(chnonew,:);
%         plot(freqlen(goodidx),yvl*sigdat(goodidx),...
%                 '-*');
        title(sprintf('wue no. %d, chno. %d (normalized, avg across trial/sess)',wueno,goodelectrodes(chnonew)));
        legend({'pre','reach','grasp','reverse'});grid on
        xlabel(['Hz']);ylabel(['PSD avg across sessions']);
        print([imgstoreloc sprintf('wueno%d',wueno) filesep...
            sprintf('chno_%d ',goodelectrodes(chnonew))],'-djpeg');close all;
    end
end



% %% plot variations within each sub epoc
% wueno=3;chnonew=17;
% for t=2:4; 
%     col1=hsv(t(end));
%     for h=1:(size(spectdata_SUBEPOC(wueno).data,5))  
%         
%         for g=1:(size(spectdata_SUBEPOC(wueno).data,4))              
%             plot(freqlen(goodidx),((spectdata_SUBEPOC(wueno).data(chnonew,(goodidx),t,g,h))),'--','Color',col1(t,:)); hold on; 
%         end
%     end
%     plot(freqlen(goodidx),mean(mean(spectdata_SUBEPOC(wueno).data(chnonew,(goodidx),t,:,:),4),5),'Color',col1(t,:),'LineWidth',2); hold on; 
%     legend({'pre','reachon','hold','reachoff'})
% end
% %plot(freqlen(goodidx),mean(mean(spectdata_SUBEPOC(wueno).data(chnonew,(goodidx),t,:,:),4),5),'-k','LineWidth',2); hold on; 



%% Plotting the voltage data for all channels 
imgstoreloc1='C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\scripts\sasicaout\sasica\selectchans\wicaOn\raw\';
mkdir(imgstoreloc1);
names={'AllEEGData_cut_practical_filtered','AllEEGData_wICAcleaned','AllEEGData_SASICAcleaned','AllEEGData_complete'};

for nimes=1:length(names)
for wueno=1%currwue
    eval(['dataNow=' names{1,nimes}]);
    nosess=length(dataNow(wueno).sessno);
    mkdir([imgstoreloc1 filesep sprintf('wue%d',wueno)]);
    for ns=3%1:nosess 
        nochans=size(dataNow(wueno).sessno(ns).data,1);
        mx=max(max(dataNow(wueno).sessno(ns).data,[],2));
        mn=min(min(dataNow(wueno).sessno(ns).data,[],2));
        relief=-0.4*(mx-mn);
        col=hsv(nochans/6);
        col=repmat(col,6,1);
        figure();
        for nc=1:nochans     
            plot(dataNow(wueno).sessno(ns).data(nc,:)+(nc-1)*(mx-mn+relief),'Color',col(nc,:));hold on
        end
        title({sprintf('wue %d, sessno %d',wueno,ns),regexptranslate('escape',names{1,nimes})},'Interpreter','none');
        axs=gca;
        ylims=axs.YLim;
        lp=size(events_onoffsetkin_data(wueno).sessno(ns).data,2);
        for lpp=1:lp
            pt=events_onoffsetkin_data(wueno).sessno(ns).data(1,lpp);    
            line([pt pt],ylims,'Color',[0 0 0],'LineStyle',':');
        end
        %print([imgstoreloc1 filesep sprintf('wue%d',wueno) filesep sprintf('wue_%d_sessno_%d_%d',wueno,ns,nimes)],'-djpeg');close all;
        
    end
end
end

