% the final structure 'events_onoffsetkin_data.sessno.data' has time data in ms (same as indices since sampling is 1000hz)
%size is ( (no. epocs+1)=5) X (no of trials). Each is the start index of
%epoc from the ALLEEGdata_complete data. 
%% this extracts data of all KIN files in one structure: 

kinloc='C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\KIN events';
wuelist={'wue02','wue03','wue06','wue09','wue10','wue11'};
load('C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\MAT_files\diffsrtval_mod.mat');


for t=1:length(wuelist)
    nameslist=dir([kinloc filesep wuelist{t} filesep '*.mat']);
    allnames={(nameslist.name)};
    for l=1:length(allnames)
        kin_data(t).sessno(l).data=load([kinloc filesep wuelist{t} filesep nameslist(l).name]);
    end
end

 
%% extracting segments of sub epocs in every trial:
% 1) rest phase (pre epoc, baseline)
% 2) towards target phase 
% 3) at target phase 
% 4) return phase  
% ( the data times will be compatible to the cut version of session data,
% same as 'new_onsetimes' variable output from the pipeline
cntr=0;
for t=1:length(wuelist)  
    nameslist=dir([kinloc filesep wuelist{t} filesep '*.mat']);
    allnames={(nameslist.name)};
    for l=1:length(allnames)
        cntr=cntr+1;
        
        currstruct=kin_data(t).sessno(l).data;
        currtimedata=int32(1000.*[currstruct.Onset_event;currstruct.Movement_event;...
            currstruct.Return_event;currstruct.Offset_event]);
        currtimedata=[(currtimedata(1,1)-diffsrtval_mod(cntr)),(currtimedata(end,1:end-1)+1);...
            currtimedata];
        events_onoffsetkin_data(t).sessno(l).data=currtimedata-(currtimedata(1,1))+1;
%         tempdataa=[];subtracter=1000*kin_data(t).sessno(l).data.Onset_event(1,1);subtracter1=0;subtracter3=0;
%         
%         for r=1:length(kin_data(t).sessno(l).data.Trial_time)
%             
%             time2=1000*kin_data(t).sessno(l).data.Onset_event(1,r);
%             if r==1
%                 subtr=time2-diffsrtval_mod(cntr);
%             else 
%                 subtr=time2-time5;
%             end
%             time1=time2pre-subtr+1;
%             time3=1000*kin_data(t).sessno(l).data.Movement_event(1,r)-subtr+1;
%             time4=1000*kin_data(t).sessno(l).data.Return_event(1,r)-subtr+1;
%             time5=1000*kin_data(t).sessno(l).data.Offset_event(1,r)-subtr+1; 
%             time2pre=time2;
% %             if r>1
% %                 subtracter1=(time1-subtracter2-1);
% %             end
% %             subtracter3=subtracter3+subtracter1;
% %             time3=1000*kin_data(t).sessno(l).data.Movement_event(1,r);
% %             time4=1000*kin_data(t).sessno(l).data.Return_event(1,r);
% %             time5=1000*kin_data(t).sessno(l).data.Offset_event(1,r); 
% %             %time6=(1000*kin_data(t).sessno(l).data.Offset_event(1,r))+500;
%             tempdataa=[tempdataa,([time1;time2;time3;time4;time5])];
%            
%         end
%         events_onoffsetkin_data(t).sessno(l).data=int32(tempdataa);      
    end
end
clearvars nameslist allnames currstruct currtimedata cntr t l



