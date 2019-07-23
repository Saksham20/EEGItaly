wue_no=6;
wuelist={'wue02','wue03','wue06','wue09','wue10','wue11'};
kinloc='C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\KIN events\';
for k=1:wue_no
    leng=length(TimeDelay(k).Tdiff);
    for j=1:leng
        temp=load([kinloc wuelist{k} filesep 'KIN_events_grasp_' num2str(j) '.mat']);
        notrials(k).sessno(j).nos=length(temp.Trial_time);
        onsetoffset(k).sessno(j).nos=[temp.Onset_event;temp.Offset_event].*1000;
    end
end
