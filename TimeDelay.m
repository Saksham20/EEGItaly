
asdf=xlsread('C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\TimeOffsetData.xlsx');
for i=1:6
    TimeDelay(i).Tdiff=asdf(4+7*(i-1),:);
    TimeDelay(i).Tdiff(isnan(TimeDelay(i).Tdiff))=[];
end