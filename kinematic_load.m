wue=[02 03 06 09 10 11];
fileloc='C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\KinematicData';
for i = 1:length(wue)
    if i<5
        fileloc_temp=[fileloc filesep 'wue0' num2str(wue(i))];
    else
        fileloc_temp=[fileloc filesep 'wue' num2str(wue(i))];
    end
    txtname=dir([fileloc_temp filesep 'grasp*.txt']);
    for j=1:size(txtname,1)
        txt_data=importdata([fileloc_temp filesep txtname(j).name]);
        Data=txt_data.data;
        AllKinData(i).sessno(j).kindata=Data;
        %save([fileloc_temp filesep 'grasp' num2str(j)],'Data');
        %plot(Data(:,2))
        %print([fileloc filesep 'wue' num2str(wue(i)) 'grasp' num2str(j)],'-djpeg','-r300')
    end
end
save([fileloc filesep 'AllKinData'],'AllKinData');
