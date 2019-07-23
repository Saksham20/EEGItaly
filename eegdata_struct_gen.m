function eegout=eegdata_struct_gen(data, goodelectrodes)
% this function takes eeg data matrix and electrode nos as input and
% outputs the EEGlab compliant data structure. 

load('C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\MAT_files\eegexample.mat');
eegout.data=data;
eegout.pnts=size(data,2);
eegout.times=[1:eegout.pnts];
eegout.xmax=(eegout.pnts)/1000;
eegout.nbchan=length(goodelectrodes);
eegout.chanlocs=eegout.chanlocs(goodelectrodes);
end