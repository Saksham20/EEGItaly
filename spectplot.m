function []=spectplot(wuelist,comment)

imgloc='C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\scripts\sasicaout\sasica\wicaOn\pipeout\norm\';

nowues=size(wuelist,2);
goodelectrodes=1:126;
load('C:\Users\PROVOST\Google Drive (sxs1790@case.edu)\MAT_files\eegexample.mat');

for i=1:nowues
    FigH = figure('Position', get(0, 'Screensize'));
    topoplot([],eegout.chanlocs(goodelectrodes), 'style', 'blank',  'electrodes', 'numpoint', 'chaninfo', eegout.chaninfo,'headrad',0.5,'emarker',{'.','k',1,1});
    topoplot([],eegout.chanlocs(goodelectrodes), 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', eegout.chaninfo,'headrad',0.5,'emarker',{'.','k',1,1});
    topoplot([],eegout.chanlocs(wuelist(i).betachans), 'style', 'blank',  'electrodes', 'on', 'chaninfo', eegout.chaninfo,'headrad',0.5,'emarker',{'o','r',25,3});
    print([imgloc sprintf('wueno%d_betamod_%s',i,comment)],'-djpeg');close all;
end

end
    
