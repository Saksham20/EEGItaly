function eegout=eegdata_struct_gen_new(data,eegout)



eegout.data=data;
eegout.pnts=size(data,2);
eegout.times=[1:eegout.pnts];
eegout.xmax=(eegout.pnts)/1000;

end