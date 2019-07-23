function [ freqtime ] = wavelet_morlet( data,frex,wuedetails,events_onoffsetkin_data)

baselinetime=diff(events_onoffsetkin_data(wuedetails(1)).sessno(wuedetails(2)).data(1:2,:));
baselinetime=baselinetime(wuedetails(3));

srate = 1000;
wavtime = -0.5:1/srate:0.5; 



nFrex=length(frex);
nData = length(data); 
nKern = length(wavtime); 
nConv = nData+nKern-1;
halfwav=floor(nKern/2)+1;

s = linspace(8,20,nFrex) ./ (2*pi.*frex); % here we are varying the n values experiment with this
%sigX = fft(data,nConv); 
freqtime = zeros(nFrex,length(data)); 

for fi=1:nFrex
    nData_mirror=floor(3*srate/frex(fi));% mirroring 3 cycles worth of data at start and end 
    data=[data(nData_mirror:-1:1) data data(end:-1:(end-nData_mirror+1))];
    sigX = fft(data,(nConv+2*nData_mirror));
    
    cmw = exp(2*1i*pi*frex(fi)*wavtime + ... 
        -(wavtime.^2)/(2*s(fi)^2)); 
    cmwX = fft(cmw,(nConv+2*nData_mirror)); 
    cmwX = cmwX ./ max(cmwX); %normalizing each spec
    as = ifft(sigX .* cmwX);
    freqtime(fi,:) = abs(as((halfwav+nData_mirror):end-halfwav+1-nData_mirror)).^2; 
end

% dB normalizing 
baselinepwr=mean(freqtime(:,(1:baselinetime+1)),2);
baselinepwr=repmat(baselinepwr,1,nData);
freqtime=10*log10(freqtime./baselinepwr);
end

