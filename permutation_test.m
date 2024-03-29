function [ pval1,testStatistic ] = permutation_test( dat1, dat2, no_iters )
%this function computes the pvals for the difference in the dat2 by
%generating a new statistic: the difference betw each data points of dat1
%and dat2. so the sample will be 1xlength(dat1). the population against
%which it computes ANOVA (surrogate data points) is generated by finding
%differences between a random point from dat1 and dat2.. which it computes
%10k times using permutation. 

datainp=dat1-dat2;
datajoin=[dat1(:);dat2(:)];

meandiff=mean(datainp);

permdata=zeros(no_iters,1);
for i=1:no_iters
    permid=randperm(length(datajoin));
    datajoin_temp=datajoin(permid);
    permdata(i,1)=mean(datajoin_temp(1:(end/2))-datajoin_temp(((end/2)+1):end));    
end

testStatistic=permdata;clearvars permdata; 
std_testStatistic=std(testStatistic);
zscore=(meandiff-mean(testStatistic))/std_testStatistic;
pval1=normcdf(zscore);
pval1=1-pval1;

end

