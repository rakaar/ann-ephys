function [bdata]=baselinezeroing_1ch(datachannel,base_samps)

[niters,nsamps]=size(datachannel);
bdata=zeros(niters,nsamps);
for jj=1:niters
    bdata(jj,:)=datachannel(jj,:)-mean(datachannel(jj,base_samps));
end
% %%%% design the notch filter at 50 hz
% fspec=fdesign.notch('N,F0,Q',8,50/(Fsi/2),25);
% d=design(fspec);
% %%%%
% fdata=filter(d,datachannel')';
end