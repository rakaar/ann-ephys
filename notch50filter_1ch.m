function [fdata]=notch50filter_1ch(datachannel,Fsi)

[niters,nsamps]=size(datachannel);
%%%% design the notch filter at 50 hz
fspec=fdesign.notch('N,F0,Q',8,50/(Fsi/2),25);
d=design(fspec);
%%%%
fdata=filter(d,datachannel')';
clear datachannel
end