
function [cref,fcref]=rng_ref(nfft,fs,pulsedur,slope)
%
% routine to compute ERS chirp and its fourier
% transform
%
% input
% fs - sampling frequency, ts=1./fs
% pulsedur - pulse duration
% slope - chirp slope
%
% set the constants and make npts be odd
%
npts=floor(fs*pulsedur);
ts=1./fs;
if(mod(npts,2.0) == 0.0)
npts=npts+1;
end
%
% compute the reference function
%
npt2=floor(npts/2.);
t=ts*(-npt2:npt2);
phase=pi*slope*t.*t;
cref1=exp(i*phase);
%
% pad the reference function to nfft
%
cref=[cref1,zeros(1,nfft-npts)]';
%
% compute the fourier transform
%
fcref=fft(cref)/nfft;
end