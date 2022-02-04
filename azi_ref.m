function [cazi,fcazi]=azi_ref(nazi,PRF,fdc,fr)
%
% routine to compute ERS azimuthal chirp and its
% fourier transform
%
% input
% nazi - number of points in azimuth
% PRF - pulse repitition frequency, ts=1./fs
% fdc - doppler centriod frequency
% fr - doppler frequency rate
%
% set the constants and make npts be odd
%
npts=min(nazi-1,1296);
ts=1./PRF;
if(mod(npts,2.0) == 0.0)
    npts=npts+1;
end
%
% compute the azimuth chirp function
%
npt2=floor(npts/2.);
t=ts*(-npt2:npt2);
phase=-2.*pi*fdc*t+pi*fr*t.*t;
cazi1=exp(i*phase);
%
% pad the function to nfft
%
cazi=[cazi1(npt2:npts),zeros(1,nazi-npts),cazi1(1:npt2-1)]';
%
% compute the fourier transform
%
fcazi=fft(cazi)/nazi;
end
