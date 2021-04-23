%EASYFFT Calculate and format the FFT centered at zero and with Proper w
%vec
%
% Syntax:
%   [fftOut,omega] = easyFFT(timeSig,fs,Npt)
% 
% Inputs:
% 
%  Requried:
%   timeSig: (double vector)
%     Signal to take the fft of
% 
%   Npt : (double Scalar)
%   	Size of the FFT. Must be larger than length of signal + 1
%      Higher values create better graphs
% 
% 
% 
% Outputs:
%   
%
% Written 7/29/19 Joe Coffin

function [fftOut,f] = easyFFT(timeSig,fs,Npt)
if nargin<3
    Npt = 2^nextpow2(length(timeSig));
end
fftOut = fftshift(fft(timeSig,Npt));
% omega = linspace(-pi,pi,Npt);

nSamples = length(fftOut);
df = fs / nSamples;
f = 0 : df : df*(nSamples-1);
f = fftshift(f);
ind = f>=fs/2-eps(fs/2);
f(ind) = f(ind)-fs;

end

