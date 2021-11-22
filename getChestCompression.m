%GETCHESTCOMPRESSION - generate a new MC realization of a chest compression
%with the standard parameters set forth by Joe
%
% Syntax:
%   [ct,vt] = getChestCompression(fs,T,N,sigN2,snrFlag)
%
% Inputs:
%
%  Requried:
%   None
%
%
%
% Optional:
%   fs :(double scalar) - default=100
%       Sampling frequency
%
%   N : (double Scalar) - default = -1
%       Defines the number of elements in the output signals
%
%   T : (double scalar) - default = 1
%       Defines the length of time the signal ct and vt will run to.
%       If N has an integervalue greater than 0 then N will be used to
%       define the length of the signal instead of T.
%
%   sigN2: (double scalar) - default = 0
%       Defines the power of the noise to be added to the compression and
%       velocity vectors. Noise is defiend as sqrt(sigN2)*randn(size)
% 
%   snrFlag : (logical scalar) - defulat = 0
%       If snrFlag is 1, then sigN2 will be interpreted as the desired SNR
%       level where 20 menas the signal is 20 dB louded than the noise.
%       This is done by setting leaving the amplitude of A as 1 but setting
%       the noise to power -20 dB.
%
% Outputs:
%
%
%   ct : (double Mx1 vector)
%       The chest compression displacement signal with noise added
% 
%   vt : (double Mx1 vector)
%       The chest compression velocity signal with noise added
%
% Written 4/24/2021 Joe Coffin


function [t,ct,vt,fr,fh] = getChestCompression(fs,T,N,sigN2,snrFlag,ctOffset,respHeight,heartRatio,t)



arguments
    fs (1,:) {mustBeNumeric,mustBeReal} = 100
    T (1,:) {mustBeNumeric,mustBeReal} = 1
    N (1,:) {mustBeInteger,mustBeReal} = -1
    sigN2(1,:) {mustBeNumeric,mustBeReal} = 0
    snrFlag(1,:)= 0
    ctOffset(1,:) {mustBeNumeric,mustBeReal} = 0
    respHeight (1,:) {mustBeNumeric,mustBeReal} = 1
    heartRatio (1,:) {mustBeNumeric,mustBeReal} = .3
    t (1,:) {mustBeNumeric,mustBeReal} = []
end

a = 15/60;
b = 18/60;
fr = a + (b-a).*rand;
%Set up fr (respiration rate) to be a random number uniformly generated between the bounds
% a = 15/60 to b = 18/60. This is the range 15 to 18 breaths per minute
% converted to Hz




a = 1;
b = 100/60;
fh = a + (b-a).*rand;
% set up fh (heart rate) to be a random number unifromly generated between
% a = 1 and b = 100/60. This sets the heart beat ragnge from 60 to 100 bpm
% (1 to 100/60) Hz. This represents a wide range of standing heart beats


ar = respHeight;
% ah = respHeight*heartRatio;
ah = heartRatio;
%ar is the amplitude of the respitory siganl
%ah is the ampllitude of the heart beat signal.
% ar and ah are chosen such that ah is much lower than ar

%Get Random phase between 0 and 360 degrees
a=0; %Phase min
b=2*pi;%phase max
respPhase = a + (b-a).*rand; %Random Phase for respitory signal
if isempty(t)
    if N == -1
        t = 0:1/fs:T;
    else
        t = 0:1/fs:(N-1)/fs;
    end
end
%Set up the sampling rate and the time vector. I choose fs as 100 since we
%are sampling a heart beat which is much less than 100Hz (ideally). Fs
%could probably be lowered in practice.

%set up noise vector to be added to ct and vt;
%if the SNR flag is up then that means sigN2 is a dB value that determines
%the power of the noise. 20 means the noise power should be -20.
%So we convert sigN2 to be the power in linear scale instead of dB. 
if snrFlag
%    sigN2 =  10^(-(sigN2/20)); %this is to find the sigma value, not the sigma squared value

   sigN2 = 10^(-sigN2/10); % because I am square rooting sigN2 on the next bit of code, I need to use /10 not /20. 
end
nt = sqrt(sigN2)*randn(size(t));

ct = ar*sin(2*pi*fr*t+respPhase) + ah*sin(2*pi*fh*t) + nt + ctOffset;

%For some simulations we need the velocity of the compression signal.
% Simple dearitive will work here
vt = 2*pi*fr*ar*cos(2*pi*fr*t) + ah*2*pi*fh*cos(2*pi*fh*t)+ nt;

end

