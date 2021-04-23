function Hd = getRespFIlter(Fs)
%GETRESPFILTER Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.10 and DSP System Toolbox 9.12.
% Generated on: 19-Apr-2021 14:17:43

% Butterworth Lowpass filter designed using FDESIGN.LOWPASS.

% All frequency values are in Hz.
% Fs = 100;  % Sampling Frequency

N  = 6;    % Order
Fc = 0.3;  % Cutoff Frequency

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.lowpass('N,F3dB', N, Fc, Fs);
Hd = design(h, 'butter');

% [EOF]
