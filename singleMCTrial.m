function [frError,fhError,rangeRes] = singleMCTrial(fc,bw,rxRad,respHeight,heartRatio)


%% Setting up parameters
% rng(2017);
% fc = 5e9; %% 2 Ghz
c = 3e8;
lambda = c/fc;
range_max = 50;
tm = 5.5*range2time(range_max);
% rangeRes  = .3;
% bw = range2bw(rangeRes,c); %0.5 Ghz
% bw = .18*10^9
sweepSlope = bw/tm;
fr_max = range2beat(range_max,sweepSlope,c);
v_max = 10*1000/3600;
fd_max = speed2dop(2*v_max,lambda);
fb_max = fr_max+fd_max;
fst = max(2*fb_max,bw);


fs = ceil(tm*fst)/tm; %Set fs to be interger multiple of 

showPlot = 0;

waveform = phased.FMCWWaveform('SweepTime',tm,'SweepBandwidth',bw,...
    'SampleRate',fs);


rangeRes = bw2range(bw); %Calculate the Range Resolution for the input bandwidth
%% Array
% Model the antenna element
antElmnt = phased.IsotropicAntennaElement('BackBaffled',true);

% Construct the receive array (UCA of two elements. One at -.5 the other at
% +.5. SO this is centered arround the TX array.
Ne = 1;
% rxRad = 0.5;
recvArray = phased.UCA('Element',antElmnt,'NumElements',2,...
    'Radius',rxRad);

% Half-power beamwidth of the receive array
% hpbw = beamwidth(rxArray,fc,'PropagationSpeed',c);

antAperture =6.06e-4;                        % Antenna aperture (m^2)
antGain = aperture2gain(antAperture,lambda);  % Antenna gain (dB)

txPkPower = db2pow(5)*1e-3;                   % Tx peak power (W)
txGain = antGain;                             % Tx antenna gain (dB)

rxGain = antGain;                             % Rx antenna gain (dB)
rxNF = 4.5;                                   % Receiver noise figure (dB)

% Waveform transmitter
transmitter = phased.Transmitter('PeakPower',txPkPower,'Gain',txGain);

% Radiator for single transmit element
radiator = phased.Radiator('Sensor',antElmnt,'OperatingFrequency',fc);

% Collector for receive array
collector = phased.Collector('Sensor',recvArray,'OperatingFrequency',fc);

% Receiver preamplifier
receiver = phased.ReceiverPreamp('Gain',rxGain,'NoiseFigure',rxNF,...
    'SampleRate',fs);

%Set up radar transciever with Tx as single Dipole and Rx as UCA. 
radar = radarTransceiver('Waveform',waveform,'Transmitter',transmitter,...
    'TransmitAntenna',radiator,'ReceiveAntenna',collector,'Receiver',receiver);

%Set up free space Propogation
channel = phased.FreeSpace(...
    'PropagationSpeed',physconst('LightSpeed'), ...
    'OperatingFrequency',fc,'TwoWayPropagation',false, ...
    'SampleRate',fs,'MaximumDistanceSource','Property','MaximumDistance',range_max);



%% Target Def

target = phased.RadarTarget('Model','Nonfluctuating','MeanRCS',0.5, ...
    'PropagationSpeed',physconst('LightSpeed'),'OperatingFrequency',fc);

antennaplatform = phased.Platform('InitialPosition',[0;0;0],'Velocity',[0;0;0]);


initTgtPos = [10; 0; 0]; %The initial position of the target

beamformer = phased.PhaseShiftBeamformer('SensorArray',recvArray,...
    'PropagationSpeed',c,'OperatingFrequency',fc,'Direction',[0;0]);


%% Target movement setup

x0 = 10;
y0 = 0;
z0 = 0;   
fsChest = 2000;                   % samples per second
dt = 1/fsChest;                   % seconds per sample
StopTime = 20; 
t = (0:dt:StopTime-dt)'; 
y = y0*ones(size(t));
z = z0*ones(size(t));

offset = initTgtPos(1);

sigN2 = 0; % db SNR

%Get the chest compression signal. We The frequencies are random each time
%and are drawn from a pool uniformly. We need to save them to compare the
%estimated freqeuncies to the actual ones. This is our metric!
[~,ct,vt,frActual,fhActual] = getChestCompression(fsChest,[],length(t),sigN2,0,offset,respHeight,heartRatio);

%Put the movement data into the trajectory for the target platform. 
wpts = [t';ct; y'; z']'; 
targetplatform = phased.Platform('MotionModel','Custom','CustomTrajectory',wpts);
%% Simulation start


Nsweep = length(t);

txpos = antennaplatform.InitialPosition;

Nft = round(waveform.SampleRate*tm);

%  Allocate array for received echoes
weights = ones(Ne,1);
w = 1/Ne * weights;

% w = hamming(Ne)



TchestSample = 1/fsChest;

reflectedSig = zeros(Nft,Nsweep);
cc=1;
reflectedPhase = zeros(1,Nsweep);

for m = 1:Nsweep
    
    % Update the target position
    [tgtpos,tgtvel] = targetplatform(TchestSample);
    % Get the range and angle to the target
    [tgtrng,tgtang] = rangeangle(tgtpos,txpos);
    % Generate the pulse
    refSig = waveform();
    % Transmit the pulse. Output transmitter status
    [sig] = transmitter(refSig);
    % Radiate the pulse toward the target
    sig = radiator(sig,tgtang);
    % Propagate the pulse to the target in free space
    sig = channel(sig,txpos,tgtpos,[0;0;0],tgtvel);
    % Reflect the pulse off the target
    sig = target(sig);
    % Propagate the echo to the antenna in free space
    sig = channel(sig,tgtpos,txpos,tgtvel,[0;0;0]);
    % Collect the echo from the incident angle at the antenna
    sig = collector(sig,tgtang);
    % Receive the echo at the antenna when not transmitting
    rxsig = receiver(sig);
    rxsig = dechirp(rxsig,refSig); %This mixxes the tx signal with the RX signal. the frequency of this is what gives range
    %     y = (w'*rxsig.').'; %This is my attempt at beamforming. This is for
    %     ULA though and we are now dealing with UCA. Beamforming is slightly
    %     different and is ignored for now.
    
    reflectedSig(:,m) = beamformer(rxsig); %Beam form the response
    
    
    reflectedPhase(m) =  angle(max(fft(reflectedSig(:,m)))); %Calculate the phase of the reflected signal AT the beat freqeuncy ONLY!!
end



%% Calculate Phase and then Displacement

reflectedPhase = unwrap(reflectedPhase);
chestSig = lambda*reflectedPhase /(4*pi);



%% IQ dmodulation
% 
% phi = unwrap(atan2(imag(reflectedSig),real(reflectedSig))); %Calcuate the phase of the reflected signal.
% 
% phiSig = min(phi); %The min phase is the reflected signal.


 %The chest signal is realted to the phase of the reflection (see source 1)
PlotChest = 0; %only plot chest signal if this is a 1. This is a debug thing
if PlotChest
    
    figure
    plot(t,chestSig)
    title('Phase signal Converted to c(t) vs Time')
    ylabel('c(t)')
    xlabel('time (s)')
end
% ctEst = chestSig(tidx,:);




[fr,fh] = estFreqs(chestSig-mean(chestSig),fsChest); %Estimate the freqeuncy of the chest sig (Must be DC so zero frequency does not have as much power

frError = abs(fr-frActual);
fhError = abs(fh-fhActual);

%% Calculate Range Estimate plot (Only works with High Band width (.5 GHz)
calcRangePlot = 0; %Only cacluate range Plot if needed. Do not need to for MC trials. It is only a diagonstic.
if calcRangePlot
    L = 2048;
    freq_res = fs/L; %The freqeucny resolution is depentant on the sampling frequency and number of FFT points
    freqGrid = (0:L/2).'*freq_res; %This is the beat frequency vector. Need to turn this into range
    
    rangeVec = beat2range(freqGrid,sweepSlope); %Turn beat freqeuncies into range (see documentation for how this works, simple multiply)
    
    [~,tidx] = min(abs(rangeVec -initTgtPos(1))); %This is the idx of the range vec that corresponds to the target.
    
    
    rngDat = fft(reflectedSig,L,1); %Calculated beat Frequencies using FFT
    
    rngDat = rngDat(1:L/2+1,:); %Take only positive Freqencies
    
    
    figure
    imagesc(t,rangeVec,20*log10(abs(rngDat))) %Plot FFT data vs time and Beat range
    
    
    title('Range vs time')
    xlabel('Time (s)')
    ylabel('Range Estimate')
    ylim([0,25])
    % caxis([-10,5])
    % caxis([0,10])
end











end
