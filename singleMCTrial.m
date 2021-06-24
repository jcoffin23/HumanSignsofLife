function [frError,fhError,rangeRes] = singleMCTrial(fc,bw,rxRad,respHeight,heartRatio,beatError,clutterRCS,NumeleMents,thetaDist)
%% Constants

feet2meter  = 0.3048;% 1 foot is this many meters.
frError = nan
fhError = nan
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
sendElmnt = phased.IsotropicAntennaElement('BackBaffled',true);

% antenna = phased.CosineAntennaElement;
% fc = 1e9;
% pattern(antenna,fc,[-30:0.1:30],0,'Type','efield', ...
%     'CoordinateSystem','polar')



directionalInPut = 10;
antElmnt = phased.CosineAntennaElement('CosinePower',[1,directionalInPut]);
%




%


% Construct the receive array (UCA of two elements. One at -.5 the other at
% +.5. SO this is centered arround the TX array.
Ne = 1;
% rxRad = lambda/2;

rxRad = .5*feet2meter;

NumeleMents = 2;

recvArray = phased.ULA('Element',antElmnt,'NumElements',NumeleMents,...
    'ElementSpacing',rxRad);


% %
[pat,az,elv] = pattern(recvArray,fc,-180:.1:180,0);

plot(az,pat)
hold on
resp = antElmnt(fc,[0;0]);
antpat=pattern(antElmnt,fc,-180:.1:180,0)
plot(az,antpat)
title(['Beam pattern 20 Elements, Element spacing = ',num2str(rxRad),' meters'])
legend('Beamformed Response','Single antenna Response')
xlabel('Degrees')
ylabel('Response Magnitude (dB)')


% 
% 
% if rxRad ==.15
%     
%    jct.util.SavetoPng('HalffootSpacing2ELE.png') 
% end
% 
% if rxRad ==.03
%     
%    jct.util.SavetoPng('HalfWavelengthSpacing2ELE.png') 
% end
% 
%         if rxRad ==0.99 
%             jct.plotting.gif('spacingGif2Elemtns.gif','DelayTime',1,'frame',gcf)
%         else
%             jct.plotting.gif     
%         end







% return

% recvArray = phased.ULA('Element',antElmnt,'NumElements',6,...
%     'ElementSpacing',rxRad);
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
radiator = phased.Radiator('Sensor',sendElmnt,'OperatingFrequency',fc);

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


%% AntennaPlatform

antennaplatform = phased.Platform('InitialPosition',[0;0;0],'Velocity',[0;0;0]);



%% Clutter Platforms
targetRCS = .5;
numClutter = 2;

xyzPoses = [[10,0,0] ; [10,10,0]];
[clutter,clutterPlatform]=setupBackClutter(xyzPoses,numClutter,clutterRCS,fc);


%% Target movement setup
initTgtPos = [10; 0; 0]; %The initial position of the target

target = phased.RadarTarget('Model','Nonfluctuating','MeanRCS',targetRCS, ...
    'PropagationSpeed',physconst('LightSpeed'),'OperatingFrequency',fc);

x0 = 10;
y0 = 0;
z0 = 0;
fsChest = 2000;                   % samples per second
dt = 1/fsChest;                   % seconds per sample
StopTime = 8;
t = (0:dt:StopTime-dt)';
y = y0*ones(size(t));
z = z0*ones(size(t));

offset1 = initTgtPos(1);

sigN2 = 0; % db SNR

%Get the chest compression signal. We The frequencies are random each time
%and are drawn from a pool uniformly. We need to save them to compare the
%estimated freqeuncies to the actual ones. This is our metric!
[~,ct,vt,frActual,fhActual] = getChestCompression(fsChest,[],length(t),sigN2,0,offset1,respHeight,heartRatio);

%Put the movement data into the trajectory for the target platform.
wpts = [t';ct; y'; z']';
targetplatform = phased.Platform('MotionModel','Custom','CustomTrajectory',wpts);
offset2 = 17;


<<<<<<< Updated upstream
=======
%Put the movement data into the trajectory for the target platform.




phi = (180-thetaDist)/2; %base angle of Isocolese Triangle
alpha = 2*10*cosd(phi); %Base of Iso Triangle, distance between the two pints

xPoint = (alpha^2-200)/-20;
yPoint = sqrt(100-xPoint^2);
y2 = y + yPoint;

[~,ct2,vt2,frActual2,fhActual2] = getChestCompression(fsChest,[],length(t),sigN2,0,xPoint,respHeight,heartRatio);
wpts2 = [t';ct2; y2'; z']';
targetplatform2 = phased.Platform('MotionModel','Custom','CustomTrajectory',wpts2);

>>>>>>> Stashed changes

%% New Target Setup
 tgt1 = struct('Position',[0 5e3 0],'Velocity',[0 0 0]);tgt2 = struct('Position',[10e3 0 0],'Velocity',[0 0 0]);tgt = [tgt1 tgt2];





%% BeamFormer

beamformer= phased.PhaseShiftBeamformer('SensorArray',recvArray,...
    'PropagationSpeed',c,'OperatingFrequency',fc,'Direction',[0;0],'WeightsOutputPort',true);

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
% txpos = [-5,0,0]'

L = 1334;
freq_res = fs/L; %The freqeucny resolution is depentant on the sampling frequency and number of FFT points
freqGrid = (0:L-1).'*freq_res; %This is the beat frequency vector. Need to turn this into range

rangeVec = beat2range(freqGrid,sweepSlope)';
idmaxSave = zeros(1,Nsweep);







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
    
    
    if sum(clutterRCS) ~= 0
        [summedSig] = updateClutter(clutter,clutterPlatform,txpos,waveform,transmitter,radiator,channel,collector,receiver);
        tSig = rxsig + summedSig;
    else
        tSig = rxsig;
    end
    
    [beamfromedResponse,beamWeights]=beamformer(tSig);
    reflectedSig(:,m) = beamfromedResponse; %Beam form the response
    
    
    
    [tpfft]=fft(beamformer(rxsig));
    if m<100
    [maxFFT,idMax] = max(abs(tpfft));
    end
    
    
    
    
    if m==100
        idRunAvg =  mean(idmaxSave(1:m-1));
        idMax = round(idRunAvg);
    end
    
    idmaxSave(m) = idMax;
    reflectedFFT = fft(reflectedSig(:,m));
    
    
    
    idMax = idMax + beatError;
    if idMax<1
        idMax = 1;
    end
    if idMax>length(reflectedFFT)
        idMax = length(reflectedFFT);
    end
    reflectedPhase(m) =  angle(reflectedFFT(idMax)); %Calculate the phase of the reflected signal AT the beat freqeuncy ONLY!!
    
    
    
    
    
    
    
    %     plot(rangeVec,abs(reflectedFFT)/max(abs(reflectedFFT)))
    %     title('FFT of Response')
    %     xlabel('Distance')
    %     ylabel('Magnitude')
    %     xlim([0,25])
  
    %     if m==1
    %         jct.plotting.gif('fftGif.gif','DelayTime',1/8,'frame',gcf)
    %     elseif  mod(m,50) ==0;
    %         plot(rangeVec,abs(reflectedFFT)/max(abs(reflectedFFT)))
    %         title('FFT of Response')
    %         xlabel('Distance')
    %         ylabel('Magnitude')
    %         xlim([0,25])
    %         jct.plotting.gif
    %         if m == (50*200)
    %             q=1
    %         end
    %     end
    
    
    %Calculate Range Angle

    
    
    
end


% rngdopresp = phased.RangeDopplerResponse('PropagationSpeed',c,...
%     'DopplerOutput','Speed','OperatingFrequency',fc,'SampleRate',fs,...
%     'RangeMethod','FFT','SweepSlope',sweepSlope,...
%     'RangeFFTLengthSource','Property','RangeFFTLength',2048,...
%     'DopplerFFTLengthSource','Property','DopplerFFTLength',256);
%
% clf;
% plotResponse(rngdopresp,reflectedSig(:,1:200));
% ylim([0,25])



%% Calculate Phase and then Displacement

reflectedPhase = unwrap(reflectedPhase,[],2);
chestSig = lambda*reflectedPhase /(4*pi);



%% IQ dmodulation
%
% phi = unwrap(atan2(imag(reflectedSig),real(reflectedSig))); %Calcuate the phase of the reflected signal.
%
% phiSig = min(phi); %The min phase is the reflected signal.


%The chest signal is realted to the phase of the reflection (see source 1)
PlotChest = 0; %only plot chest signal if this is a 1. This is a debug thing
if PlotChest
    cnt=1;
    str = {};
    figure
<<<<<<< Updated upstream
    plot(t,chestSig)
=======
    for plotIdx = numTgt:-1:1
    plot(t,chestSig(plotIdx,:))
    hold on
    str{cnt}=['Target ',num2str(cnt)]
    cnt=cnt+1;
    end
>>>>>>> Stashed changes
    title('Phase signal Converted to c(t) vs Time')
    ylabel('c(t)')
    xlabel('time (s)')
    
<<<<<<< Updated upstream
=======
  legend(str)
>>>>>>> Stashed changes
end
% ctEst = chestSig(tidx,:);




[fr3,fh3] = estFreqs(chestSig(3,:)-mean(chestSig(3,:)),fsChest); %Estimate the freqeuncy of the chest sig (Must be DC so zero frequency does not have as much power

% 
% [fr2,fh2] = estFreqs(chestSig(2,:)-mean(chestSig(2,:)),fsChest);
% [fr1,fh1] = estFreqs(chestSig(1,:)-mean(chestSig(1,:)),fsChest);
% frTotal = [fr1,fr2,fr3];
% fhTotal = [fh1,fh2,fh3];
frError = min(abs(fr3-frActual));
fhError = min(abs(fh3-fhActual));

%% Calculate Range Estimate plot (Only works with High Band width (.5 GHz)
calcRangePlot = 0; %Only cacluate range Plot if needed. Do not need to for MC trials. It is only a diagonstic.
if calcRangePlot
    L = size(reflectedSig,1);
    freq_res = fs/L; %The freqeucny resolution is depentant on the sampling frequency and number of FFT points
    freqGrid = (0:L/2).'*freq_res; %This is the beat frequency vector. Need to turn this into range
    
    rangeVec = beat2range(freqGrid,sweepSlope); %Turn beat freqeuncies into range (see documentation for how this works, simple multiply)
    
    [~,tidx] = min(abs(rangeVec -initTgtPos(1))); %This is the idx of the range vec that corresponds to the target.
    
    
    rngDat = fft(reflectedSig,L,1); %Calculated beat Frequencies using FFT
    
    rngDat = rngDat(1:L/2+1,:); %Take only positive Freqencies
    
%     
    figure
    imagesc(t,rangeVec,20*log10(abs(rngDat))) %Plot FFT data vs time and Beat range
    
    
    title(['Range vs time (',num2str(NumeleMents), 'Antenna Elements)'])
    xlabel('Time (s)')
    ylabel('Range Estimate')
    ylim([0,25])
    % caxis([-10,5])
    % caxis([0,10])
<<<<<<< Updated upstream
end



=======
% jct.util.SavetoPng(['Range vs time (',num2str(NumeleMents), ' Antenna Elements).png'])

%% Build Doppler Info 
showDop = 0;
>>>>>>> Stashed changes








end

