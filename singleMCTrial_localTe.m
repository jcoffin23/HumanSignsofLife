function [frError,fhError,rangeRes,chestSig,peakFFT,tgtStructOut,recvPow] = singleMCTrial_localTe(inStruct,tgtStruct)
%% Constants
fc = inStruct.fc;
bw = inStruct.bw;
recvArray = inStruct.recvArray;
usingArray = inStruct.usingArray;
sendElmnt = inStruct.sendElmnt;
respHeight = inStruct.respHeight;
heartHeight = inStruct.heartHeight;
beatError = inStruct.beatError;
clutterRCS = inStruct.clutterRCS;
bistatic = inStruct.bistatic;
stv = inStruct.stv;
% targetOffset = inStruct.targetOffset;
addNoise = inStruct.addNoise;
numRecv = inStruct.numRecv;
calcRangePlot = inStruct.enablePic;
miliPow = inStruct.miliPow;
rxRad = inStruct.rxRad;
numTgt = tgtStruct.numTgt;
feet2meter  = 0.3048;% 1 foot is this many meters.


numHuTgt = sum( tgtStruct.human ==1);%Find out how many targets are human (how many phase points to extract)
%% Setting up parameters
% rng(2017);
% fc = 5e9; %% 2 Ghz
c = 3e8;
lambda = c/fc;
range_max = 50;
tm = 5.5*range2time(range_max);
% tm = 4*10^-3
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





%% Waveform setup

showPlot = 0;

waveform = phased.FMCWWaveform('SweepTime',tm,'SweepBandwidth',bw,...
    'SampleRate',fs);
sig = waveform();

rangeRes = bw2range(bw); %Calculate the Range Resolution for the input bandwidth
%% Set up range Vector

L = size(sig,1);
freq_res = fs/L; %The freqeucny resolution is depentant on the sampling frequency and number of FFT points
freqGrid = (0:L-1).'*freq_res; %This is the beat frequency vector. Need to turn this into range

rangeVec = beat2range(freqGrid,sweepSlope);

[theta,range]=cart2pol(tgtStruct.offset,tgtStruct.yoffset); %Calculate range for targetIndex calcualtion
[~,targetIndex] = min(abs(rangeVec - range));
% targetIndex = [1:numHuTgt;targetIndex(tgtStruct.human ==1)]

%% Radar Setup
antAperture =6.06e-4;                        % Antenna aperture (m^2)

% frq = 77e9;
% lambda = c/frq;
% antGain = aperture2gain(antAperture,lambda);  % Antenna gain (dB)
% 
% antGain = 3.2538;
antGain = inStruct.gain;


txPkPower = db2pow(5)*1e-3;                   % Tx peak power (W)
txPkPower = miliPow*1e-3;
txGain = antGain;                             % Tx antenna gain (dB)

rxGain = antGain;                             % Rx antenna gain (dB)
rxNF = 4.5;                                   % Receiver noise figure (dB)

% Waveform transmitter
transmitter = phased.Transmitter('PeakPower',txPkPower,'Gain',txGain);

if bistatic
    % Radiator for single transmit element
    radiator = phased.Radiator('Sensor',sendElmnt,'OperatingFrequency',fc,'Polarization','Combined');
    
    % Collector for receive array
    collector = phased.Collector('Sensor',recvArray,'OperatingFrequency',fc,'Polarization','Combined');
else
    % Radiator for single transmit element
    radiator = phased.Radiator('Sensor',sendElmnt,'OperatingFrequency',fc);
    
    % Collector for receive array
    collector = phased.Collector('Sensor',recvArray,'OperatingFrequency',fc);
    
end

% Receiver preamplifier
receiver = phased.ReceiverPreamp('Gain',rxGain,'NoiseFigure',rxNF,...
    'SampleRate',fs);

% %Set up radar transciever with Tx as single Dipole and Rx as UCA.
% radar = radarTransceiver('Waveform',waveform,'Transmitter',transmitter,...
%     'TransmitAntenna',radiator,'ReceiveAntenna',collector,'Receiver',receiver);

%Set up free space Propogation
channel = phased.FreeSpace(...
    'PropagationSpeed',physconst('LightSpeed'), ...
    'OperatingFrequency',fc,'TwoWayPropagation',false, ...
    'SampleRate',fs,'MaximumDistanceSource','Property','MaximumDistance',range_max);

txchannel = phased.FreeSpace('SampleRate',fs,...
    'OperatingFrequency',fc,'PropagationSpeed',c );
rxchannel = phased.FreeSpace('SampleRate',fs,...
    'OperatingFrequency',fc,'PropagationSpeed',c);

%% AntennaPlatform

rxPlatform = phased.Platform('InitialPosition',[0;-rxRad;0],'Velocity',[0;0;0],'OrientationAxesOutputPort',true);
txPlatform = phased.Platform('InitialPosition',[0;rxRad;0],'Velocity',[0;0;0],'OrientationAxesOutputPort',true);


%% Target movement setup




fsChest = 500;    
dt = 1/fsChest;  
% seconds per sample
StopTime = 8;
t = (0:dt:StopTime-dt)';

[targets,targetplatform,tgtStruct] = TargetTool(tgtStruct,respHeight,heartHeight,fc,dt,t,fsChest,bistatic);
frActual = tgtStruct.frActual;
fhActual = tgtStruct.fhActual;


tgtStruct.t = t;
%% BeamFormer

if usingArray
    beamformer= phased.PhaseShiftBeamformer('SensorArray',recvArray,...
        'PropagationSpeed',c,'OperatingFrequency',fc,'WeightsOutputPort',true,'DirectionSource','Input port');
else
    beamformer = false;
end
%% Simulation start

% 
% pedest = backscatterPedestrian( 'Height',1.8, ...
%     'OperatingFrequency',fc,'InitialPosition',[10;0;0], ...
%     'InitialHeading',180,'WalkingSpeed',0.01);
% 


%% Simulation Loop
Nsweep = length(t);

if inStruct.collectAllTargetPhase ==1 %This may be missing in some test files, just add the field at set to 0
    %             if only phase from humans is wanted. Otherwise set to 1.
    numPhaseExtractPoints = tgtStruct.numTgt;
else
    numPhaseExtractPoints = numHuTgt;
end

reflectedPhase = zeros(numPhaseExtractPoints,Nsweep);

recvPow = zeros(1,Nsweep);
if numRecv ==1
    reflectedSig = zeros(length(sig),Nsweep);
    
else
    reflectedSig = zeros(length(sig),numTgt,numRecv);
end
phaseExtractPoints = zeros(numHuTgt,1);
peakFFT= zeros(numHuTgt,1);
for m = 1:Nsweep
    
    if bistatic
        [rxSig,~,~,~] = bistaticChirp(targets,transmitter,receiver,txchannel,rxchannel,radiator,collector,beamformer,targetplatform,rxPlatform,txPlatform,dt,waveform,numTgt);
    else
        [rxSig,~,~,~] = monostaticChirp(targets,transmitter,receiver,txchannel,rxchannel,radiator,collector,beamformer,targetplatform,rxPlatform,txPlatform,dt,waveform,numTgt,tgtStruct);
    end
%     targetplatform(dt);
    recvPow(m) = (norm(rxSig)^2) / length(rxSig);
    if 1
        tpResponse = fft(rxSig);
        
        maxtp = max(abs(tpResponse));
        noiseMult = db2pow(-stv)*maxtp;
        
        rxSig =  rxSig + noiseMult*randn(size(rxSig));
    end
    
    
     %Beam form the response
    if numRecv ==1
        reflectedSig(:,m) = rxSig;
        
    else
       reflectedSig(:,:,m) = rxSig;
    end
    reflectedFFT = fft(rxSig);
    
%     if m<100
%         [maxFFT,idMax] = max(abs(reflectedFFT));
%     end

    

    
        phaseExtractPoints = targetIndex + inStruct.binError; %Induce the range error. No need to do this every loop but its pretty here
        reflectedPhase(:,m) =  angle(reflectedFFT(phaseExtractPoints)); %Calculate the phase of the reflected signal AT the beat freqeuncy ONLY!!

    
    
    
end


genFig = 0;
if genFig
    
    L = size(reflectedSig,1);
    freq_res = fs/L; %The freqeucny resolution is depentant on the sampling frequency and number of FFT points
    freqGrid = (0:L-1).'*freq_res; %This is the beat frequency vector. Need to turn this into range
   
    rangeVec = beat2range(freqGrid,sweepSlope); %Turn beat freqeuncies into range (see documentation for how this works, simple multiply)
    figure
    plot(rangeVec,abs(reflectedFFT))
    title('FFT of Recieved Signal')
    ylabel('Magnitude')
    xlabel('Range (m)')
    xlim([0,40])
end



recvPow = mean(recvPow);

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
    for plotIdx = numHuTgt:-1:1
        plot(t,chestSig(plotIdx,:))
        hold on
        str{cnt}=['Target ',num2str(cnt)]
        cnt=cnt+1;
    end
    
    title(' Chest Compression c(t) vs Time')
    ylabel('c(t)')
    xlabel('time (s)')
    
    legend(str)
    
    
end
% ctEst = chestSig(tidx,:);



whichTarget = 1;
 %Which target to compare the frequencies
 
[fr3,fh3] = estFreqs(chestSig(whichTarget,:)-mean(chestSig(whichTarget,:)),fsChest); %Estimate the freqeuncy of the chest sig (Must be DC so zero frequency does not have as much power


tgtStructOut = tgtStruct;
% tgtStructOut.muFFT = muFFT;




frError = min(abs(fr3-frActual(whichTarget)));
fhError = min(abs(fh3-fhActual(whichTarget)));

if fhError > .5
   qq= 21; 
else
    qq=22;
end

%% Calculate Range Estimate plot (Only works with High Band width (.5 GHz)
% calcRangePlot = 0; %Only cacluate range Plot if needed. Do not need to for MC trials. It is only a diagonstic.
if calcRangePlot
    
    L = size(reflectedSig,1);
    freq_res = fs/L; %The freqeucny resolution is depentant on the sampling frequency and number of FFT points
    freqGrid = (0:L-1).'*freq_res; %This is the beat frequency vector. Need to turn this into range
    
    rangeVec = beat2range(freqGrid,sweepSlope); %Turn beat freqeuncies into range (see documentation for how this works, simple multiply)
    
%     q = 100000
%     rangeVec=rangeVec(1:200);
%     reflectedSigTemp = reflectedSig(1:q,:);
    reflectedSigTemp = reflectedSig;
    
%     [~,tidx] = min(abs(rangeVec -targetOffset)); %This is the idx of the range vec that corresponds to the target.
    
    win = hamming(size(reflectedSigTemp,1));
    rngDat = fft(reflectedSigTemp.*win,10000,1); %Calculated beat Frequencies using FFT
    
%     rngDat = rngDat(1:L,:); %Take only positive Freqencies
    
%     titleStr = ['clutter located at - (',num2str(tgtx),',',num2str(tgty),')'];
    
    figure
    imagesc(t,rangeVec,20*log10(abs(rngDat))) %Plot FFT data vs time and Beat range
    title(['Range vs time'])
    xlabel('Time (s)')
    % xlabel('Angle to Second Target')
    ylabel('Range Estimate')
    ylim([0,range_max/2])
%      ylim([0,10])
%     jct.util.SavetoPng(['.\temp\Range vs time ',titleStr,'.png'])

    %% Build Doppler Info
    showDop = 0;
    
    if showDop
        
        ds = 4000;
        N = 10*ds;
        TF = ds*tm;
        
        dopMax = lambda/(4*tm);
        dopRes = lambda/(2*N*tm);
        % Se page 40 of TI pdf that is linked below
        % https://training.ti.com/sites/default/files/docs/mmwaveSensing-FMCW-offlineviewing_4.pdf
        dopVec = -dopMax:dopRes:dopMax-dopRes;
        toDop = rngDat(:,1:ds);
        
        
        
        dopResp = fftshift(fft(toDop,N,2),2);
        
        
        figure
        imagesc(dopVec,rangeVec,abs(dopResp)) %Plot FFT data vs time and Beat range
        title('Range vs Doppler')
        xlabel('Doppler Speed (m/s)')
        ylabel('Range Estimate')
%         ylim([0,range_max])
%         xlim([-100,100])
    end
end

end

