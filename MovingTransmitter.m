
%% 
c = 3e8;

feet2meter  = 0.3048;% 1 foot is this many meters.
%% Radar Params
fc = 5e9; %Carrier Freqeuncy

lambda = c/fc; %Wave Length
bw = .5e9; %BandWidth
miliPow = 10; %Tx Power

rxRad = .5*feet2meter; %How far apart Tx and Rx are


%Respitory Height is about 4 - 12 mm
% Heart heightis about .2-.5mm
heartHeight = .5*10^-3; %Heart rate is 30 percent of the respitory rate.
respHeight = 10 * 10^-3;

clutterRCS = .5; % Mean RCS of each background target
bistatic = 0

% NumeleMents = 2;
%
%
%

directionalInPut = 10;
antElmnt = phased.CosineAntennaElement('CosinePower',[1,directionalInPut]);

recvArray = antElmnt;
% recvArray = phased.ULA('Element',antElmnt,'NumElements',NumeleMents,...
%     'ElementSpacing',rxRad);

% recvArray = phased.ShortDipoleAntennaElement('AxisDirection','Y');
% recvArray=phased.CrossedDipoleAntennaElement;
usingArray = 0; %Set to 0 because array is element not actual array
%




% 
% [pat,az,elv] = pattern(recvArray,fc,-180:.1:180,0);
% % figure
% plot(az,pat)
% hold on
% % resp = antElmnt(fc,[0;0]);
% antpat=pattern(sendElmnt,fc,-180:.1:180,0);
% plot(az,antpat)
% title(['Beam pattern 20 Elements, Element spacing = ',num2str(rxRad),' meters'])
% legend('Beamformed Response','Single antenna Response')
% xlabel('Degrees')
% ylabel('Response Magnitude (dB)')


%%

numRecv = 1 %Number of rxs
directionalInPut = 10
sendElmnt = phased.CosineAntennaElement('CosinePower',[1,directionalInPut]);
% sendElmnt = phased.ShortDipoleAntennaElement('AxisDirection','Y');

%% Target Def
numTgt = 1;


numSamps = 10000;
tgtStruct.numTgt = numTgt;

tgtStruct.human = [1];
tgtStruct.offset = [10];
tgtStruct.yoffset= [0];
tgtStruct.zoffset= [1];
tgtStruct.fhActual=zeros(numSamps,numTgt);
tgtStruct.frActual=zeros(1,numTgt);


meanHumanRCS = 0.1992;


%% Constants
% fc = inStruct.fc;
% bw = inStruct.bw;
% recvArray = inStruct.recvArray;
% usingArray = inStruct.usingArray;
% sendElmnt = inStruct.sendElmnt;
% respHeight = inStruct.respHeight;
% heartHeight = inStruct.heartHeight;
% beatError = inStruct.beatError;
% clutterRCS = inStruct.clutterRCS;
% bistatic = inStruct.bistatic;
% stv = inStruct.stv;
% % targetOffset = inStruct.targetOffset;
% addNoise = inStruct.addNoise;
% numRecv = inStruct.numRecv;
% calcRangePlot = inStruct.enablePic;
% miliPow = inStruct.miliPow;
% rxRad = inStruct.rxRad;
% numTgt = tgtStruct.numTgt;


%% Setting up parameters
% rng(2017);
% fc = 5e9; %% 2 Ghz




range_max = 50;
tm = 5.5*range2time(range_max);
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


[~,targetIndex] = min(abs(rangeVec - tgtStruct.offset));


%% Radar Setup
antAperture =6.06e-4;                        % Antenna aperture (m^2)
% 
% frq = 77e9
% lambda = c/frq
% antGain = aperture2gain(antAperture,lambda);  % Antenna gain (dB)

antGain = 3.2538;

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
    'OperatingFrequency',fc,'PropagationSpeed',c);
rxchannel = phased.FreeSpace('SampleRate',fs,...
    'OperatingFrequency',fc,'PropagationSpeed',c);

%% AntennaPlatform

rxPlatform = phased.Platform('InitialPosition',[0;-rxRad;0],'Velocity',[0;0;.1],'OrientationAxesOutputPort',true);
txPlatform = phased.Platform('InitialPosition',[0;rxRad;0],'Velocity',[0;0;.1],'OrientationAxesOutputPort',true);


%% Target movement setup




fsChest = 500;    
dt = 1/fsChest;  
% seconds per sample
StopTime = 20;
t = (0:dt:StopTime-dt)';

[targets,targetplatform,tgtStruct] = TargetTool(tgtStruct,respHeight,heartHeight,fc,dt,t,fsChest,bistatic);
frActual = tgtStruct.frActual;
fhActual = tgtStruct.fhActual;



%% BeamFormer

if usingArray
    beamformer= phased.PhaseShiftBeamformer('SensorArray',recvArray,...
        'PropagationSpeed',c,'OperatingFrequency',fc,'WeightsOutputPort',true,'DirectionSource','Input port');
else
    beamformer = false;
end
%% Simulation start


pedest = backscatterPedestrian( 'Height',1.8, ...
    'OperatingFrequency',fc,'InitialPosition',[10;0;0], ...
    'InitialHeading',180,'WalkingSpeed',0.01);



%% Simulation Loop
Nsweep = length(t);
reflectedPhase = zeros(numTgt,Nsweep);

recvPow = zeros(1,Nsweep);



reflectedSig = zeros(length(sig),numRecv);

for m = 1:Nsweep
    
    if bistatic
        [rxSig,~,~,~] = bistaticChirp(targets,transmitter,receiver,txchannel,rxchannel,radiator,collector,beamformer,targetplatform,rxPlatform,txPlatform,dt,waveform,numTgt);
    else
        [rxSig,~,~,~] = monostaticChirp(targets,transmitter,receiver,txchannel,rxchannel,radiator,collector,beamformer,targetplatform,rxPlatform,txPlatform,dt,waveform,numTgt);
    end
%     targetplatform(dt);
    recvPow(m) = (norm(rxSig)^2) / length(rxSig);
    
    reflectedSig(:,m) = rxSig; %Beam form the response
    
    reflectedFFT = fft(rxSig);
    
    if m<100
        [maxFFT,idMax] = max(abs(reflectedFFT));
         idmaxSave(m) = idMax;
    end

    
    if m==100
        idRunAvg =  mean(idmaxSave(1:m-1));
        idMax = round(idRunAvg);
        
        meanfft = fft(mean(reflectedSig(:,1:m),2));
        fftAtTarget = abs(meanfft(targetIndex));
        muFFT = mean(abs(meanfft));
        [pks,locs] = findpeaks(abs(meanfft));
        [pks,sortIdx]=sort(pks);
        locs=locs(sortIdx);
        
        pks = pks(end-2:end);
%         phaseExtractPoints = locs(end-1:end);
        phaseExtractPoints = locs(end);
        peakFFT = pks(end);
    end
    
    if m>100
        reflectedPhase(:,m) =  (angle(reflectedFFT(phaseExtractPoints))); %Calculate the phase of the reflected signal AT the beat freqeuncy ONLY!!
    end
    
    
    
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
    plot(t(101:end),chestSig(101:end))
    for plotIdx = numTgt:-1:1
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





frError = min(abs(fr3-frActual(whichTarget)));
fhError = min(abs(fh3-fhActual(whichTarget)));

if fhError > .5
   qq= 21; 
else
    qq=22;
end

%% Calculate Range Estimate plot (Only works with High Band width (.5 GHz)
calcRangePlot = 0; %Only cacluate range Plot if needed. Do not need to for MC trials. It is only a diagonstic.
if calcRangePlot
    L = size(reflectedSig,1);
    freq_res = fs/L; %The freqeucny resolution is depentant on the sampling frequency and number of FFT points
    freqGrid = (0:L-1).'*freq_res; %This is the beat frequency vector. Need to turn this into range
    
    rangeVec = beat2range(freqGrid,sweepSlope); %Turn beat freqeuncies into range (see documentation for how this works, simple multiply)
    
%     [~,tidx] = min(abs(rangeVec -targetOffset)); %This is the idx of the range vec that corresponds to the target.
    
    
    rngDat = fft(reflectedSig,L,1); %Calculated beat Frequencies using FFT
    
    rngDat = rngDat(1:L,:); %Take only positive Freqencies
    
%     titleStr = ['clutter located at - (',num2str(tgtx),',',num2str(tgty),')'];
    
    figure
    imagesc(t,rangeVec,20*log10(abs(rngDat))) %Plot FFT data vs time and Beat range
    title(['Range vs time'])
    xlabel('Time (s)')
    % xlabel('Angle to Second Target')
    ylabel('Range Estimate')
    ylim([0,range_max + 10])

%     jct.util.SavetoPng(['.\temp\Range vs time ',titleStr,'.png'])

    %% Build Doppler Info
    showDop = 0;
    
    if showDop
        
        ds = 5000;
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
        ylim([0,range_max])
        xlim([-100,100])
    end
end

%%

pltYes = 0;
if pltYes
    to = txPlatform.InitialPosition;
    ro = rxPlatform.InitialPosition;
    
    posmat = targetplatform.CustomTrajectory';
    originPos = targetplatform.CustomTrajectory(1,2:end)
    figure(1)
    plot3(posmat(2,:),posmat(3,:),posmat(4,:),'b.')
    hold on
    plot3(ro(1),ro(2),ro(3),'rx')
    
    plot3(to(1),to(2),to(3),'rx')
    
    plot3(originPos(1),originPos(2),originPos(3),'ro')
    xlabel('X (m)')
    ylabel('Y (m)')
    zlabel('Z (m)')
    grid
    hold off
end

