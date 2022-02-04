
%% 
c = 3e8;

feet2meter  = 0.3048;% 1 foot is this many meters.
%% Radar Params
fc = 26e9; %Carrier Freqeuncy

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
numTgt = 4;


numSamps = 10000;
tgtStruct.numTgt = numTgt;

tgtStruct.human = [1,0,0,0];
tgtStruct.offset = [10,15,12,20];
tgtStruct.yoffset= [0,0,0,0];
tgtStruct.zoffset= [1,0,0,0];
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

prf = 1/tm;
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

radarVelocity = 0;

rxPlatform = phased.Platform('InitialPosition',[0;-rxRad;-2],'Velocity',[0;0;radarVelocity],'OrientationAxesOutputPort',true);
txPlatform = phased.Platform('InitialPosition',[0;rxRad;-2],'Velocity',[0;0;radarVelocity],'OrientationAxesOutputPort',true);


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



reflectedSig = zeros(length(sig),Nsweep);

for m = 1:Nsweep
    
    if bistatic
        [rxSig,~,~,~] = bistaticChirp(targets,transmitter,receiver,txchannel,rxchannel,radiator,collector,beamformer,targetplatform,rxPlatform,txPlatform,dt,waveform,numTgt);
    else
        [rxSig,~,~,~] = monostaticChirp(targets,transmitter,receiver,txchannel,rxchannel,radiator,collector,beamformer,targetplatform,rxPlatform,txPlatform,dt,waveform,numTgt);
    end
%     targetplatform(dt);
    recvPow(m) = (norm(rxSig)^2) / length(rxSig);
    
    reflectedSig(:,m) = rxSig; %Beam form the response
    
%     reflectedFFT = fft(rxSig);   
    
end


%% Calculate Range Estimate plot (Only works with High Band width (.5 GHz)
calcRangePlot = 0; %Only cacluate range Plot if needed. Do not need to for MC trials. It is only a diagonstic.
if calcRangePlot
    
    pulseCompression = phased.RangeResponse('RangeMethod', 'FFT', 'PropagationSpeed', c, 'SampleRate', fs,'SweepSlope',sweepSlope,'ReferenceRangeCentered',0);
%     [cdata, rnggrid] = pulseCompression(reflectedSig);
    
    
    L = size(reflectedSig,1);
    freq_res = fs/L; %The freqeucny resolution is depentant on the sampling frequency and number of FFT points
    freqGrid = (0:L-1).'*freq_res; %This is the beat frequency vector. Need to turn this into range
    
    rangeVec = beat2range(freqGrid,sweepSlope); %Turn beat freqeuncies into range (see documentation for how this works, simple multiply)
    

    
    rngDat = fft(reflectedSig,L,1); %Calculated beat Frequencies using FFT
    
    rngDat = rngDat(1:L,:); %Take only positive Freqencies

    
    figure
    imagesc(t,rangeVec,20*log10(abs(rngDat))) %Plot FFT data vs time and Beat range
    title(['Range vs time'])
    xlabel('Time (s)')

    ylabel('Range Estimate')
    ylim([0,range_max + 10])



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


range_max = 50;
tm = 5.5*range2time(range_max);
sweepSlope = bw/tm;
fr_max = range2beat(range_max,sweepSlope,c);
v_max = 10*1000/3600;
fd_max = speed2dop(2*v_max,lambda);
fb_max = fr_max+fd_max;
fst = max(2*fb_max,bw);

initSig = reflectedSig
%% 
reflectedSig = initSig;



reflectedSig = reflectedSig.';
ncol = width(reflectedSig);
nrow = height(reflectedSig);

[cref,fcref]=rng_ref(ncol,fs,tm,sweepSlope);
%
% take the fft of the SAR data
%
fcdata=fft(reflectedSig);
%
% multiply by the range reference function
%
cout=0.*fcdata;
for k=1:nrow
    ctmp=fcdata(:,k);
    ctmp=fcref.*ctmp;
    cout(:,k)=ctmp;
end
clear cdata
%
% now take the inverse fft
%
odata=ifft(cout.');


figure
imagesc(t,rangeVec,20*log10(abs(odata.'))) %Plot FFT data vs time and Beat range
title(['Range vs time'])
xlabel('Time (s)')

ylabel('Range Estimate')


%% Az correction

fdc=.1;
fr=2*radarVelocity*radarVelocity/(10*lambda);

[cazi,fcazi]=azi_ref(nrow,prf,fdc,fr);
%
% take the column-wise fft of the range-compressed
%
fcdata=fft(odata.');
%
% multiply by the azimuth reference function
%
cout=0.*fcdata;
for k=1:ncol
    ctmp=fcdata(:,k);
    ctmp=fcazi.*ctmp;
    cout(:,k)=ctmp;
end
%
% now take the inverse fft and plot the data
%
odata=ifft(cout);

figure
imagesc(t,rangeVec,20*log10(abs(odata.'))) %Plot FFT data vs time and Beat range
title(['Range vs time'])
xlabel('Time (s)')

ylabel('Range Estimate')

%%
Rc = 10
fastTime = (0:1/fs:(length(sig)-1)/fs);

crossRangeResolution = 10
data = helperBackProjection(reflectedSig,rnggrid,fastTime,fc,fs,prf,radarVelocity,crossRangeResolution,c)


