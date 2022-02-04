
% fc = 4e9;
c = physconst('LightSpeed');
rangeResolution = .1;
crossRangeResolution = .1;

bw = c/(2*rangeResolution);


fs = 150e6;
c = physconst('LightSpeed');
fc = 26e9;


prf = 10000; 
aperture = 4;  
tpd = 3*10^-6; 
fs = 120*10^6;



waveform = phased.LinearFMWaveform('SampleRate',fs, 'PulseWidth', tpd, 'PRF', prf,...
    'SweepBandwidth', bw);
sig = waveform();

speed = 1;
flightDuration = 4;
radarPlatform  = phased.Platform('InitialPosition', [0;0;0], 'Velocity', [0; speed; 0]);
slowTime = 1/prf;
numpulses = flightDuration/slowTime +1;
numpulses = round(numpulses);
maxRange = 50;
truncrangesamples = ceil((2*maxRange/c)*fs);
fastTime = (0:1/fs:(truncrangesamples-1)/fs);
% Set the reference range for the cross-range processing.
Rc = 1000;

antenna = phased.CosineAntennaElement('FrequencyRange', [1e9 6e9]);
% antennaGain = aperture2gain(aperture,c/fc);
antennaGain = 5;
transmitter = phased.Transmitter('PeakPower', 50e3, 'Gain', antennaGain);
radiator = phased.Radiator('Sensor', antenna,'OperatingFrequency', fc, 'PropagationSpeed', c);

collector = phased.Collector('Sensor', antenna, 'PropagationSpeed', c,'OperatingFrequency', fc);
receiver = phased.ReceiverPreamp('SampleRate', fs, 'NoiseFigure', 30);

channel = phased.FreeSpace('PropagationSpeed', c, 'OperatingFrequency', fc,'SampleRate', fs,...
    'TwoWayPropagation', true);
%% Targets
targetpos= [800,0,0;1000,0,0; 1300,0,0]';

pedest = backscatterPedestrian( 'Height',1.8, ...
    'OperatingFrequency',fc,'InitialPosition',[10;0;0], ...
    'InitialHeading',0,'WalkingSpeed',0.00001);
[targetpos,bpvel,bpax] = move(pedest,10,0);




targetvel = [0,0,0;0,0,0; 0,0,0]';
targetvel = zeros(3,16);

target = phased.RadarTarget('OperatingFrequency', fc, 'MeanRCS', ones(1,16));
pointTargets = phased.Platform('InitialPosition', targetpos,'Velocity',targetvel);
% The figure below describes the ground truth based on the target
% locations.
figure;h = axes;plot3(targetpos(1,:),targetpos(2,:),targetpos(3,:),'*b');
set(h,'Ydir','reverse');

% xlim([-10 10]);ylim([700 1500]);


title('Ground Truth');ylabel('Range');xlabel('Cross-Range');
%%






refangle = zeros(1,size(targetpos,2));
rxsig = zeros(truncrangesamples,numpulses);
for ii = 1:numpulses
   % Update radar platform and target position
    [radarpos, radarvel] = radarPlatform(slowTime);
    [targetpos,targetvel] = pointTargets(slowTime);
    
    % Get the range and angle to the point targets
    [targetRange, targetAngle] = rangeangle(targetpos, radarpos);
    
    % Generate the LFM pulse
    sig = waveform();
    % Use only the pulse length that will cover the targets.
    sig = sig(1:truncrangesamples);
    
    % Transmit the pulse
    sig = transmitter(sig);
    
    % Define no tilting of beam in azimuth direction
    targetAngle(1,:) = refangle;
    
    % Radiate the pulse towards the targets
    sig = radiator(sig, targetAngle);
    
    % Propagate the pulse to the point targets in free space
    sig = channel(sig, radarpos, targetpos, radarvel, targetvel);
    
    % Reflect the pulse off the targets
    sig = target(sig);
    
    % Collect the reflected pulses at the antenna
    sig = collector(sig, targetAngle);
    
    % Receive the signal  
    rxsig(:,ii) = receiver(sig);
    ii
end

imagesc(real(rxsig));title('SAR Raw Data')
xlabel('Cross-Range Samples')
ylabel('Range Samples')

pulseCompression = phased.RangeResponse('RangeMethod', 'Matched filter', 'PropagationSpeed', c, 'SampleRate', fs);
matchingCoeff = getMatchedFilter(waveform);
[cdata, rnggrid] = pulseCompression(rxsig, matchingCoeff);

imagesc(real(cdata));title('SAR Range Compressed Data')
xlabel('Cross-Range Samples')
ylabel('Range Samples')


Rc = 10
rma_processed = helperRangeMigration(cdata,fastTime,fc,fs,prf,speed,numpulses,c,Rc);


bpa_processed = helperBackProjection(cdata,rnggrid,fastTime,fc,fs,prf,speed,crossRangeResolution,c);

figure(1);
imagesc((abs((rma_processed(1700:2300,600:1400).'))));
title('SAR Data focused using Range Migration algorithm ')
xlabel('Cross-Range Samples')
ylabel('Range Samples')





