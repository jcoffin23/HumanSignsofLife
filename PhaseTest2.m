
%% This Wrapper sets up and runs 10 simulations under the same conditons and then changes the conditions. This allows us to get a curve that
% shows how the system reacts to different conditons.

%The idea is to run several simulations (10) and average the error results
%to get a good estimate of the error at a set of conditions.

%So Far I am only testin for Bandwidth, Carrer Freq, Array Width, the
%amplitude of the chest compression, and the ratio of the chest compression
%ampltidue that is the heart beat.

% Defaults:
%
% Carrier Freq (fc) = 5Ghz
% Bandwidth(bw) = 180 MHz
% Rx Array width (rxRad) = 0.5
% Resp Amplitude (respHeight) = 1
% Heart Rate Ratio(heartRatio) = 0.3
c = 3e8;
% fc = 5e9;
% lambda = c/fc;

feet2meter  = 0.3048;% 1 foot is this many meters.

%% Set Up Variables to be used as test inputs
% Preallocate the arrays for results


numSamps = 4000; %this is known because of internal decisions. bad code?




%% Set up Defaults
fc = 24e9;

lambda = c/fc;
bw = .18e9;
bw = .5e9;
beatError = 0;

%Respitory Height is about 4 - 12 mm
% Heart heightis about .2-.5mm
heartHeight = .5*10^-3; %Heart rate is 30 percent of the respitory rate.
respHeight = 10 * 10^-3;
clutterRCS = .5; % Mean RCS of each background target
stv = 20;
addNoise = 0; % Should noise be added




%% Recieve array setup

rxRad = .5*feet2meter;

NumeleMents = 1;

recvArray=phased.CrossedDipoleAntennaElement;
usingArray = 0; %Set to 0 because array is element not actual array



%% Send Array Setup

sendElmnt = phased.ShortDipoleAntennaElement('AxisDirection','Y');

%% Target Def


tgtStruct.human = [1,0];
tgtStruct.offset = [15,5];


numTgt = length(tgtStruct.human);


tgtStruct.yoffset= zeros(1,numTgt);
tgtStruct.zoffset = zeros(1,numTgt);
tgtStruct.fhActual=zeros(numSamps,numTgt);
tgtStruct.frActual=zeros(1,numTgt);

tgtStruct.numTgt = numTgt;
meanHumanRCS = 0.1992;
tgtStruct.scatterMatt = [2*meanHumanRCS meanHumanRCS;meanHumanRCS 2*meanHumanRCS];
%% For Loop



%% Loop same conditions over 10 trials
% Loop over the same conditions 10 times and store the error results.

%         tic
inputStruct.fc = fc;
inputStruct.bw = bw;
inputStruct.recvArray = recvArray;
inputStruct.usingArray = usingArray;
inputStruct.sendElmnt = sendElmnt;
inputStruct.respHeight = respHeight;
inputStruct.heartHeight = heartHeight;
inputStruct.beatError = beatError;
inputStruct.clutterRCS = clutterRCS;
inputStruct.stv = stv;
inputStruct.addNoise = addNoise;
inputStruct.numRecv = (NumeleMents * usingArray) + (1*(1-usingArray)); % If using array, then numelements will be equal to number of elemnnts, otherwise it will be 1
inputStruct.enablePic = 0;
inputStruct.rxRad=rxRad;
inputStruct.bistatic = 0;


inputStruct.miliPow = 100;


inputStruct.bw = bw;


clutterOffsetVec = [10:5:30];
lenTest = length(clutterOffsetVec);
huPhase = zeros(numSamps,lenTest);

clutPhase = zeros(numSamps,lenTest);
for k = 1:lenTest
   
    %Paramters that are different for this test are below
tgtStruct.human = [1,0];
tgtStruct.offset = [5,clutterOffsetVec(k)];
inputStruct.collectAllTargetPhase = 1;

[frError,fhError,rangeRes,chestSig,peakFFT,tgtStructOut,recvPow] = singleMCTrial(inputStruct,tgtStruct);


huPhase(:,k) = chestSig(1,:);

clutPhase(:,k) = chestSig(2,:);


    
end

figure
plot(tgtStructOut.t,huPhase)
xlabel('Time(s)')
ylabel('c(t)')
title('Chest Compression vs time')
% legend('Clutter: 15m','Clutter: 20m','Clutter: 25m','Clutter: 30m')

figure
plot(tgtStructOut.t,clutPhase)
xlabel('Time(s)')
ylabel('Clutter Phase')
title('Clutter Phase vs time')
% legend('Clutter: 15m','Clutter: 20m','Clutter: 25m','Clutter: 30m')

save('pt2.mat')
