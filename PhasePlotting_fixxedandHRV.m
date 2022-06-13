
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


tgtStruct.human = [1];
tgtStruct.offset = [15];


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
inputStruct.gain = 18;
inputStruct.collectAllTargetPhase = 0;
inputStruct.miliPow = 20;


inputStruct.bw = bw;

hrVector = [0,1];

lenTest = length(hrVector);
huPhase = zeros(numSamps,lenTest);
heartActual= zeros(numSamps,lenTest);
clutPhase = zeros(numSamps,lenTest);
for k = 1:lenTest
   tgtStruct.legacyHR = hrVector(k)
[frError,fhError,rangeRes,chestSig,peakFFT,tgtStructOut,recvPow] = singleMCTrial(inputStruct,tgtStruct);

heartActual(:,k) = tgtStructOut.fhActual;
huPhase(:,k) = chestSig(1,:);




    
end
heartActual(:,2) = heartActual(1,2)*ones(numSamps,1);

figure
plot(tgtStructOut.t,huPhase)
xlabel('Time(s)')
ylabel('c(t)')
title('Chest Compression vs time')
legend('With HRV simulation','Without HRV simulation')
saveas(gcf,[pathS,'\chestcompression.png'])
figure
plot(tgtStructOut.t,heartActual)
xlabel('Time(s)')
ylabel('Heart Rate (Hz)')
title('Heart Rate comparison')
legend('With HRV simulation','Without HRV simulation')
saveas(gcf,[pathS,'\compressionHeartRate.png'])

pathS='C:\Users\Joe\OneDrive - George Mason University - O365 Production\RTL\Pictures\Journal'

save('pt2.mat')
