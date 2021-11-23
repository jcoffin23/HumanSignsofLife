
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


offSetVec = linspace(5,49,40);


powerVec = linspace(1,20,50) ;


numVar = length(offSetVec); %Find out how long it is
numPow = length(powerVec);




avgNum = 5; %How many MC trials to average over per varVec value. Total number of runs is numMCtrials * numVar
respRoc = zeros(1,numVar);%This is the finial ROC curve for respitory error
heartRoc = zeros(1,numVar);%Finial ROC curve for heart rate error


numSamps = 4000;
FHA = zeros(numPow,numVar,numSamps);
chestSig = zeros(numPow,numVar,numSamps);


%% Set up Defaults
fc = 5e9;
% bw = 2*10^9;
lambda = c/fc;
bw = .5e9;
% rxRad = 0.5;
beatError = 0;

%Respitory Height is about 4 - 12 mm
% Heart heightis about .2-.5mm
heartHeight = .5*10^-3; %Heart rate is 30 percent of the respitory rate.
respHeight = 10 * 10^-3;
clutterRCS = .5; % Mean RCS of each background target
stv = 20;
addNoise = 0; % Should noise be added




%% Recieve array setup
% rxRad = lambda/2;

rxRad = .5*feet2meter;

NumeleMents = 1;

recvArray=phased.CrossedDipoleAntennaElement;
usingArray = 0; %Set to 0 because array is element not actual array
%



%% Send Array Setup

sendElmnt = phased.ShortDipoleAntennaElement('AxisDirection','Y');

%% Target Def
numTgt = 1;

tgtStruct.legacyHR=1; %Set to use constant heart rate
tgtStruct.human = [1];
tgtStruct.offset = [10];
tgtStruct.yoffset= [0];
tgtStruct.zoffset = [0];
tgtStruct.fhActual=zeros(1,numTgt);
tgtStruct.frActual=zeros(1,numTgt);

tgtStruct.numTgt = length(tgtStruct.human);
meanHumanRCS = 0.1992;
tgtStruct.scatterMatt = [2*meanHumanRCS meanHumanRCS;meanHumanRCS 2*meanHumanRCS];
% tgtStruct.scatterMatt = cat(3,tgtStruct.scatterMatt,tgtStruct.scatterMatt)
%% For Loop
peakFFTMat = zeros(length(powerVec),length(offSetVec));
stdMat = zeros(size(peakFFTMat));
respMat=zeros(size(peakFFTMat));
heartMat=zeros(size(peakFFTMat));
targetFFTMat=zeros(size(peakFFTMat));
heartActualMat=zeros(size(peakFFTMat));
respActualMat=zeros(size(peakFFTMat));


for k = 1:numVar
    
%     RespErrVec = zeros(1,numMCTrials); %Error Vector to be averaged over
%     HeartErrVec = zeros(1,numMCTrials);
    %     tic
    
    %% Loop same conditions over 10 trials
    % Loop over the same conditions 10 times and store the error results.
    
    for trialNum = 1:length(powerVec)
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
        
         inputStruct.collectAllTargetPhase = 0;
        inputStruct.miliPow = powerVec(k);
        
        avgResp = zeros(1,avgNum);
        avgHeart = zeros(1,avgNum);
        avgfftAtTarget = zeros(1,avgNum);
        avgmuFFT = zeros(1,avgNum);
        avgpeakFFT = zeros(1,avgNum);
        for avg = 1:avgNum
            inputStruct.bw = bw;
            [respErr,heartErr,~,chestSig(trialNum,k,:),peakFFT,tgtStructOut,recvPow] = singleMCTrial(inputStruct,tgtStruct);
            
            avgResp(avg) = respErr;
            avgHeart(avg) =heartErr;
            avgpeakFFT(avg) =peakFFT;
            avgmuFFT(avg) = tgtStructOut.muFFT;
        end
        
        peakFFTMat(trialNum,k) = mean(avgpeakFFT);
        stdMat(trialNum,k) =mean( avgmuFFT);
        targetFFTMat(trialNum,k) = mean(avgfftAtTarget);
        respMat(trialNum,k) =mean(avgResp);
        heartMat(trialNum,k) = mean(avgHeart);
        heartActualMat(trialNum,k) = tgtStructOut.fhActual;
        respActualMat(trialNum,k) = tgtStructOut.frActual;
        %         toc
    end
    %     toc
    %     heartRoc(k) = mean(HeartErrVec);%Average over the error to get the error estimate for this variable value
    %     respRoc(k) = mean(RespErrVec);
    
end

imagesc(powerVec,offSetVec,respMat)
xlabel('Offset(m)')
ylabel('TxPower(mW)')
title('Respiratory error vs Tx Power and Offset')

imagesc(powerVec,offSetVec,heartMat)
xlabel('Offset(m)')
ylabel('TxPower(mW)')
title('Heart rate error vs Tx Power and Offset')

