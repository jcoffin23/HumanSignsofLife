
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

% varVectory = [[-5,-2,-1,1,2,5] , zeros(1,6)];
%
% varVectorx = [10*ones(1,6),[15,18,19,21,22,25]];


% offSetVec = linspace(5,100,41);


% powerVec = linspace(1,20,40) * 1e9;



offSetVec = linspace(5,49,40);


powerVec = linspace(1,77,50) * 1e9;





a= 5;
b=15;
% offSetVec =  a + (b-a).*rand(100,1);


% powerVec = 5;
numVar = length(offSetVec); %Find out how long it is
numPow = length(powerVec)




numMCTrials = 5; %How many MC trials to average over per varVec value. Total number of runs is numMCtrials * numVar
respRoc = zeros(1,numVar);%This is the finial ROC curve for respitory error
heartRoc = zeros(1,numVar);%Finial ROC curve for heart rate error


numSamps = 4000
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
%
%
%
%
% directionalInPut = 10;
% antElmnt = phased.CosineAntennaElement('CosinePower',[1,directionalInPut]);
%
%
% recvArray = phased.ULA('Element',antElmnt,'NumElements',NumeleMents,...
%     'ElementSpacing',rxRad);

recvArray = phased.ShortDipoleAntennaElement('AxisDirection','Y');
recvArray=phased.CrossedDipoleAntennaElement;
usingArray = 0; %Set to 0 because array is element not actual array
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


%% Send Array Setup
% sendElmnt = phased.IsotropicAntennaElement('BackBaffled',true);
% sendElmnt = phased.CosineAntennaElement('CosinePower',[1,directionalInPut]);
sendElmnt = phased.ShortDipoleAntennaElement('AxisDirection','Y');

%% Target Def
numTgt = 1;


tgtStruct.human = [1,1,1];
tgtStruct.offset = [7,10,14];
tgtStruct.yoffset= [0,0,0];
tgtStruct.zoffset = [0,0,0]
tgtStruct.fhActual=zeros(numSamps,numTgt);
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
    
    RespErrVec = zeros(1,numMCTrials); %Error Vector to be averaged over
    HeartErrVec = zeros(1,numMCTrials);
%     tic
    
    %% Loop same conditions over 10 trials
    % Loop over the same conditions 10 times and store the error results.
    
    for trialNum = 1:length(powerVec)
%         tic
        inputStruct.fc = powerVec(trialNum);
        inputStruct.bw = bw;
        inputStruct.recvArray = recvArray;
        inputStruct.usingArray = usingArray;
        inputStruct.sendElmnt = sendElmnt;
        inputStruct.respHeight = respHeight;
        inputStruct.heartHeight = heartHeight;
        inputStruct.beatError = beatError;
        inputStruct.clutterRCS = clutterRCS;
        inputStruct.stv = stv;
        % inputStruct.targetOffset = targetOffset;
        inputStruct.addNoise = addNoise;
        inputStruct.numRecv = (NumeleMents * usingArray) + (1*(1-usingArray)); % If using array, then numelements will be equal to number of elemnnts, otherwise it will be 1
        inputStruct.enablePic = 0;
        inputStruct.rxRad=rxRad;
        
        inputStruct.bistatic = 0;
        
%         tgtStruct.offset =offSetVec(k);
%         tgtStruct.offset = [20];
%         inputStruct.miliPow = powerVec(trialNum);
         inputStruct.miliPow = 20
        %     if trialNum ==1
        %         inputStruct.enablePic = 1;
        %     end
        
        avgNum=1;
        avgResp = zeros(1,avgNum);
        avgHeart = zeros(1,avgNum);
        avgfftAtTarget = zeros(1,avgNum);
        avgmuFFT = zeros(1,avgNum);
        avgpeakFFT = zeros(1,avgNum);
        for avg = 1:avgNum
            inputStruct.bw = bw;
            [respErr,heartErr,~,chestSig(trialNum,k,:),fftAtTarget,muFFT,peakFFT,tgtStructOut,recvPow] = singleMCTrial(inputStruct,tgtStruct);
            
            FHA(trialNum,k,:) = tgtStructOut.fhActual;
            
            RPMat(trialNum,k) = recvPow; %Average recieve power for the test conditions. 
            
            avgResp(avg) = respErr;
            avgHeart(avg) =heartErr;
            avgfftAtTarget(avg) =fftAtTarget;
            avgmuFFT(avg) =muFFT;
            avgpeakFFT(avg) =peakFFT;
        end
        
        peakFFTMat(trialNum,k) = mean(avgpeakFFT);
        stdMat(trialNum,k) =mean( avgmuFFT);
        targetFFTMat(trialNum,k) = mean(avgfftAtTarget);
        respMat(trialNum,k) =mean(avgResp);
        heartMat(trialNum,k) = mean(avgHeart);
%         heartActualMat(trialNum,k) = tgtStructOut.fhActual;
         respActualMat(trialNum,k) = tgtStructOut.frActual;
%         toc
    end
%     toc
    heartRoc(k) = mean(HeartErrVec);%Average over the error to get the error estimate for this variable value
    respRoc(k) = mean(RespErrVec);
    
end




%  save('rawWorkspaceHRVPHase1.mat')


% datOut=jct.util.threedto2d(chestSig(:,:,:)).';
% heartActualMat  = heartActualMat(1:k-1);
% respActualMat = respActualMat(1:k-1); 

% offSetVec = offSetVec(1:k-1);

% save('ML_phaseFull_HRV','datOut','heartActualMat','respActualMat','offSetVec','-v7.3')
   save('tpdatBWGain20.mat') 
   datOut = chestSig;
 save('C:\Users\Joe\Desktop\ML_phaseFull_2dROC_fcOffset_20mw_fixbgain','datOut','FHA','respActualMat','offSetVec','-v7.3')

% preAmble = '2dROC2\2dROC_bistaic_updatedSMat\';
% 
% percentH = abs(heartMat)./heartActualMat;
% percentR = abs(respMat)./respActualMat;

%% Resmat Plotting




offSetVec ;


powerVec ;


x = offSetVec;
y = powerVec/1e9
figure
imagesc(x,y,respMat)
title('2D - ROC Resp Error')
xlabel('Offset')
ylabel('Carrier Freq')

% jct.util.SavetoPng([preAmble, 'RespError.png'])
figure
imagesc(percentR)
title('2D - ROC Resp Error')
xlabel('Offset')
ylabel('Tx Power')

figure
imagesc(heartMat)
title('2D - ROC Heart Error')
xlabel('Offset')
ylabel('Tx Power')

% jct.util.SavetoPng([preAmble, 'HeartError.png'])

figure
imagesc(20*log10(peakFFTMat))
title(['Value of Peak FFT (dB)'])
xlabel('Offset (m)')
ylabel('Tx Power (mW)')
% jct.util.SavetoPng([preAmble, 'peakFFT.png'])


figure
imagesc(20*log10(targetFFTMat))
title(['Value of FFT at Target Location (dB)'])
xlabel('Offset (m)')
ylabel('Tx Power (mW)')

% jct.util.SavetoPng([preAmble, 'targetFFT.png'])
figure

confFact = 20;
imagesc(~(targetFFTMat>confFact*stdMat))
title(['FFT value of Target above ',num2str(confFact),' times the average FFT'])
xlabel('Offset (m)')
ylabel('Tx Power (mW)')

% jct.util.SavetoPng([preAmble, 'targetThreshold.png'])
figure

confFact = 20;
imagesc(~(peakFFTMat>confFact*stdMat))
title(['FFT value of Target above ',num2str(confFact),' times the average FFT'])
xlabel('Offset (m)')
ylabel('Tx Power (mW)')


% jct.util.SavetoPng([preAmble, 'PEAKThreshold.png'])
%% Plotting
% %Plot the two curves, the labels and title need to be changed to match the
% %variable that was changed.
% figure
% subplot 211
% plot(varVector,heartRoc)
% title('Heart Rate Error vs Clutter Offset')
% ylabel('Error (Hz)')
% xlabel('Clutter Offset (m)')
% % ylim([0,1])
% 
% 
% subplot 212
% 
% plot(varVector,respRoc)
% title('Respiratory Rate Error vs Clutter Offset')
% ylabel('Error (Hz)')
% % ylim([0,.015])
% xlabel('Clutter Offset (m)')