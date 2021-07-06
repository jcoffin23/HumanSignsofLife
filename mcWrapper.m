
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
fc = 5e9;
lambda = c/fc;


%% 
% varVector = linspace(.1,.0001,20); %This is the vector that contains the different parameters for the tests. THis will be the xaxis in the finial curve of Error vs (Parameter)
% varVector = logspace(-1,-6,20);

% varVector = lambda/2
varVector = 20:-1:10;
% varVector = varVector(2:end);

% varVector = -20:1:20
% varVector = linspace(.1e9,2e9,20);
numVar = length(varVector); %Find out how long it is
numMCTrials = 5; %How many MC trials to average over per varVec value. Total number of runs is numMCtrials * numVar
respRoc = zeros(1,numVar);%This is the finial ROC curve for respitory error
heartRoc = zeros(1,numVar);%Finial ROC curve for heart rate error


%% Set up Defaults 


fc = 5e9;
bw = .18*10^9;

% bw = .5e9;
rxRad = 0.5; 
beatError = 0;

%Respitory Height is about 4 - 12 mm 
% Heart heightis about .2-.5mm
heartRatio = .5*10^-3; %Heart rate is 30 percent of the respitory rate. 
respHeight = 10 * 10^-3;
clutterRCS = .5;
%% Variable Loop
% This outer loop will change the variable being tested over. 


for k = 1:numVar
k
RespErrVec = zeros(1,numMCTrials); %Error Vector to be averaged over
HeartErrVec = zeros(1,numMCTrials);
tic

%     bw = varVector(k);
%     fc=varVector(k)*1e9;
%     respHeight = varVector(k); 
% rxRad = varVector(k);
% heartRatio = varVector(k);
% clutterRCS = varVector(k);
stv = varVector(k);
% beatError = varVector(k);
% clutterRCS = 0;
%% Loop same conditions over 10 trials
% Loop over the same conditions 10 times and store the error results. 
for trialNum = 1:numMCTrials
   

   
    [RespErrVec(trialNum),HeartErrVec(trialNum),~] = singleMCTrial(fc,bw,rxRad,respHeight,heartRatio,beatError,clutterRCS,stv);
    
    
end
toc
heartRoc(k) = mean(HeartErrVec);%Average over the error to get the error estimate for this variable value 
respRoc(k) = mean(RespErrVec);

end

%% Plotting
%Plot the two curves, the labels and title need to be changed to match the
%variable that was changed. 
plot(varVector,heartRoc)
figure
title('Heart Rate Error vs SNR')
ylabel('Error (Hz)')
xlabel('SNR (dB)')
% ylim([0,1])


figure


plot(varVector,respRoc)
title('Respiratory Rate Error vs SNR')
ylabel('Error (Hz)')
% ylim([0,.015])
xlabel('SNR (dB)')