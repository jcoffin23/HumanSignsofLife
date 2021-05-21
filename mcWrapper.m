
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



%% 
varVector = linspace(1,.001,20); %This is the vector that contains the different parameters for the tests. THis will be the xaxis in the finial curve of Error vs (Parameter)


numVar = length(varVector); %Find out how long it is
numMCTrials = 10; %How many MC trials to average over per varVec value. Total number of runs is numMCtrials * numVar
respRoc = zeros(1,numVar);%This is the finial ROC curve for respitory error
heartRoc = zeros(1,numVar);%Finial ROC curve for heart rate error

%% Set up Defaults 


fc = 5e9;
bw = .18*10^9;
rxRad = 0.5; 
heartRatio = .3; %Heart rate is 30 percent of the respitory rate. 
respHeight = 1;
%% Variable Loop
% This outer loop will change the variable being tested over. 


for k = 1:numVar

RespErrVec = zeros(1,numMCTrials); %Error Vector to be averaged over
HeartErrVec = zeros(1,numMCTrials);
tic

%     bw = varVector(k)*10^9;
%     fc=varVector(k)*1e9;
    respHeight = varVector(k); 

%% Loop same conditions over 10 trials
% Loop over the same conditions 10 times and store the error results. 
for trialNum = 1:numMCTrials
   

   
    [RespErrVec(trialNum),HeartErrVec(trialNum),~] = singleMCTrial(fc,bw,rxRad,respHeight,heartRatio);
    
    
end
toc
heartRoc(k) = mean(HeartErrVec);%Average over the error to get the error estimate for this variable value 
respRoc(k) = mean(RespErrVec);

end

%% Plotting
%Plot the two curves, the labels and title need to be changed to match the
%variable that was changed. 

plot(varVector,heartRoc)
title('Heart Rate Error vs Amplitude of respiratory signal')
ylabel('Error (Hz)')
xlabel('Amplitude of respiratory signal (m)')

plot(varVector,respRoc)
title('Respiratory Rate Error vs Amplitude of respiratory signal')
ylabel('Error (Hz)')
xlabel('Amplitude of respiratory signal (m)')