fs = 20000;
fsChest = fs;
Tmax = 20;
heartHeight = .5*10^-3; 
respHeight = 10 * 10^-3;
t = 0:1/fs:Tmax;


numDat = 100000;%20k data samples (training and testing and Validation);

a = .1;
b = 2;
x1 = a + (-a+b)*rand(numDat,1); %Actual Frequency data

datOut = nan(length(t),numDat);

a = 0;
b = 2*pi;
phi =  a + (-a+b)*rand(numDat,1); %Actual Frequency data
phi = zeros(size(x1));
for k = 1:numDat
    
    datOut(:,k) = sin(2*pi*x1(k)*t + phi(k))' + .01*randn(length(t),1);
    
    
end


save('ML_singleFreq2.mat','datOut','x1','phi','-v7.3')