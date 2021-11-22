
fs = 20000;
fsChest = fs;
Tmax = 2;
heartHeight = .5*10^-3; 
respHeight = 10 * 10^-3;
t = 0:1/fs:Tmax;




numRuns = 10000;

frActual = nan(numRuns,1)';
fhActual = nan(numRuns,1)';
datOut = nan(length(t),numRuns); 
for k = 1:numRuns
    
    
    [~,datOut(:,k),~,frActual(k),fhActual(k)] = getChestCompression(fsChest,[],length(t),0,0,0,respHeight,heartHeight,t);
k

end


% save('C:\Users\Joe\Desktop\ML_genTrain_fs20k_1.mat','datOut','frActual','fhActual','t','-v7.3')