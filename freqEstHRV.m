
fs=500;
lenn = 1500;
numTest = 10;
kVec = 1:10:length(datOut)-lenn;
fr = zeros(length(kVec),numTest);
fh = zeros(length(kVec),numTest);

MSEMat = zeros(numTest,2);
for testNum = 1:numTest
    tic
    cc=1;
    for k = kVec
        ct = datOut(k:k+lenn);
        [fr(cc,testNum),fh(cc,testNum)] = estFreqs(ct- mean(ct),fs);
        cc=cc+1;
    end
    toc
    
    MSEMat(testNum,1) = mean(  (fh(:,testNum) - frActual(testNum)).^2 );
    MSEMat(testNum,2) = mean((fr(:,testNum) - fhActual(kVec,testNum)).^2);
    
    testNum
    
end

