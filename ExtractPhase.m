function [phaseExtractPoints,muFFT] = ExtractPhase(rxSig,m,targetIndex,numHuTgt,numRecv,reflectedSig)

if numRecv ~= 1
    % this is done because it is very important to extract phase in the
    % right FFT range bins. There can be some movement of the peak FFT
    % magnitude but averaging this helps.
    meanfft = fft(mean(reflectedSig(:,:,1:m),3)); %Take the mean of the FFT over 100 time samples
    
    
    targetIndex = sub2ind(size(meanfft),targetIndex(2,:),targetIndex(1,:)); %Convert target index to indcies instead of row,col
    fftAtTarget = abs(meanfft(targetIndex)); %The targetIndex is the
    %          actual range bin of the target. if TargetIndex ==
    %          phaseExtractPoints, then this will be the same as peakFFT.
    muFFT = mean(abs(meanfft)); %The average FFT over freq (Not important used for debugging)
    
    for fftcnt = 1:numHuTgt
        [pks,locs] = findpeaks(abs(meanfft(:,fftcnt)));
        [pks,sortIdx]=sort(pks); %Find and sort the peaks of the fft
        locs=locs(sortIdx); %apply the sorting to the range bins
        
        pks = flip(pks);
        locs = flip(locs); %Flip the sorting to be from greatest to least
        
        
        phaseExtractPoints(fftcnt) = locs(1); %Set the phase extract points to be the largest fft peaks
        peakFFT(fftcnt) =pks(1); %value of largest FFT peak, used for debugging
        
        
        
    end
    phaseExtractPoints = sub2ind(size(reflectedFFT),phaseExtractPoints',1:numHuTgt);
    %Convert target index from row,col to indcies
    %This is needed because with multiple targets the fft will be
    %size MxL where M is len of pulse and L is num of targets.
    
    
    
else
    
        % this is done because it is very important to extract phase in the
        % right FFT range bins. There can be some movement of the peak FFT
        % magnitude but averaging this helps.
        
        meanfft = fft(mean(reflectedSig(:,1:m),2)); %Take the mean of the FFT
        fftAtTarget = abs(meanfft(targetIndex)); %The targetIndex is the
%          actual range bin of the target. if TargetIndex ==
%          phaseExtractPoints, then this will be the same as peakFFT.
        muFFT = mean(abs(meanfft)); %The average FFT 
        [pks,locs] = findpeaks(abs(meanfft));
        [pks,sortIdx]=sort(pks); %Find and sort the peaks of the fft
        locs=locs(sortIdx); %apply the sorting to the range bins 

        pks = flip(pks);
        locs = flip(locs); %Flip the sorting to be from greatest to least


        phaseExtractPoints = locs(1:numHuTgt); %Set the phase extract points to be the largest fft peaks
        peakFFT =pks(1:numHuTgt); %value of largest FFT peak, used for debugging
    
    
    
end


end

