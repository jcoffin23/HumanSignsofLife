function [fr,fh] = estFreqs(ct,fs)

%This function is used to estimate the heart rate and respitory rate of a
%given input chest compression signal.


% The following lines of code are used to determine the effective
% ness of the filter.

%This is the FFT size to take. We need a large value to get a good
%enough picture of the freqeuncies (1024 will probably work)
Npt = 2^nextpow2(length(ct));

%       Get the second order filter defined in the function getRespFilter
%       It is a butterworth filter that has a freqeucny response that will
%       extract only the respirtory information from the signals.
Hd = getRespFilter(fs);
sos = Hd.sos;

%Filter the Signal ct
%Then normalize it since the filter has a very high gain
ctResp = sosfilt(sos,ct);
ctResp = ctResp/max(ctResp);

%       Get the second order filter defined in the function getHeartFilter
%       It is a butterworth filter that has a frequency response that will
%       extract only the heart beat information from the signals.
Hd = getHeartFilter(fs);
sos = Hd.sos;

%Filter the noisey Refected signal ct to have only the heart beat
%in it
%Then normalize it since the filter has a very high gain
ctHeart = sosfilt(sos,ct);
ctHeart = ctHeart/max(ctHeart);

%Take the FFT of the resp rate filtered signals
[fftResp,f] =easyFFT(ctResp,fs,Npt);
[fftHeart,~] = easyFFT(ctHeart,fs,Npt);

[fftTotal,f] = easyFFT(ct,fs,Npt);



[fftTotal,f] = easyFFT(ct,fs,Npt);
genFig = 0
if genFig
    figure
    fftTotal=abs(fftTotal);
    fftTotal=normalize(fftTotal,'range');
    plot(f,abs(fftTotal))
    title('Total Chest Signal FFT')
    xlabel('Freq (Hz)')
    ylabel('Mag')
    xlim([-5,5])
    
    
    figure
    fftHeart1=abs(fftHeart);
    fftHeart1=normalize(fftHeart1,'range');
    plot(f,abs(fftHeart1))
    title('Chest Signal FFT (Post Heart Filter) ')
    xlabel('Freq (Hz)')
    ylabel('Mag')
    xlim([-5,5])
    
    
end


%Find the index of the highest power part of the FFT
% This should be the to the respitory signal only
% then find the freqeuncy that correspondes to it it
% This will be the estimated respitory rate in Hz.
% This will get done for the raw reflected signal and the LMS
% filtered signal

%         dy = diff(fftTotal);
%         dx=diff(f);
%         der = abs(dy./dx); %Find the Deritive of the FFT
%
%         [pks,locs,w,p] = findpeaks(abs(fftTotal)); %Find all the peaks
%
%         [p,pSortIdx]=sort(p,'descend'); %sort the prominate peaks
%         pks = pks(pSortIdx); %Sort the pks by how prominate they are
%         locs = locs(pSortIdx);%Sort the locs indexs to match
%

%         fr = abs(f(locs(1))); %Take the most prominate peak as the respiration rate
[~,maxIdx] = max(abs(fftResp));
fr = abs(f(maxIdx));

[~,maxIdx] = max(abs(fftHeart));
fh = abs(f(maxIdx));

%         fh = abs(f(locs(4))); %Take the second most prominate peak as the Heart rate
end


%{
Ignore this pleae
  
        
        %Use the find peaks tool to find all of the peaks
        [pks,locs] = findpeaks(abs(fftTotal));
        
        %Sort by by the biggest peaks
        [pks,sortIDx]=sort(pks);
        
        %Remove any duplicates (force single side band)
        pks = unique(pks);
        
        
        fSort = locs(sortIDx);
        fSort=f(fSort);
        fSort = fSort(fSort>0);
        fSort = fSort(end-5:end);
        pks = pks(end-5:end);
        
        pks = pks(fSort>.1);
        
        %Fr is the highest amplitude FFT value so the freqeuncy value that
        %coresponds to the highest peak should be the Fr result.
        fr = f((abs(fftTotal) == pks(end)));
        fr = abs(fr(end));


        %fh should be the second highest FFT value, so the freq value that
        %coresponds to the second highest peak should be the Fh result.
        fh = f((abs(fftTotal) == pks(end-1)));
        fh = abs(fh(end));
%
%         locs=locs(sortIDx);

%}
