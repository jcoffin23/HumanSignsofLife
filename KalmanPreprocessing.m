% [t,val]=plotATM('chF01m');


% save('chF01_datExtracted.mat')

% load('chF01_datExtracted.mat')
% load('Sim_HrvDat.mat')
load('Sim_NO_HrvDat.mat')
specLimDb = specNohrv';

t = tNohrv;
f = fNohrv;
%%




% [s,f,t]  = spectrogram(val(1,:),5000,4000,10000,fs,'yaxis');

% specLimited = s(1:75,:);
% flimited = f(1:75);

% specLimDb = 20*log(abs(specLimited));


% imagesc(tHRV,fhrv,specLimDb)

specLimDb = specLimDb + abs(min(min(specLimDb)));
noiseEst = mean(mean(specLimDb(100:end,:)));

imagesc(t,f,specLimDb)
title('Spectrogram of Non')
xlabel('Time (s)')
ylabel('Frequency (Hz)')


threshSpec = (specLimDb>1.5*noiseEst);

figure
imagesc(t,f,threshSpec)
title('Thresholded Spectrogram')
xlabel('Time (s)')
ylabel('Frequency (Hz)')

filtFreq = 1:size(threshSpec,1);
filtFreq = filtFreq>15;

filtFreq = filtFreq'.*ones(size(threshSpec));

threshSpecCutof = threshSpec .* logical(filtFreq);
imagesc(t,f,threshSpecCutof)
title('Thresholded Spectrogram edited')
xlabel('Time (s)')
ylabel('Frequency (Hz)')

dIn = threshSpecCutof;

% save('nOhrvKalDat.mat','dIn','f','t')



% Compare results after EKF ing

d=(fNohrvtrue-hrEst') %hrEst from EKF

mseEKF = mean(d(1:800).^2)
pe = mean(abs(fNohrvtrue-hrEst')./fNohrvtrue) %percent error.
