% [t,val]=plotATM('chF01m');


% save('chF01_datExtracted.mat')

load('chF01_datExtracted.mat')

%%


freqOut = []
for k = 1:length(RR)
    
    numSampHB = RR(k) * 250;
    freqHB = 1/RR(k);
    freqOut = [freqOut;freqHB*ones(numSampHB,1)];
    
end



numPlot = 20*250;
figure;
yyaxis left
plot(t(1:numPlot),val(1,1:numPlot))

hold on

yyaxis right

plot(t(1:numPlot),freqOut(1:numPlot));
ylim([1,1.5])



% save('Ml_EKGDat_fs250.mat','freqOut','val','t')


[s,f,t]  = spectrogram(val(1,:),5000,4000,10000,fs,'yaxis');

specLimited = s(1:75,:);
flimited = f(1:75);

specLimDb = 20*log(abs(specLimited));


noiseEst = mean(mean(specLimDb(60:end,:)));


imagesc(t,flimited,specLimDb)
title('Spectrogram of EKG')
xlabel('Time (s)')
ylabel('Frequency (Hz)')


threshSpec = (specLimDb>2.5*noiseEst);
imagesc(t,flimited,threshSpec)
title('Thresholded Spectrogram')
xlabel('Time (s)')
ylabel('Frequency (Hz)')

filtFreq = 1:size(threshSpec,1);
filtFreq = filtFreq>35;

filtFreq = filtFreq'.*ones(size(threshSpec));

threshSpecCutof = threshSpec .* logical(filtFreq);
imagesc(t,flimited,threshSpecCutof)
title('Thresholded Spectrogram edited')
xlabel('Time (s)')
ylabel('Frequency (Hz)')

dIn = threshSpecCutof;

save('ekgKalmanINdat.mat','dIn','flimited','t')



% Compare results after EKF ing
t_afterSpec = q;
q=freqOut(t_afterSpec*250)

d=(q-hrEst') %hrEst from EKF

mseEKF = mean(d(1:800).^2)
pe = mean(abs(q-hrEst')./q) %percent error.
