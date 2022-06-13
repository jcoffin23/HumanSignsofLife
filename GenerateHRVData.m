
fs = 500;
fsChest = fs;
Tmax = 50;
heartHeight = .5*10^-3; %Heart rate is 30 percent of the respitory rate.
respHeight = 10 * 10^-3;
t = 0:1/fs:Tmax;




numRuns = 10;

frActual = nan(numRuns,1)';
fhActual = nan(length(t),numRuns);
datOut = nan(length(t),numRuns); 
for k = 1:numRuns
    
    
    [~,datOut(:,k),~,frActual(k),fhActual(:,k)] = hrvChestComp(fsChest,[],length(t),0,0,0,respHeight,heartHeight,Tmax);
k

end


dOut = datOut(:,k);
fhA = fhActual(:,k);
frA = frActual(k);

 save('EKFDat_HRV_long.mat','fhA','frA','dOut','fs')
% save('C:\Users\Joe\Desktop\ML_genTrain_HRV_full_fs500_long_1.mat','datOut','frActual','fhActual','t','-v7.3')

save('C:\Users\Joe\Desktop\ML_HRV_full_fs500_long_1.mat','datOut','frActual','fhActual','t','-v7.3')
% 
% 
% 
% % q=datOut(:,1);
% % % 
% % qq = fftshift(fft(q(1:500)));
% % % 
% figure
%  spectrogram(q,2048,2000,10000,fs,'yaxis');
% ylim([0,3.5])
% hold on 
% plot(t,fhActual(:,1))

%% One Possible solution is to use built in Matlab Stuff
[~,f,t,p]  = spectrogram(q,2048,2000,10000,fs,'yaxis');

%f will have (nfft/2) +1 elements
[fridge,~,lr] = tfridge(p,f);

% hold on
% plot3(t,fridge,abs(p(lr)),'LineWidth',4)
plot(t,x(round(t*fs)))
hold off

%%
% STFT
stft_options.sr = 500;
stft_options.M1 = 512;
stft_options.M2 = 2^12;
stft_options.hop = 32;
stft_options.win = 'hann';
stft_options.N = 2^14;

% PLCA
plca_options.max_iter = 10;
plca_options.dic_size = 3;
plca_options.seg_size = stft_options.sr;


[super_spec, spec_t, spec_f, w, zf, zt, h] = super_spectrogram(q, stft_options, plca_options);
t = (0:size(spec_t,2)-1)*stft_options.hop/stft_options.sr;
f = stft_options.sr*[0:size(spec_t,1)-1]/stft_options.N;

max_val = max(max(db(abs(super_spec))));
figure
imagesc(t, f, db(abs(super_spec)),[max_val-60 max_val]); axis xy; %colorbar;
title('Super-resolution spectrogram','fontsize',12);
xlabel('Time [sec.]','fontsize',12); ylabel('Frequency [kHz]','fontsize',12);
set(gca,'fontsize',12,'fontname','times');

ylim([0,4])


%% 

%  save('C:\Users\Joe\Desktop\ML_HRV_honly_fs20k_1.mat','datOut','frActual','fhActual','-v7.3')

% save('C:\Users\Joe\Desktop\ML_genTrain_HRV_honly_fs500_1.mat','datOut','frActual','fhActual','t','-v7.3')

%  save('C:\Users\Joe\Desktop\ML_HRV_honly_fs500_1.mat','datOut','frActual','fhActual','t','-v7.3')