

numAvg = 5; %number of trials for each noise power
numNoise = 100; %Number of different noise powers

%initalize the error Matrices to be empty
nre = nan(numNoise,numAvg);
nhe= nan(numNoise,numAvg);
fre= nan(numNoise,numAvg);
fhe= nan(numNoise,numAvg);

sigPow = linspace(.000001,2,numNoise); %Set up noise power vector


for nl = 1:numNoise %loop over the noise vector
    
    
    
    for aIDX = 1:numAvg %for each noise vector run this inner part numAvg times
        
        a = 15/60;
        b = 18/60;
        fr = a + (b-a).*rand;
        %Set up fr (respiration rate) to be a random number uniformly generated between the bounds
        % a = 15/60 to b = 18/60. This is the range 15 to 18 breaths per minute
        % converted to Hz
        
        
        
        
        a = 1;
        b = 100/60;
        fh = a + (b-a).*rand;
        % set up fh (heart rate) to be a random number unifromly generated between
        % a = 1 and b = 100/60. This sets the heart beat ragnge from 60 to 100 bpm
        % (1 to 100/60) Hz. This represents a wide range of standing heart beats
        
        
        ar = 1;
        ah = .3;
        %ar is the amplitude of the respitory siganl
        %ah is the ampllitude of the heart beat signal.
        % ar and ah are chosen such that ah is much lower than ar
        
        
        fs = 100;
        t = 0:1/fs:100;
        %Set up the sampling rate and the time vector. I choose fs as 100 since we
        %are sampling a heart beat which is much less than 100Hz (ideally). Fs
        %could probably be lowered in practice.
        
        ct = ar*sin(2*pi*fr*t) + ah*sin(2*pi*fh*t);
        %Set up chest compression signal. Right now this are assuming this is the
        %phase information and will have the heart beat inside of it. The source
        %for the Adaptive Line Enhancer that this code is based off of uses only
        %the magnitude information.
        
        
        % sigPow = .001
        nt = sqrt(sigPow(nl))*randn(size(ct));
        rt = ct + nt;
        
        %These three lines create noise with power equal to sigPow. This is to be
        %adjusted to get some kind of ROC.
        %% Run LMS
        
        
        
        numIt = 1000;
        %This defines how many itterations to run the LMS for. If the entire signal
        %is to be filtered, set numIt equal to the length of the signal.
        % if numIt is too there is not enough data to FFT to get the resp rate. It
        % seems to need to be atleast 1000 data samples to get an accurate estimate
        % for the heart rate we need much less, like 400 will work. (400 at fs=100
        % is 4 seconds)
        
        N = 10;
        %N is the number of filter taps we are using for the LMS algorithm.
        
        W = zeros(N,numIt);
        y = zeros(1,numIt);
        %Initialize W and y to be empty vectors. W is the vector that will contain
        %the weights of the LMS filter. y is the vector that will contain the
        %filtered output of the LMS
        
        
        freqEst = nan(1,numIt);
        
        
        
        
        W(:,N) = 1;
        %Initialize the weight vector to 1. This is needed otherwise the
        
        
        pk=.1;
        mu = .001;
        lam = .99;
        gamma = .01;
        roa = .99;
        muMax = .01;
        muMin = .00001;
        % Initialize the parameters of the LMS filter
        % mu is the initial step size (.001) is good
        % lambda is the forgetting factor typacially between 0.9 and 0.99 (must be less than 1)
        % gamma is a factor that helps with convergence. Set to 0.01 ish
        % roa is a weighting factor, set it and forget it. (0.99)
        % muMax is the maxium step size mu can take. Bigger than around .01 is bad
        % muMin is the min step size. Should not matter much but .000001 works
        
        
        xn = rt;
        % The LMS taks xn as an input. I want the input to be the noisy
        % reflected signal
        
        
        Nd = 1;
        %The ALE algorithm requires the desired signal to be a delayed
        %version of the input signal xn. Set the number of dealys to 1
        
        
        dn = delayseq(xn',Nd);
        % Set the desired signal to the input signal but delayed by Nd
        % samples
        
        
        
        %For loop has to start at N+1 so that we can populate the Xk vector
        %below. We need to always be multiplying N values of xn by W
        for k = N+1:numIt
            %This line of code is a convenence. The last N values of xn get
            %multiplied by W vector
            
            Xk = [xn(k:-1:k-(N-1))]';
            
            
            %Filtering:
            % Filter the input xn values by our current weight vector to
            % get an estimate of the signal y at time k
            y(k) =  W(:,k-1)'*Xk;
            
            
            % Error Update:
            %find the error between the desired input and y. If y is equal
            %to the desired input then this will be 0
            %The LMS algorithm attempts to set W to minimize this value
            e(k) = dn(k) - y(k);
            
            % Covaranice update:
            %pk is a term that represents the covarance estimate of y. It
            %should be very small if y is converging to dn.
            pk = roa*pk + (1-roa)*e(k)*e(k-1);
            
            % Step Size Update:
            %Update the step size based on the covarance and the previous
            %step size
            mu = lam*mu + gamma * pk*conj(pk);
            
            %Also check if the new mu value is too big or too small. If it
            %is to big or small set it equal to the min/max value demending
            if mu> muMax
                mu = muMax;
            end
            
            if mu< muMin
                mu = muMin;
            end
            
            % Weight update:
            %Up date the weight vector based on previous weight vector and
            %the stepsize / error
            W(:,k) = W(:,k-1) + 2*mu*e(k) * y(k);
            
            
        end
        
        
        % The following lines of code are used to determine the effective
        % ness of the filter.
        
        %This is the FFT size to take. We need a large value to get a good
        %enough picture of the freqeuncies (1024 will probably work)
        Npt = 4096;
        
        %       Get the second order filter defined in the function getRespFilter
        %       It is a butterworth filter that has a freqeucny response that will
        %       extract only the respirtory information from the signals.
        Hd = getRespFilter(fs);
        sos = Hd.sos;
        
        %Filter the noisey Refected signal rt
        %Then normalize it since the filter has a very high gain
        rtResp = sosfilt(sos,rt);
        rtResp = rtResp/max(rtResp);
        
        %Filter the LMS estimate signal y to contain only the respirtory
        %signal
        %then normalize it since the filter has a high gain
        yResp = sosfilt(sos,y);
        yResp = yResp/max(yResp);
        
        
        %       Get the second order filter defined in the function getHeartFilter
        %       It is a butterworth filter that has a frequency response that will
        %       extract only the heart beat information from the signals.
        Hd = getHeartFilter;
        sos = Hd.sos;
        
        %Filter the noisey Refected signal rt to have only the heart beat
        %in it
        %Then normalize it since the filter has a very high gain
        frHeart = sosfilt(sos,rt);
        frHeart = frHeart/max(frHeart);
        
        %Filter the LMS estimate signal y to have only the heart beat in it
        %Normalize it
        yHeart = sosfilt(sos,y);
        yHeart = yHeart/max(yHeart);
        
        
        %Take the FFT of the resp rate filtered signals
        [frR,f] =easyFFT(rtResp,fs,Npt);
        [fyR,~] = easyFFT(yResp,fs,Npt);
        
        % Take the FFT of the heart rate filtered signals
        [frH,~] = easyFFT(frHeart,fs,Npt);
        [fyH,~] = easyFFT(yHeart,fs,Npt);
        
        %Find the index of the highest power part of the FFT
        % This should be the to the respitory signal only
        % then find the freqeuncy that correspondes to it it
        % This will be the estimated respitory rate in Hz.
        % This will get done for the raw reflected signal and the LMS
        % filtered signal
        [~,maxIdx] = max(abs(frR));
        noiseRespRate = abs(f(maxIdx));
        
        [~,maxIdx] = max(abs(fyR));
        filtRespRate = abs(f(maxIdx));
        
        
        %Find the index of the highest power part of the FFT
        % This should be the to the heart rate signal only
        % then find the freqeuncy that correspondes to it it
        % This will be the estimated heart rate in Hz.
        % This will get done for the raw reflected signal and the LMS
        % filtered signal
        [~,maxIdx] = max(abs(frH));
        noiseHeartRate = abs(f(maxIdx));
        
        [~,maxIdx] = max(abs(fyH));
        filtHeartRate = abs(f(maxIdx));
        
        
        %Save the error in a matrix
        %nre is the error of the respitory rate for the noisy rt signal
        %nhe is the error of the heart rate for the nosiy rt signal
        nre(nl,aIDX) = fr-noiseRespRate;
        nhe(nl,aIDX) = fh-noiseHeartRate;
        
        
        % fre is the error of teh respitory rate for the filtered signal y
        % fhe is the error of the heart rate for the filtered signal y
        fre(nl,aIDX) = fr-filtRespRate;
        fhe(nl,aIDX) = fh-filtHeartRate;
        
        
        
        
    end
end

%%
%These two lines calculated the average error. The result will be a sigle
%number for each sigPow. We can plot sigPow vs avgNH to give us the average
%error in heart rate for different noise inputs.
avgNH = mean(nre,2);
avgFH = mean(fre,2);

plot(sigPow,avgFH)
hold on
plot(sigPow,avgNH)
