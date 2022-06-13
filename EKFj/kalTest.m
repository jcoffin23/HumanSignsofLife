

runRealData =1;
if runRealData
  
  
    
  %%
clear 
close all

% load('EKFDat_HRV.mat')

%Get Spectrogram


% load('ekgKalmanINdat.mat')
load('hrvKalDat.mat')
%After expert filtering
flimited = f
corMaskTimeSplit=dIn;
[r,c,~]=find(corMaskTimeSplit == 1);
xAll = t(c);
yAll = flimited(r);
%%




end




%Sort the data to have increasing azRamp for all data
[timeAxis,srtIdx] = sort(xAll);
timeDelay = yAll(srtIdx);

%Determine the AzVector by getting the Unique theta vals from the data.
[ta,ia] = unique(timeAxis);

%Organize the data to be a matrix nxm where m is length of the az vector
%and n is the maxium number of estimates for any az vector. 
numEstimates = diff(ia);
thirdDem = max(numEstimates);
rngEst = nan(thirdDem,length(t));



verfiyDropouts = ismember(t,ta);
timeVecFull = nan(size(t));
timeVecFull(verfiyDropouts) = t(verfiyDropouts);
idx = 1;

%Put the esitames along the n axis so that for every az there is a number
%of possible range estimates. 
for outerIdx = 1:length(verfiyDropouts) -1
  if verfiyDropouts(outerIdx)
    if idx == length(ia)
%       rngEst(1:numEstimates(idx),outerIdx) = timeDelay(ia(idx) +1);
        doNothing = 1;
    else
      
      rngEst(1:numEstimates(idx),outerIdx) = timeDelay(ia(idx):ia(idx+1) - 1)   ;
      idx = idx +1;
    end
  end
  
end




datApend = timeDelay(ia(idx):end);
rngEst(:,end) = [datApend;nan(thirdDem-length(datApend),1)] ;


% TWEAKE PARAMS HERE

T = 4;
thetadot = 4;

dt = mean(T);

innFactor =  1 * 10e2;


%% Run Kalman Approach with PDAF

%Set initial Range to be based on first measurmetns. (This could be tuned)
xri = mean(rngEst(:,1),'omitnan');
%Set Initial Conditions
xInitial = [ xri,0,ta(1),thetadot]';
pInitial = diag([1,1,1,1]) ;




%Set up Q and R matricies
%Justy got to find the nums that work bb

rq = 10;
sigq = .001;



q = [rq,dt^3/2;dt^3/2,dt^2] .* sigq;
qq =[dt^4/4,dt^3/2;dt^3/2,dt^2] .* sigq;

Q = [q,zeros(2);zeros(2),qq];
R = diag([1,1]) ;

%Run Klamn filtering as a single function
[xperdict,xcov,~]=kalLoop(rngEst,t,xInitial,pInitial,R,Q,T,innFactor);
%Output grossAdjFeet as the range estimate of Kalman Filter
hrEst = xperdict(1,:);



if 1
  
  p = threedto2d(xcov);
  xx=xperdict([1,3],:);
  figure
  
  plot(xAll,yAll,'.')
  xlabel('Time')
  ylabel('Frequency Estimate')
  title('Cloud of Inputs vs Kal Solution')
  hold on
  plot(t,xx(1,:),'.')
  
%   plot(azDeg,y1)
%   plot(azDeg,dp)
  legend "Inputdata" "KalmanSolution"
  
  figure
  subplot 211
  plot(t,xx(1,:))
  subplot 212
  plot(t,p(:,1))
  
end







function [xp,pp,inn] = kalLoop(rngEst,azVec,xInitial,pInitial,R,Qin,T,innovationFactor)

Q = Qin;
xp =nan(4,length(azVec)); %Allocate Mat's for the kalman variables
xm = nan(4,length(azVec));
pp = nan(4,4,length(azVec));
pm = nan(4,4,length(azVec));
inn = zeros(2,length(azVec));
xp(:,1) = xInitial; %Set Initial Values
pp(:,:,1) = pInitial;



nz=1;
cNz= 2;

dropOutDet=0;
qHold=0;
%Probabilities
pd = .999;
pg = .99842;
Rin =  diag([1,1]);

beta0 = 0;
innFact = innovationFactor;


H = [1,0,0,0;
     0,0,1,0];


phi = [1,T,0,0;
       0,1,0,0;
       0,0,1,T;
       0,0,0,1];

  
for k = 2:length(azVec) 
  


  
  
  %% update based on previous postiotrito
  
  %Run State/Error Estimate Extrapolation

  xm(:,k) = phi*xp(:,k-1);
  pm(:,:,k) = phi*pp(:,:,k-1)*phi' + Q;
  sk = H * pm(:,:,k) * H' + R;
  
  
  
  
  %Predict measuremnt / calculate innovations
  zestMat = abs(rngEst(:,k) -   xm(1,k));
  zAll = rngEst(:,k);
  zAll = zAll(~isnan(zAll));
  zAll = [zAll, azVec(k) * ones(size(zAll))];
  
  vAll = zAll - (H * xm(:,k))';
  
  %Gate estiamtes based on innovation scores
  zi = zAll(abs(vAll(:,1)) < innFact,:);
  vi =  zi - (H * xm(:,k))';
  numZ = size(zi,1);
  
  
  
  if any(isempty(zAll)) %No measuremnts for this AZ, this is a dropout
    kGain= 0;
    xp(:,k) = xm(:,k);
    pp(:,:,k) = pm(:,:,k);
    dropOutDet = 1;
    innFact = .8;
%     Q = Q*100;
  elseif any(isempty(zi)) %Gate is empty, all measurements do not fit model 
    xp(:,k) = xm(:,k);
    pp(:,:,k) = beta0 * pm(:,:,k);
    
  else

    % If there is a dropout, turn the innovation higher.
    if dropOutDet ==1 && 0
      
      dropOutDet = 0;
      disp('Drop Out detected, Back on Track')
      qHold = 50;
    end
    if qHold >0
      qHold = qHold-1;
    else
      qHold =0;
      innFact = innovationFactor; 
      Q = Qin;
    end
    
 
    
    
    %% State Update
    
    %Kalculate Kalman gain
    kGain = pm(:,:,k)*H' * inv(sk);
    
    
    
    
    %Kalculate the PDAF parameters (see documention for explination)
    ei = exp(-.5 *vi * inv(sk)*vi'); 
    b = (2*pi/innovationFactor)^(nz/2) * numZ * (1/cNz) * ((1-pd*pg)/pd);
    
    if numZ >0
      beta = diag(ei)/(b+sum(diag(ei)));
      beta0 = b/(b+sum(diag(ei)));
      beta = [beta0;beta];
      
    else %NumZ is 0 and this measuremnt is not a dropout (hadled above with If)
      beta0=0;
      beta = 0; %since Numz is 0 b = 0 so beta = 0
    end
    
    
    
    combineInno = sum(beta(2:end).*vi,1)';
    inn(:,k) = combineInno;
    
    %Calculate the Updated State estiamte for this time sample
    xp(:,k) = xm(:,k) + kGain*combineInno;
    
    
    
    
    
    %% Calculate the  updated Covariance for this time sample 
    pc = pm(:,:,k) - kGain*sk*kGain';
    sumTerm = 0;
    
    for q = 1:numZ
      sumTerm = sumTerm + beta(q+1) .* vi(q,:)' * vi(q,:);
    end

    pTilda = kGain * (sumTerm - combineInno*combineInno') * kGain';
    
    pp(:,:,k) = beta(1)*pm(:,:,k) + (1-beta(1))*pc + pTilda;
   
  end
end


%% Run Backwards Pass of RTS smoother
N = length(xp);
xs = zeros(size(xp));
xs(:,N) = xp(:,N);


for k=N-2:-1:1
  A = pp(:,:,k) *phi'* inv(pm(:,:,k+1));
  xs(:,k) = xp(:,k) + A*(xs(:,k+1) - xm(:,k+1));
end

xp = xs;

end
