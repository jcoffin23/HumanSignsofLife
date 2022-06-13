function [xOut,xCov,innO, de ] = kalLoop(zMe,xInitial,pInitial,H,R,Q,tv)


xOut = zeros(4,length(zMe));
xCov = zeros(4,4,length(zMe));
wasjustNan = 0;

xPlus = xInitial;
pPlus = pInitial;

innO = zeros(2,1,length(zMe) - 4);
%Kalman update Loop
FF = 0;
innovationFactor = 100;
for kIdx = 1:length(zMe) - 4
  
  T = tv;
  
  rang = xPlus(1);

 
  
  F = [1,T,0,0;
    0,1 ,0,0;
    0,0,1,T;
    0,0,0,1];
  
  
  
  zk = zMe(kIdx,:)';
  
  xHatMinus = F * xPlus; 
  
  %   Q = randn(size(Q));
  %   Q = zeros(size(Q));
  
  pMinus = F * pPlus * F' + Q;
  
   
  
  zHat =  H * xHatMinus; % This should be  the same as the lines below
  
  
  inn = zk- zHat;
  
  
  
  
  if any(isnan(zk))
    xPlus = xHatMinus;
    pPlus = pMinus;
    inn = 0;
    
    wasjustNan = 1;
  elseif  abs(inn(1))> innovationFactor
    xPlus = xHatMinus;
    pPlus = pMinus;
    inn = 0;
    
    
  else
    
    
    
    
    %% UPdate
    
    
    
    sk = H * pMinus * H' + R;
    
    
    de = inn' * inv(sk) * inn;
    
    
    kMat = (pMinus * H') * inv(sk) ;
    %     Kmat = sk /(pMinus * H');
    
    
    
    
    pPlus = (eye(4) - kMat*H) * pMinus;
    
    xPlus = xHatMinus + kMat*inn;
    wasjustNan = 0;
  end
  %%
  xOut( :, kIdx ) = xPlus;
  
  xCov( :, :, kIdx ) = pPlus;
  
  innO (:,:,kIdx) = inn;
  
end



end

