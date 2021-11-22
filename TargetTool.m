%This Target tool creates the targets as specified by the input structure.
%an Example is shown below. 
% 
% numTgt = 2; Number of targets
% tgtStruct.numTgt = numTgt; 
% 
% tgtStruct.human = [1,0]; %Is each target breathing? needs to be 1xnumTgt
% tgtStruct.offset = [5,7]; %Xoffset of each target
% tgtStruct.yoffset= [0,0];% Yoffset of each target
% tgtStruct.fhActual=[]; %Emtpy vectors to pass the resulting freqeuncies into
% tgtStruct.frActual=[];
% 
% tgtStruct.scatterMatt = cat(3,[1 0;0 1],[0 1;1 0]); %Scattering matriceis, need to be one 2x2 mat per target. 

function [targets,targetplatform,tgtStruct,t] = TargetTool(tgtStruct,respHeight,heartHeight,fc,dt,t,fsChest,bistatic)

numTgt = tgtStruct.numTgt;







               % samples per second
               
               
        

targets = cell(numTgt,1);
wpts = [];

for k = 1:numTgt
    
    
    if tgtStruct.human(k)
        z = tgtStruct.zoffset(k)*ones(size(t'));
        
        %Get the chest compression signal. We The frequencies are random each time
        %and are drawn from a pool uniformly. We need to save them to compare the
        %estimated freqeuncies to the actual ones. This is our metric!
%         [~,ct,vt,tgtStruct.frActual(k),tgtStruct.fhActual(k)] = getChestCompression(fsChest,[],length(t),0,0,tgtStruct.offset(k),respHeight,heartHeight);
        
       [~,ct,~,tgtStruct.frActual(k),tgtStruct.fhActual(:,k)] = hrvChestComp(fsChest,[],length(t),0,0,tgtStruct.offset(k),respHeight,heartHeight,1);
        wpts = cat(3,wpts,[t';ct; tgtStruct.yoffset(k)*ones(size(t')); z]');%appeend this targets path to all targets paths
        
      
    else
        wpts = cat(3,wpts,[t';tgtStruct.offset(k)*ones(size(t')); tgtStruct.yoffset(k)*ones(size(t')); z']'); %append the stationary target to the wpts to make the platform
        tgtStruct.frActual(k) = 0;
        tgtStruct.fhActual(k) = 0;
    end
    


    
%This line creates the actual target object. The only thing that is
%differeent per target is the scattering matrix.
if bistatic
    targets{k} = phased.RadarTarget('Model','Nonfluctuating','Mode','Bistatic',...
        'PropagationSpeed',physconst('LightSpeed'),'OperatingFrequency',fc,'EnablePolarization',true,'ScatteringMatrix',tgtStruct.scatterMatt(:,:,k));
    
else
        
    meanHumanRCS = 0.1992; %From Paper in HSOL doc
    targets{k} = phased.RadarTarget('Model','Nonfluctuating','MeanRCS',meanHumanRCS,...
        'PropagationSpeed',physconst('LightSpeed'),'OperatingFrequency',fc,'EnablePolarization',false);
end

    
    
end
    oMat = cat(3,azelaxes(0,0),azelaxes(0,0));% Initial Orentation of each target, (0,0) is stairing down X axis I think
    oMat = azelaxes(0,0);
    targetplatform = phased.Platform('MotionModel','Custom','CustomTrajectory',wpts,'InitialOrientationAxes',oMat,'OrientationAxesOutputPort',true);
%     targetplatform.InitialPosition = zeros(3,2);
%     targetplatform.InitialVelocity = zeros(3,2);




end




