function [clutterOut,clutterPlatformOut] = setupBackClutter(xyzPos,numTgt,rcsVector,fc)

clutterOut = cell(numTgt,1);
clutterPlatformOut = cell(numTgt,1);



if size(xyzPos,1) == 1
    xyzPos = repmat(xyzPos,numTgt,1);
end

if size(rcsVector,1) ~= numTgt
    rcsVector = rcsVector*ones(1,numTgt);
end

for k = 1:numTgt
    xyz = xyzPos(k,:);

    cluterTemp = phased.RadarTarget('Model','Nonfluctuating','MeanRCS',rcsVector(k), ...
        'PropagationSpeed',physconst('LightSpeed'),'OperatingFrequency',fc);
    
    clutterPlatTemp = phased.Platform('InitialPosition',xyz', ...
        'Velocity',[0;0;0]);
    
    clutterOut{k} = cluterTemp;
    clutterPlatformOut{k} = clutterPlatTemp;
    
end

end

