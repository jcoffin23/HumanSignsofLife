function [summedSig] = updateClutter(clutter,clutPlat,txpos,waveform,transmitter,radiator,channel,collector,receiver)
summedSig = 0;
for k = 1:length(clutPlat)
    platForm = clutPlat{k};
    clutTgt = clutter{k};
    
    % Update the target position
    clutpos = platForm.InitialPosition;
    clutvel = platForm.InitialVelocity;
    
%     [tgtpos,tgtvel] = platForm(TchestSample);
    % Get the range and angle to the target
    [tgtrng,tgtang] = rangeangle(clutpos,txpos);
    % Generate the pulse
    refSig = waveform();
    % Transmit the pulse. Output transmitter status
    [sig] = transmitter(refSig);
    % Radiate the pulse toward the target
    sig = radiator(sig,tgtang);
    % Propagate the pulse to the target in free space
    sig = channel(sig,txpos,clutpos,[0;0;0],clutvel);
    % Reflect the pulse off the target
    sig = clutTgt(sig);
    % Propagate the echo to the antenna in free space
    sig = channel(sig,clutpos,txpos,clutvel,[0;0;0]);
    % Collect the echo from the incident angle at the antenna
    sig = collector(sig,tgtang);
    % Receive the echo at the antenna when not transmitting
    rxsig = receiver(sig);
    rxsig = dechirp(rxsig,refSig);
    
    summedSig = summedSig + rxsig;
end

end

