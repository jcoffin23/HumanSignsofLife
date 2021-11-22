% MONOSTATICCHIRP simulation of a single FMCW Chrip for a MONOstatic radar
% setup
%
% Syntax:
%  [rxSig] = monostaticChirp(transmitter,receiver,channel,radiator,collector,beamformer,tgt,txMotion,rxMotion,dt)
%
% Inputs:
%
%  Required:
%     transmitter: Radar Transmitter object
%
%     receiver: Radar Receiver object
%
%     channel: Channel object (Freespace likly)
%
%     radiator: Radar Radiator Object
%
%     collector: Radar Colelctor Object
%
%     beamformer: Beamformer (likly CBF)
%
%     tgt: (numTgt x 1 Array)
%         Array that contains the phased.platform objects for each target
%
%     txMotion: Phased.Platform Object
%
%     rxMotion: Phased.Platform Object
%
%     dt: (sclar)
%         amount of time that passes between each chirp
%
% Outputs:
%
%   rxSig (Nx1 double)
%     Recieved radar chrip wher N is the length of the chirp
%
% Written 7/6/21 Joe Coffin


function [rxSig,tgtPlat,rxPlat,txPlat] = monostaticChirp(targets,transmitter,receiver,txchannel,rxchannel,radiator,collector,beamformer,tgtPlat,rxPlat,txPlat,dt,waveform,numTgt,tgtStruct)



%Update Rx and TX movement (likly unmoving, but could move in future)
[txPos,txVel] = txPlat(dt);
[rxPos,rxVel] = rxPlat(dt);

%Update target movement for EACH target
[tgtPos,tgtVel] = tgtPlat(dt);
refSig = waveform();
txSig = transmitter(refSig); %Simulate propagation of puse

[tgtRng,tgtAng]=rangeangle(tgtPos,txPos); %Get angle as seen by transmiter

rxSig = radiator(txSig,tgtAng); %Radiate the targets
rxSig = txchannel(rxSig,txPos,tgtPos,txVel,tgtVel); %Simulate freepath loss
% sigtgt = zeros(length(txSig),numTgt);

for k = 1:numTgt
    %Angle from Tx to target
    [~,fwAng] = rangeangle(txPos,tgtPos(:,k));
    
    %Angle From Rx to target
    [rxRange,bkAng] = rangeangle(rxPos,tgtPos(:,k));
    %     fwAng(1) = 180-fwAng(1);
    %     bkAng(1) = 180+bkAng(1);
    %Calculate target reflections
    sigtgt(:,k) = targets{k}(rxSig(:,k));
    
end
%Simulate prop from target to recv
rxSig = rxchannel(sigtgt,tgtPos,rxPos,tgtVel,rxVel);

[~,rxAng]=rangeangle(tgtPos,rxPos);
rxSig = collector(rxSig,rxAng); %collect the signal

%Check if there should be beamforming, for bistatic with 1 recv this wont
%be true, with 2 recv this will be needed
if beamformer ~=0
    [theta,range]=cart2pol(tgtStruct.offset,tgtStruct.yoffset); %Calculate theta for beamforming
    
    rxSig = beamformer(rxSig,[rad2deg(theta);zeros(1,length(theta))]); %Beamform looking at each target angle and 0 eleveation
end

rxSig = receiver(rxSig); %Recieve the signal

rxSig = dechirp(rxSig,refSig); %Mix signals together to get beat freq signal

% rxSig = rxSig(:,1); %Only take from the first element to simulate Bistatic conditions.

end

