%replicaz generate  a  matrix  of  replica  vectors  for  a  linear  arrayoriented along thez-axi
%   Take in zspacing of sensors and kz freqs to comput v(kz)
%   Output V(kz)

function [V,kzOut] = replicaz(zn,spatialIn,selSpatial,lam,d)
arguments
    zn (1,:) {mustBeNumeric,mustBeReal}
    spatialIn (1,:) {mustBeNumeric,mustBeReal}
    selSpatial (1,:) char {mustBeMember(selSpatial,{'kz','psi','u','theta'})} = 'kz'
    lam = 0
    d = 0
end

n = 0:length(zn) -1;
d = zn(2) - zn(1);
V = nan(length(zn),length(spatialIn));

switch selSpatial
    case 'psi'
        if d == 0
            error('Spacing not defined')
        end
        kz = -d^-1 * spatialIn;
    case 'u'
        if lam == 0
            error('Wavelength not defined')
        end
        kz = -(2*pi)/lam * spatialIn; 
    case 'theta'
        kz = -2*pi/lam * cosd(spatialIn);
        
    case 'kz'
        kz = spatialIn;
end

kzOut.kz = kz;

if d~=0                    % only fill in psi if d exists   
    kzOut.psi=-kzOut.kz*d;
end
kzOut.uz=-lam/(2*pi)*kzOut.kz;
kzOut.theta=acosd(kzOut.uz);

% Set the values of kzarg.theta outside the visible region equal to NaN
% so they won’t be plotted
notvisibleind=find(imag(kzOut.theta) ~= 0);
kzOut.theta(notvisibleind)=NaN;

pz = d*n;
V=exp(-j*pz'*kz);                        % compute replicas

end

