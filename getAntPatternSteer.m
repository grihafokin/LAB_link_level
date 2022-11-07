% (c) Fokin G.A., Volgushev D.B., SPbSUT, 2022.
% function to calculate antenna array (AA) phase distribution vector
function w = getAntPatternSteer(antElPos, f, azAng, elAng)
% antElPos - array of AA elements coordinates [x,y,z], m
% f -        carrier frequency, Hz
% azAng, elAng - azimuth, elevation direction in which 
%                AA element phase distribution is calculated, degrees
 c = physconst('LightSpeed'); % speed of light, m/s
 % vector, that specifies the direction of signal arrival at the AA
 incidentDir = [-cosd(elAng).*cosd(azAng);...
                -cosd(elAng).*sind(azAng);...
                -sind(elAng)];
% calculation of the projection of the direction vector on the vectors, 
% connecting the origin coordinates and AA element position, 
% recalculation of this projection into delay
tau = antElPos*incidentDir/c;
% phase distribution calculation (in complex coefficient format)
% based on delay tau and operating frequency f
w = exp(-1i*2*pi*f*tau);
end