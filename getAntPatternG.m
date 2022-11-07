% (c) Fokin G.A., Volgushev D.B., SPbSUT, 2022.
% function to calculate antenna array (AA) gain
function g = getAntPatternG(antElPos, f, azAng, elAng, wS, backLobe)
% antElPos -     array of AA element coordinates [x,y,z], m
% f -            carrier frequency, Hz
% azAng, elAng - direction in azimuth, elevation, in which AA gain 
%                is calculated, degrees 
% wS -           beam steering vector of coefficients
% backLobe -     use backlobe suppression
% g -            AA gain in direction [azAng, elAng]

% beaforming - applying wS to coefficients, describing AA elements 
% phase distribution for a given direction  [azAng, elAng]
% w = wS'*getAntPatternSteer(antElPos, f, azAng, elAng);
w = sum(conj(wS).*getAntPatternSteer(antElPos, f, azAng, elAng));
% AA gain (amplitude)
g = abs(w);
% back-lobe suppression through reduction of AA gain 
% in the sector of azimuth angles from 90 to 270 degrees
if (backLobe == 1)
    if (azAng > 90 && azAng < 270)
        azAng = wrapTo360(azAng) - 90;
        % attenuation of AA gain in the sector of angles from 90 to 270
        % degrees is modeled through radius of ellipse with eccentricity
        % equal to 0.9999; the closer angle to 180, the higher suppression
        g = g*sqrt(cosd(azAng)^2 + (sqrt(1-0.9999^2)*sind(azAng))^2);
    end
end
end