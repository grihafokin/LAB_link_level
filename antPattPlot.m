% (c) Fokin G.A., Volgushev D.B., SPbSUT, 2022.
% function to calculate and plot antenna array (AA) pattern
% antElPos -   array of coordinates [x,y,z] for AA element, m
% f -          carrier frequency, Hz
% nodeStruct - structure with gNB or UE link level parameters 
% angStep -    angle grid step in azimuth and elevation, degree
% scl -        scaling factor for AA pattern visualization
% backLobe -   use backlobe suppression (1-yes, 0-no)
% [x, y, z] -  coordinates for AA patter visualization, dimensionless 
function [x, y, z] = antPattPlot(antElPos, f, nodeStruct, angStep, scl, backLobe)
% grid of angles in azimuth and elevation for calculating AA pattern
azA = 0:angStep:360;  % azimuth
elA = -90:angStep:90; % elevation 
% array initialization to store AA pattern values
azN = length(azA);
elN = length(elA);
x = zeros(elN, azN);
y = zeros(elN, azN);
z = zeros(elN, azN);
for i=1:elN % loop through array of elevation angles
    for j=1:azN % cycle through array of azimuth angles
        % calculation of the AA pattern value for the i-th elevation angle 
        % and the j-th azimuth angle, taking into account vector of
        % AA steering coefficients (nodeStruct.Steer)
        p = getAntPatternG(antElPos, f, azA(j), elA(i), ...
            nodeStruct.Steer, backLobe)*scl;
        % recalculation of the AA pattern value from polar to rectangular 
        % coordinates, considering AA orientation (nodeStruct.AntOrient)
        r = p.*cosd(elA(i));
        xyz = nodeStruct.AntOrient*[r.*cosd(azA(j)); ...
            r.*sind(azA(j)); p.*sind(elA(i))];
        x(i,j) = xyz(1);
        y(i,j) = xyz(2);
        z(i,j) = xyz(3);
    end
end
end