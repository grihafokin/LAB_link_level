% (c) Fokin G.A., Volgushev D.B., SPbSUT, 2022.
% function to create eNB structure
function eNB = createNB(coords, antDir)
% coords - eNB [x,y,z] cordinates, m
% antDir - eNB antenna array (AA) orientation,
%          (axis of symmetry direction), [azimuth, tilt], degrees
eNB.Coords = coords.'; % eNB coordinates
eNB.AntDir = antDir.'; % eNB AA orientation
% rotation matrix according to antDir values, used for recalculation
% direction vectors from global coordinates to eNB AA local coordinates
eNB.AntOrient = rotz(antDir(1))*roty(-antDir(2));
eNB.Steer = 1;
end