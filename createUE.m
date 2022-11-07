% (c) Fokin G.A., Volgushev D.B., SPbSUT, 2022.
% function to create UE structure
function ue = createUE(xPnts, yPnts, zPnts, serveNB, v, T, antDir)
% xPnts - initial and final values of the UE trajectory along x axis, m
% yPnts - initial and final values of the UE trajectory along y axis, m
% zPnts - constant along all trajectory UE height, m
% serveNB - serving eNB number, index in array of eNB structure
% v - UE speed, m/s 
% T - UE coordinate measurement period, s
% antDir - UE antenna array (AA) orientation 
%          (axis of symmetry direction), [azimuth, tilt], degrees
ue.ServeNB = serveNB; % serving eNB number
% creating an array of UE motion trajectory coordinates
trajPnts = [xPnts, yPnts];
trajPnts(:,3) = zPnts;
ue.Trajectory = getTrajectory(trajPnts, v, T);    
% initialization of the array of signal-to-interference ratio (SIR) values,
% SIR is calculated at each point of the UE trajectory
ue.SNR = zeros(length(ue.Trajectory(:,1)), 1);    
% rotation matrix according to antDir values; is used for recalculation
% direction vectors from global coordinates to UE AA locals coordinates
ue.AntOrient = rotz(antDir(1))*roty(-antDir(2));
% vector of AA directional coefficients (by default, AA is omnidirectional)
ue.Steer = 1;
end