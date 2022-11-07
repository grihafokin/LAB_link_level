% (c) Fokin G.A., Volgushev D.B., SPbSUT, 2022.
% function to calculate UE trajectory
function xyUN = getTrajectory(pnts, v, T)
% v - UE velocity, m/s
% total UE motion time, s
times = [0; sqrt(sum((pnts(2:end, :) - pnts(1:end-1, :)).^2,2))/v];
% time at each reference point/points, describing UE motion path
elapsedTime = zeros(1, length(times));    
for i=1:length(times)
    elapsedTime(i) = sum(times(1:i));
end
% built-in functions for generating coordinates
% of the trajectory on the given reference points
ts = trackingScenario('UpdateRate', 1/T);
target = platform(ts);
traj = waypointTrajectory('Waypoints', pnts, ....
    'TimeOfArrival', elapsedTime, 'SampleRate', 1/T);
target.Trajectory = traj;
r = record(ts);
posUN = [r(:).Poses];
xyUN = vertcat(posUN.Position); % coordinates of UE trajectory points
end