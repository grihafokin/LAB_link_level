% (c) Fokin G.A., Volgushev D.B., SPbSUT, 2022.
clear all; close all; clc;
% Link Level Simulation Model for Interference Evaluation in 5G 
% Millimeter-Wave Ultra-Dense Network with Location-Aware Beamforming;
% Part I. Two eNB and two moving UE
c = physconst('LightSpeed');
f = 30e9;      % carrier frequency, Hz
da = 0.5*c/f;  % distance betweeb antenna array (AA) elements, m 
snrThr = 10;   % signal-ro-noise-ratio (SNR) threshold value, dB
anim = 0;      % 1 - enable animation
stdCoords = 0; % standard deviation of UE [x,y,z] coordinates estimation, m
useAntUE = 0;  % enable use of AA at the UE; 1-yes, 0-no
backLobe = 1;  % enable backlobe suppression of AA pattern (useful for URA)
Ta=0; % beamforming (BF) period, s; if Ta=0, BF performed every time step
AAtype=1;      % choose AA type: 1-URA, 2-ULA, 3-UCA
Nel = 4;       % number of AA elements in one dimension
antElPos = createAnt(AAtype, Nel, da); % creation of AA elements positions
% chose scenario of simulation:
    % case1: eNBs are on the same side of UE movement 
    % trajectory; UEs move in parallel at a distance d
    % case 2: eNBs are located on opposite sides of the UE 
    % movement trajectory; UEs move in parallel at a distance d
    % case 3: eNBs are located on the same side of the UE 
    % movement trajectory; UEs move one after another at a distance d
[eNB, ueNode, d, T, v] = createScenarion(3); % chose scenario of simulation
Nd = length(d); % 
Nnb = length(eNB);    % eNB number
Nue = length(ueNode); % UE number
N = length(ueNode(1).Trajectory(:,1)); % number of calculation points 
                              % (number of UE trajectory coordinate points)
trajArray = [ueNode.Trajectory];       % UEs coordinates array [N x 3*Nue]
eNBcoords = [eNB(:).Coords].';         % eNBs coordinates array[Nnb x 3]

if (anim == 1)
    antPattScl = 0.8;   % scalling coefficient for eNB AA pattern animation 
    antPattSclUe = 0.6; % scalling coefficient for UE AA pattern animation 
    fg = figure(5); 
    fg.WindowState = 'maximized'; 
    grid on; hold on;
    indStart = N/Nd*9;
    % display eNB AA pattern (see in antPattPlot)
    eNBptrnPlot = gobjects(1, 2);
    for i=1:Nnb % loop through eNB number
        [x, y, z] = antPattPlot(antElPos,...
            f, eNB(i), 6, antPattScl, backLobe);
        eNBptrnPlot(i) = surf(x+eNB(i).Coords(1),...
                        y+eNB(i).Coords(2),...
                        z+eNB(i).Coords(3), 'FaceColor', '#4DBEEE');
    end
    
    % display UE AA pattern (see in antPattPlot)
    if (useAntUE == 1)
        ueptrnPlot = gobjects(1, 2);
        for i=1:Nue
            [x, y, z] = antPattPlot(antElPos,...
                f, ueNode(i), 6, antPattSclUe, backLobe);
            ueptrnPlot(i) = surf(x+ueNode(i).Trajectory(indStart,1),...
                            y+ueNode(i).Trajectory(indStart,2),...
                            z+ueNode(i).Trajectory(indStart,3),...
                            'FaceColor', '#4DBEEE');
        end
    end

    % display UE position
    uePlot = gobjects(1, 2);
    ueText = gobjects(1, 2);
    for i=1:Nue % loop through UE number
        uePlot(i) = plot3(ueNode(i).Trajectory(1,1),...    
                          ueNode(i).Trajectory(1,2),...
                          ueNode(i).Trajectory(1,3), '^', ...
                          'MarkerSize', 10);
        ueText(i) = text(ueNode(i).Trajectory(1,1), ...
                         ueNode(i).Trajectory(1,2), ...
                         ueNode(i).Trajectory(1,3), ...
            sprintf('UE_{%i}', i), 'FontSize', 14, 'Color', '#A2142F');
    end
    
    ueDirPlot  = gobjects(1, 2);
    ueDirPlot2 = gobjects(1, 2);
    indSnoi = [2,1];
    for i=1:Nnb % loop through eNB number
        % display vector from eNB to SOI UE
        ueDirPlot(i)=plot3([eNB(i).Coords(1);ueNode(i).Trajectory(1,1)],...
                           [eNB(i).Coords(2);ueNode(i).Trajectory(1,2)],...
                           [eNB(i).Coords(3);ueNode(i).Trajectory(1,3)],...
                           'Color', '#76AB2F');
        % display vector from eNB to SNOI UE
        ueDirPlot2(i) = ...
            plot3([eNB(i).Coords(1);ueNode(indSnoi(i)).Trajectory(1,1)],...
                  [eNB(i).Coords(2);ueNode(indSnoi(i)).Trajectory(1,2)],...
                  [eNB(i).Coords(3);ueNode(indSnoi(i)).Trajectory(1,3)],...
                  'Color', '#D95319');
    end
    ueDirPlot(2).LineStyle = '--'; 
    ueDirPlot2(2).LineStyle = '--';

    for i=1:Nnb % loop through eNB number  
        text(eNB(i).Coords(1), eNB(i).Coords(2), eNB(i).Coords(3)+5, ...
            sprintf('gNB_{%i}', i), 'FontSize', 16, 'Color', '#A2142F');
    end
    
    xlabel('x, m'); ylabel('y, m'); axis equal;
    axis([-5, 155, -5, 65, 0, 30]); 
    view([0, 90]);
else
    indStart = 1;
end

% initialize arrays to store angles of departure (AOD) from each eNB to
% each UE (required to resolve the error reading a non-existent array when
% in some modes and scenarious)
azAng = zeros(Nue, Nnb);
elAng = zeros(Nue, Nnb);
azAngUE = zeros(Nue, Nnb);
elAngUE = zeros(Nue, Nnb);

% main processing loop
for i=indStart:N % loop through number of UE trajectory coordinate points
    % array of coordinates of all UEs for the i-th calculation point
    ueCoordsi = reshape(trajArray(i, :).', 3, Nue).';
    % apply error in UE coordinate estimation according to stdCoords
    ueCoordsiErr = ueCoordsi + stdCoords*randn(size(ueCoordsi));
    
    if (mod(i, round(Ta/T)) == 1 || Ta == 0)
        % initialize arrays to store AOD from each eNB to each UE
        azAng = zeros(Nue, Nnb); % azimuth AOD
        elAng = zeros(Nue, Nnb); % elevation AOD

        % initialize arrays to store angles of arrival (AOA) for each UE 
        % from each eNB; used if there is an AA on the UE, if useAntUE = 0
        azAngUE = zeros(Nue, Nnb); % azimuth AOA
        elAngUE = zeros(Nue, Nnb); % elevation AOA

        for j=1:Nue % loop through all UE 
            % vector specifying the direction from eNB 
            % to UE in global coordinate system x,y,z
            diffCoord = ueCoordsiErr(j,:) - eNBcoords;        
            for n=1:Nnb % loop through all eNB
                % vector specifying the direction from the eNB to the UE 
                % in local coordinate system of the eNB AA,
                % i.e. considering the position of the eNB antenna array
                dirVect = eNB(n).AntOrient.'*diffCoord(n,:).';
                % calculate AODs from n-th eNB to j-th UE 
                azAng(j, n) = rad2deg(atan2(dirVect(2), dirVect(1)));
                elAng(j, n) = rad2deg(atan2(dirVect(3), ...
                    sqrt(sum(dirVect(1:2).^2))));        

                % calculation of the AA direction vector coefficients 
                % for the n-th eNB, serving j-th UE; serving eNB number 
                % is indicated in the parameter ServeNB per each UE
                if (n == ueNode(j).ServeNB)
                    eNB(n).Steer = getAntPatternSteer(antElPos, f,...
                        azAng(j, n), elAng(j, n));
                end

                % calculate AOD and AA direction vector coefficients for UE
                if ( useAntUE == 1)
                    % vector, giving direction from UE to eNB in the local 
                    % coordinate system of the UE AA, 
                    % i.e. considering the position of the UE antenna array
                    dirVect = -ueNode(j).AntOrient.'*diffCoord(n,:).';
                    % calculate AODs from j-th UE to n-th eNB
                    azAngUE(j, n) = rad2deg(atan2(dirVect(2), dirVect(1)));
                    elAngUE(j, n) = rad2deg(atan2(dirVect(3), ...
                        sqrt(sum(dirVect(1:2).^2))));
                    % calculation of the AA direction vector coefficients
                    % for the j-th eNB, working with n-th UE
                    if (n == ueNode(j).ServeNB)
                        ueNode(j).Steer = getAntPatternSteer(antElPos, ...
                            f, azAngUE(j,n), elAngUE(j,n));
                    end
                end
            end
        end
    end
    
    % array for temporary storage of values 
    % of power, received by UE from each eNB
    eNBpwr = zeros(Nnb, 1);

    % calculation of the ratio of power, received from serving eNB (SOI) 
    % to power, received from neighbor eNB (SNOI) for each UE
    for j=1:Nue % loop through all UE 
        for n=1:Nnb % loop theough all eNB 
            % calculate UE receiving AA gain
            if ( useAntUE == 1 && j == 1)
                % UE gain with AA, which beam is oriented to serving eNB
                gUE = getAntPatternG(antElPos, f, ...
                    azAngUE(n), elAngUE(n), ueNode(1).Steer, backLobe).^2;
            else
                % UE gain without AA
                gUE = 1;
            end
            % calculation of power, received by j-th UE from the n-th eNB, 
            % taking into account beamforming at UE & eNB, excluding range
            eNBpwr(n) = gUE*getAntPatternG(antElPos, f,...
                azAng(j, n), elAng(j, n), eNB(n).Steer, backLobe).^2;
        end
        % calculate distance from j-th UE to each eNB
        diffCoord = ueCoordsi(j,:) - eNBcoords;
        distSpace = sqrt(sum(diffCoord.^2,2));
        % calculation of power, received by j-th UE from each eNB, taking 
        % into account range; path loss is calculated from the FSPL
        % (free-space path loss) attenuation model
        eNBpwr = pow2db(eNBpwr) - fspl(distSpace,c/f);
        % расчет отношения принимаемой мощности от обслуживающей eNB 
        % к мощности принимаемой от соседней eNB для j-й UE

        % calculation of ratio of received power from serving eNB (SOI)
        % to received power from neighbor eNB (SNOI) for j-th UE
        ueNode(j).SNR(i) = eNBpwr(ueNode(j).ServeNB) - ...
            sum(eNBpwr(1:end ~= ueNode(j).ServeNB));
    end
    
    % update the display of the position of the UE 
    % and eNB/UE AA pattern for the current sample time
    if (anim == 1)
        for ip=1:Nnb
            % update eNB AA pattern
            [x, y, z] = antPattPlot(antElPos, f,...
                eNB(ip), 6, antPattScl, backLobe);
            eNBptrnPlot(ip).XData = x + eNB(ip).Coords(1);
            eNBptrnPlot(ip).YData = y + eNB(ip).Coords(2);
            eNBptrnPlot(ip).ZData = z + eNB(ip).Coords(3);
            % update UE AA pattern
            if (useAntUE == 1)
                [x, y, z] = antPattPlot(antElPos, f, ...
                    ueNode(ip), 6, antPattSclUe, backLobe);
                ueptrnPlot(ip).XData = x + ueCoordsi(ip,1);
                ueptrnPlot(ip).YData = y + ueCoordsi(ip,2);
                ueptrnPlot(ip).ZData = z + ueCoordsi(ip,3);
            end
            % update direction vector from eNB to (SOI) UE
            ueDirPlot(ip).XData = [eNB(ip).Coords(1); ueCoordsi(ip,1)];
            ueDirPlot(ip).YData = [eNB(ip).Coords(2); ueCoordsi(ip,2)];
            ueDirPlot(ip).ZData = [eNB(ip).Coords(3); ueCoordsi(ip,3)];
            % update direction vector from eNB to (SNOI) neighbour UE
            ueDirPlot2(ip).XData = ...
                [eNB(ip).Coords(1); ueCoordsi(indSnoi(ip),1)];
            ueDirPlot2(ip).YData = ...
                [eNB(ip).Coords(2); ueCoordsi(indSnoi(ip),2)];
            ueDirPlot2(ip).ZData = ...
                [eNB(ip).Coords(3); ueCoordsi(indSnoi(ip),3)];
            % update UE position
            uePlot(ip).XData = ueCoordsi(ip,1);
            uePlot(ip).YData = ueCoordsi(ip,2);
            uePlot(ip).ZData = ueCoordsi(ip,3);
            ueText(ip).Position = [ueCoordsi(ip,1)+2, ...
                                   ueCoordsi(ip,2), ...
                                   ueCoordsi(ip,3)+5];
        end      
        pause(0.001)
    end % if (anim == 1)
end % loop through number of UE trajectory coordinate points

% prepare arrays of coordinates and array of the 
% signal-to-interference ratio (SIR) values for the 1st UE to display
X = reshape(ueNode(1).Trajectory(:,1), [], Nd).';
Y = reshape(ueNode(1).Trajectory(:,2), [], Nd).';
Z = reshape(ueNode(1).SNR, [], Nd).';

% SIR map in each 1st UE trajectory point
figure(1); surf(X, Y, Z, 'FaceColor', 'interp', 'EdgeColor','none');
grid on; xlabel('x, m'); ylabel('y, m'); view([0, 90]);
c1 = colorbar; c1.Label.String = 'SIR, dB';

% map of the 1st UE position points, 
% in which the SIR exceeds the SNR threshold snrThr
figure(2); surf(X, Y, double(Z>snrThr), 'EdgeColor','none');
grid on; xlabel('x, m'); ylabel('y, m'); view([0, 90]);
colormap(winter(2)); c3 = colorbar;
c3.Label.String = sprintf('SIR > %.0f dB', snrThr);
c3.Ticks = [0, 1]; view([0, 90]);

% SIR map, displaying eNB AA position and orientation
figure(3); surf(X, Y, Z, 'FaceColor', 'interp', 'EdgeColor','none');
grid on; hold on; 
for n=1:Nnb
    absAntCoord = (eNB(n).AntOrient*[0, -1, 0; 0, 1, 0].'*3 + eNB(n).Coords).';
    plot3(absAntCoord(:,1), absAntCoord(:,2), ...
        absAntCoord(:,3), 'Color', '#ECB01F', 'LineWidth', 3)
    text(eNB(n).Coords(1), eNB(n).Coords(2)*1.06, ...
        eNB(n).Coords(3), sprintf('Ant gNB %i', n));
end
xlabel('x, m'); ylabel('y, m'); axis equal; view([0, 90]);
c2 = colorbar; c2.Label.String = 'SIR, dB';

%% SCENARIO FUNCTION
% function to create scenario of eNB positions and UE trajectories
function [eNB, ueNode, d, T, v] = createScenarion(sceneN)
T = 0.1;  % measurement period, s 
v = 10;   % UE velocity, m/s
switch sceneN
    % case1: eNBs are on the same side of UE movement 
    % trajectory; UEs move in parallel at a distance d
    case 1 
        % distance between the motion paths of the two UEs;
        % in this scenario, UE trajectories are spaced by d in y axis
        d = 0:1:10;
        % create array of two eNB structures; see createNB for eNB params
        eNB(1) = createNB([25, 50, 5], [-90, -1]);
        eNB(2) = createNB([125, 50, 5], [-90, -1]);
        % create array of two UE structures; see createUE for UE params
        ueNode(1) = createUE([0; 150], [0; 0], 0, 1, v, T, [90, 0]);
        ueNode(2) = createUE([150; 0], [0; 0], 0, 2, v, T, [90, 0]);             
        % building of a set of UE motion trajectories for different values 
        % of d; UE2 trajectory does not change, UE1 trajectory shifts by d
        Nd = length(d); % number of d values
        trajUE1 = zeros(length(ueNode(1).Trajectory(:,1)), 3, Nd);
        trajUE2 = zeros(length(ueNode(1).Trajectory(:,1)), 3, Nd);
        for i=1:Nd
            trajUE1(:,:,i) = ueNode(1).Trajectory + [0, d(i), 0];
            trajUE2(:,:,i) = ueNode(2).Trajectory;
        end
        % group a set of trajectories into 
        % one array of coordinates for each UE
        trajUE1 = permute(trajUE1, [1 3 2]);
        trajUE1 = reshape(trajUE1, [],3 ,1);
        trajUE2 = permute(trajUE2, [1 3 2]);
        trajUE2 = reshape(trajUE2, [],3 ,1);
        % save trajectory in UE structure ueNode
        ueNode(1).Trajectory = trajUE1;
        ueNode(2).Trajectory = trajUE2;
            
    % case 2: eNBs are located on opposite sides of the UE 
    % movement trajectory; UEs move in parallel at a distance d
    case 2     
        % distance between the motion paths of the two UEs;
        % in this scenario, UE trajectories are spaced by d in y axis
        d = 0:1:10;
        % create array of two eNB structures
        eNB(1) = createNB([25, 50, 5], [-90, -1]);
        eNB(2) = createNB([125, 0, 5], [90, -1]);
        % create array of two UE structures
        ueNode(1) = createUE([0; 150], [20; 20], 0, 1, v, T, [90, 0]);
        ueNode(2) = createUE([150; 0], [20; 20], 0, 2, v, T, [90, 0]);
        % building of a set of UE motion trajectories for different values 
        % of d; UE2 trajectory does not change, UE1 trajectory shifts by d
        Nd = length(d); % number of d values
        trajUE1 = zeros(length(ueNode(1).Trajectory(:,1)), 3, Nd);
        trajUE2 = zeros(length(ueNode(1).Trajectory(:,1)), 3, Nd);
        for i=1:Nd
            trajUE1(:,:,i) = ueNode(1).Trajectory + [0, d(i), 0];
            trajUE2(:,:,i) = ueNode(2).Trajectory;
        end
        % group a set of trajectories into 
        % one array of coordinates for each UE
        trajUE1 = permute(trajUE1, [1 3 2]);
        trajUE1 = reshape(trajUE1, [],3 ,1);
        trajUE2 = permute(trajUE2, [1 3 2]);
        trajUE2 = reshape(trajUE2, [],3 ,1);
        % save trajectory in UE structure ueNode
        ueNode(1).Trajectory = trajUE1;
        ueNode(2).Trajectory = trajUE2;
  
    % case 3: eNBs are located on the same side of the UE 
    % movement trajectory; UEs move one after another at a distance d
    case 3
        % distance between the motion paths of the two UEs;
        % in this scenario, UE trajectories are spaced by d in x axis
        d = 0:2:50;
        % create array of two eNB structures
        eNB(1) = createNB([25, 50, 5], [-90, -1]);
        eNB(2) = createNB([125, 50, 5], [-90, -1]);
        % create array of two UE structures
        ueNode(1) = createUE([0; 150], [0; 0], 0, 1, v, T, [90, 0]);
        ueNode(2) = createUE([0; 150], [0; 0], 0, 2, v, T, [90, 0]);
        % building of a set of UE motion trajectories for different values 
        % of d; UE2 trajectory does not change, UE1 trajectory shifts by d
        Nd = length(d); % number of d values
        trajUE1 = zeros(length(ueNode(1).Trajectory(:,1)), 3, Nd);
        trajUE2 = zeros(length(ueNode(1).Trajectory(:,1)), 3, Nd);
        for i=1:Nd
            trajUE1(:,:,i) = ueNode(1).Trajectory + [d(i), -(i-1)/100, 0];
            trajUE2(:,:,i) = ueNode(2).Trajectory;
        end
        % group a set of trajectories into 
        trajUE1 = permute(trajUE1, [1 3 2]);
        trajUE1 = reshape(trajUE1, [],3 ,1);
        trajUE2 = permute(trajUE2, [1 3 2]);
        trajUE2 = reshape(trajUE2, [],3 ,1);
        % save trajectory in UE structure ueNode
        ueNode(1).Trajectory = trajUE1;
        ueNode(2).Trajectory = trajUE2;
end
end