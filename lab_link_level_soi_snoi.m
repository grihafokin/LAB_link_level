% (c) Fokin G.A., Volgushev D.B., SPbSUT, 2022.
% clear all; close all; clc;
% Link Level Simulation Model for Interference Evaluation in 5G 
% Millimeter-Wave Ultra-Dense Network with Location-Aware Beamforming;
% Part II. One eNB and two stationary UE
c = physconst('LightSpeed');
f = 30e9;      % carrier frequency, Hz
da = 0.5*c/f;  % distance betweeb antenna array (AA) elements, m 
snrThr = 15;   % signal-ro-noise-ratio (SNR) threshold value, dB
stdCoords = 0; % standard deviation of UE [x,y,z] coordinates estimation, m
useAntUE = 0;  % использовать антенную решетку на UE
backLobe = 0;  % enable backlobe suppression of AA pattern (useful for URA)
AAtype=1;      % choose AA type: 1-URA, 2-ULA, 3-UCA
Nel = 8;       % number of AA elements in one dimension
antElPos = createAnt(AAtype, Nel, da); % creation of AA elements positions

T = 0.5;        % measurement period, s
v = 1;          % UE velocity, m/s
d = 0:1:100;    % distance between UE and eNB, m 
alpha = 0:1:60; % angular separation between UEs, degre

% create eNB parameters structure (see createNB)
eNB(1) = createNB([0, 0, 15], [90, -1]);
% create UE parameters structure (see createUE)
ueNode(1) = createUE([0; 15], [0; 0], 0, 1, v, T, [-90, 0]);
ueNode(2) = createUE([0; 15], [0; 0], 0, 1, v, T, [0 ,0]);

% creating an array of coordinates of UE location points, where
% each value in alpha corresponds to a set of values d; for each 
% separation angle the calculation is carried out for all ranges d

Nd = length(d);     % number of calculation points w.r.t. range  d
Na = length(alpha); % number of calculation points w.r.t. angle alpha
trajUE1t = zeros(Nd, 3, Na);
trajUE2t = zeros(Nd, 3, Na);
for i=1:Na
    trajUE1t(:,1:2,i) = [( sind(alpha(i)).*d).', (cosd(alpha(i))*d).'];
    trajUE2t(:,1:2,i) = [(-sind(alpha(i)).*d).', (cosd(alpha(i))*d).'];
end

% group a set of UE location point into
% one common array of coordinates for each UE
trajUE1 = permute(trajUE1t, [1 3 2]); trajUE1 = reshape(trajUE1, [],3 ,1);
trajUE2 = permute(trajUE2t, [1 3 2]); trajUE2 = reshape(trajUE2, [],3 ,1);
% save trajectory in UE structure ueNode
ueNode(1).Trajectory = trajUE1;
ueNode(2).Trajectory = trajUE2;

Nnb = length(eNB);      % number of eNB
Nue = length(ueNode);   % number of  UE
N = length(ueNode(1).Trajectory(:,1)); % number of calculation points
trajArray = [ueNode.Trajectory];      % UE coordinates array of [N x 3*Nue]
eNBcoords = [eNB(:).Coords].';        % eNB coordinates array of [Nnb x 3]
% инициализация массива направляющий коэфф. eNB (для расчета SOI/SNOI 
% eNB имеет два вектора коэфф., для направления луча на UE1 и UE2)

% initialization of array for eNB steering vector; to calculate SOI/SNOI
% eNB has two vectors of coefficient, for directing the beam to UE1 and UE2
eNB.Steer = zeros(size(antElPos,1), Nue);

% main processing loop
for i=1:N % loop through number of UE coordinate points
    % array of coordinates of all UEs for the i-th calculation point
    ueCoordsi = reshape(trajArray(i, :).', 3, Nue).';
    % apply error in UE coordinate estimation according to stdCoords
    ueCoordsiErr = ueCoordsi + stdCoords*randn(size(ueCoordsi));
    
    % initialize arrays to store AOD from eNB to each UE
    azAng = zeros(Nue, 1); % azimuth AOD 
    elAng = zeros(Nue, 1); % elevation AOD
    
    for j=1:Nue % loop through the number of UE
        % vector specifying the direction from eNB 
        % to UE in global coordinate system x,y,z
        diffCoord = ueCoordsiErr(j,:) - eNBcoords;        
        % vector specifying the direction from the eNB to the UE 
        % in local coordinate system of the eNB AA,
        % i.e. considering the position of the eNB antenna array
        dirVect = eNB.AntOrient.'*diffCoord.';
        % calculate AODs from n-th eNB to j-th UE 
        azAng(j) = rad2deg(atan2(dirVect(2), dirVect(1)));
        elAng(j) = rad2deg(atan2(dirVect(3), sqrt(sum(dirVect(1:2).^2))));      
        % calculation of eNB AA direction vector coefficients for j-th UE
        eNB.Steer(:,j) = getAntPatternSteer(antElPos, f, azAng(j), elAng(j));
        % calculate AODs and UE AA direction vector coefficients
        if ( useAntUE == 1 && j == 1)
            % vector, giving direction from UE to eNB in the local 
            % coordinate system of the UE AA, 
            % i.e. considering the position of the UE antenna array
            dirVect = -ueNode(j).AntOrient.'*diffCoord.';
            % calculate AODs from j-th UE to eNB 
            azAngUE=rad2deg(atan2(dirVect(2), dirVect(1)));
            elAngUE=rad2deg(atan2(dirVect(3),sqrt(sum(dirVect(1:2).^2))));
            % calculation AA direction vector 
            % coefficients of the j-th UE working with eNB
            ueNode(j).Steer = ...
                getAntPatternSteer(antElPos, f, azAngUE, elAngUE);
        end
    end % for j=1:Nue
    % array for temporary storage of values 
    % of power, received by UE from each eNB
    eNBpwr = zeros(1, Nue);
    % calculate UE receive AA gain
    if (useAntUE == 1)
        % UE gain with antenna array
        gUE = getAntPatternG(antElPos, f,...
            azAngUE, elAngUE, ueNode(1).Steer, backLobe).^2;
    else
        % UE gain without antenna array
        gUE = 1;
    end
       
    % calculation of power, received by j-th UE from eNB, 
    % taking into account beamforming at UE & eNB
    for j=1:Nue
        eNBpwr(j) = pow2db(gUE*getAntPatternG(antElPos, f,...
            azAng(1), elAng(1), eNB.Steer(:,j), backLobe).^2);
    end
    % calculate SOI/SNOI ratio
    ueNode(1).SNR(i) = eNBpwr(1) - eNBpwr(2);
end % for i=1:N

% prepare arrays of coordinates and array of the SOI/SNOI
% signal-to-interference ratio (SIR) values for the 1st UE to display
X = repmat(d, Na, 1);
Y = repmat(alpha, Nd, 1).';
Z = reshape(ueNode(1).SNR, Nd, Na).';
Z(Z>50) = 50; % limit max value for better visualization

% SIR map for all values of range d and angle separation alpha
figure(1); surf(X, Y, Z, 'FaceColor', 'interp', 'EdgeColor','none');
grid on; xlabel('d, m'); ylabel('\alpha, \circ');
view([0, 90]); c1 = colorbar; c1.Label.String = 'SIR, dB';

% map of the 1st UE position points, 
% in which the SIR exceeds the SNR threshold snrThr
figure(2); surf(X, Y, double(Z>snrThr), 'EdgeColor','none');
grid on; xlabel('d, m'); ylabel('\alpha, \circ'); 
view([0, 90]); colormap(winter(2)); c3 = colorbar;
c3.Label.String = sprintf('SIR > %.0f dB', snrThr);
c3.Ticks = [0, 1];

% slice (срез) of SIR dependency for given d=ds and all alpha values
figure(3)
ds = 60; % range value for which slice is plotted
% нахождение ближайшего к ds значений доступных d
% find d, closest to ds, for plotting 
[~,ind] = min(abs(X(1,:) - ds)); ds = X(1,ind);
plot(Y(:, ind), Z(:, ind)); grid on;
ylabel('SIR, dB'); xlabel('\alpha, \circ');
title(sprintf('SIR, dB for d=%0.0f, m', ds)); hold on;

% visualize simulated scenario of SOI and SNOI 
figure(4);
alphPlt = 20; % angular separation for visualization
[~,alphInd] = min(abs(alpha - alphPlt));
plot(trajUE1t([1,Nd],1,alphInd), trajUE1t([1,Nd],2,alphInd)); hold on;
plot(trajUE2t([1,Nd],1,alphInd), trajUE2t([1,Nd],2,alphInd));
plot([0;0], [0; d(end)], '--');

plot(0, 0, '^', 'MarkerSize', 10); text(3, 3, 'gNB');
plot(trajUE1t(Nd,1,alphInd), trajUE1t(Nd,2,alphInd), 'o', 'MarkerSize',10);
text(trajUE1t(Nd,1,alphInd)-10, trajUE1t(Nd,2,alphInd), 'UE_1');
plot(trajUE2t(Nd,1,alphInd), trajUE2t(Nd,2,alphInd), 'o', 'MarkerSize',10)
text(trajUE2t(Nd,1,alphInd)+5, trajUE2t(Nd,2,alphInd), 'UE_2');
text(trajUE1t(Nd,1,alphInd)/1.3-5,...
    trajUE1t(Nd,2,alphInd)/1.3, 'd', 'FontSize', 16);
text(trajUE2t(Nd,1,alphInd)/1.3+5,...
    trajUE2t(Nd,2,alphInd)/1.3, 'd', 'FontSize', 16);
alphArc = 90-alpha(alphInd):90;
xArc = cosd(alphArc)*d(end-10);
yArc = sind(alphArc)*d(end-10);
plot(xArc, yArc); text(mean(xArc), max(yArc), '\alpha', 'FontSize', 16);

% plot two rays for eNB AA patterb, directed to UE
azAng = -90:90; % plot AA pattern in azimuth (only front hemisphere)
gg1 = zeros(size(azAng));
gg2 = zeros(size(azAng));
gg = zeros(size(azAng));
% vector of eNB AA beamsteering coefficients to direction of UE1
w1 = getAntPatternSteer(antElPos, f, -alpha(alphInd), 0);
% vector of eNB AA beamsteering coefficients to direction of UE2;
% UE2 mirrors UE1
w2 = getAntPatternSteer(antElPos, f,  alpha(alphInd), 0);
% calculate eNB AA pattern for UE1 (gg1) and UE2 (gg2), 
% accounting steering beamsteering vectors w1 and w2
for i=1:length(gg1)
    gg1(i) = getAntPatternG(antElPos, ...
        f, azAng(i), 0, w1, backLobe)/size(antElPos,1);
    gg2(i) = getAntPatternG(antElPos, ...
        f, azAng(i), 0, w2, backLobe)/size(antElPos,1);
end
% calculate coordinates of AA pattern for visulization 
alphPatt = azAng + 90; % rotated 90 degrees for better visualization
xPatt1 = cosd(alphPatt).*gg1*mean(d);
yPatt1 = sind(alphPatt).*gg1*mean(d);
xPatt2 = cosd(alphPatt).*gg2*mean(d);
yPatt2 = sind(alphPatt).*gg2*mean(d);
set(gca,'ColorOrderIndex',1); plot(xPatt1, yPatt1);
plot(xPatt2, yPatt2); xlabel('x, m'); ylabel('y, m'); grid on; axis equal;