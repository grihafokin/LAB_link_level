% (c) Fokin G.A., Volgushev D.B., SPbSUT, 2022.
% function to create antenna array (AA)
function antPos = createAnt(antType, Nel, d)
% antType - AA type 
% Nel     - number of AA elements
% d       - distance between AA elements, m
% antPos  - array of coordinates [x,y,z] of AA elements, m 
% dimension [Nel x 3] for uniform linear array (ULA) and
%           [Nel^2 x 3] for uniform rectangular/planar array (URA)
switch antType
    case 1 % URA (square)
        % number of AA elements on one side of a square array
        % total number of elements Nel^2
        yCoords = repmat((-(Nel-1)/2:1:(Nel-1)/2).'*d, Nel, 1);
        zCoords = reshape(repmat((-(Nel-1)/2:1:(Nel-1)/2)*d,Nel,1), [], 1);
        zCoords = zCoords(:);
        xCoords = zeros(size(yCoords));
        antPos = [xCoords, yCoords, zCoords];        
    case 2     % ULA 
        % ULA of Nel elements
        antLen = (Nel-1)*d;
        yCoords = (-antLen/2:d:antLen/2).';
        zCoords = zeros(size(yCoords));
        xCoords = zeros(size(yCoords));
        antPos = [xCoords, yCoords, zCoords];
    case 3 % Uniform Circular Array (UCA)
        Nel = Nel*2; % total number of elements 2*Nel
        % calculation of the radius of a UCA through the formula for the 
        % length of the chord of the arc circle; the chord length is equal 
        % to the given distance d between the AA elements 
        r = d/2/sind(360/Nel/2);
        antDph = (0:360/Nel:360 - 360/Nel).';
        xCoords = r*cosd(antDph);
        yCoords = r*sind(antDph);
        zCoords = zeros(size(yCoords));
        antPos = [xCoords, yCoords, zCoords];  
end
end