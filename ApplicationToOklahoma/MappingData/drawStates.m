function [ axHandle ] = drawStates( stateNames, figHandle )
%drawState draws the US state boundaries specified in cell stateNames
%   The function uses the mapping toolbox in MATLAB.

if nargin <= 1
    figHandle = figure;
end

%% Enter folder location of shapefiles
FOLDER_PATH = 'Shapefiles';
shapeFile = 'statep010';
colorStates = [237 237 237]./255;  % All states will be drawn in this color
frameColor = [206 205 210]./255;
boundaryColor = [237 237 237]./255;

%% Read shapefile and store in variable

folderSelf = fileparts(mfilename('fullpath'));  % Get the folder location of the current file
states = shaperead(fullfile(folderSelf, FOLDER_PATH, shapeFile), 'UseGeoCoords', true);

successes = 0;   % Stores the number of states successfully read from the shapefile

%% Find extents of map to nearest whole latitude and longitude
lengthStates = length(stateNames);
bBoxStates = [180, 90; -180, -90];
successStates = cell(lengthStates, 1);
indStates = cell(lengthStates, 1);
for i = 1:lengthStates
    ind = find(ismember({states.STATE}, stateNames{i}));
    if ~isempty(ind)
        successes = successes + 1;
        successStates{successes} = stateNames{i};
        indStates{successes} = ind;
        for j = 1:length(ind)
            bBox = states(ind(j)).BoundingBox;
            if bBoxStates(1, 1) > bBox(1, 1), bBoxStates(1, 1) = bBox(1, 1); end
            if bBoxStates(1, 2) > bBox(1, 2), bBoxStates(1, 2) = bBox(1, 2); end
            if bBoxStates(2, 1) < bBox(2, 1), bBoxStates(2, 1) = bBox(2, 1); end
            if bBoxStates(2, 2) < bBox(2, 2), bBoxStates(2, 2) = bBox(2, 2); end
        end
    end
end

mapBounds = zeros(2, 2);
mapBounds(1, :) = bBoxStates(1, :) - 0.5;    % Rounding to lower .5 longitude and latitude
mapBounds(2, :) = bBoxStates(2, :) + 0.5;

% Draw the Map with Lambert Conic Projection
figure(figHandle);
axHandle = axesm('lambert', 'MapLonLimit', mapBounds(:, 1), ...
    'MapLatLimit', mapBounds(:, 2), 'frame', 'on', 'grid', 'off',...
    'MeridianLabel', 'on', 'MLabelLocation', 2, 'ParallelLabel', 'on',...
    'PLabelLocation', 2, 'FFaceColor', frameColor, 'FontSize', 12,...
    'grid', 'on', 'mlinelocation', 1, 'plinelocation', 1, 'gcolor', [0.7, 0.7, 0.7]);
set(axHandle, 'Visible', 'off')

% Draw each state
for i = 1:successes
    geoshow(axHandle, states(indStates{i}), 'FaceColor', colorStates, 'EdgeColor', boundaryColor);
end

end

