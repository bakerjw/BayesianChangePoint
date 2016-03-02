function [ lats, longs ] = getStatePolygon( stateName )
%returns the latitude and longitude polygon associated with a state
%   The function uses the mapping toolbox in MATLAB.

%% Enter folder location of shapefiles
FOLDER_PATH = 'Shapefiles';
shapeFile = 'statep010';

%% Read shapefile and store in variable
folderSelf = fileparts(mfilename('fullpath'));  % Get the folder location of the current file
states = shaperead(fullfile(folderSelf, FOLDER_PATH, shapeFile), 'UseGeoCoords', true);

%% Find polygons for state 
ind = find(ismember({states.STATE}, stateName), 1);
if ~isempty(ind)
    lats = states(ind).Lat;
    longs = states(ind).Lon;
end
