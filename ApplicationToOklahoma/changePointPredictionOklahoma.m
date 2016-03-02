% This file performs the change point calculations as described in the
% in-review paper - Gupta and Baker, Estimating Spatially Varying Event 
% Rates with a Change Point using Bayesian Statistics: Application to Induced Seismicity

% This script divides the state into equally spaced grid points and
% implements change point at each grid point. 
% There is a stop date defined. The stop date is for the end of training
% catalog. The end date is for the end of test catalog. 
% It calculates the probability gain 
% using the methodology specified in Werner et al, 2011, with reference to 
% a uniform rate model.

%% Specify input parameters
stateName = {'Oklahoma'}; % State name used to draw state on map

fileIn = 'complete.csv';    % earthquake catalog file
magMin = 3;
magMax = 8;

% Oklahoma bounding box
latRange = [33.5 37];
longRange = [-103 -94.5];

gridStep = 0.1; % grid over which rates are calculated
radkm = 25;    % Radius of circular region in km

% Parameters for no change model
k0 = 0.5;
theta0 = Inf;
% Parameters for single change model
k1 = 0.5;
k2 = 0.5;
theta1 = Inf;
theta2 = Inf;

startDate = datenum(1974, 1, 1); % January 01, 1974 for the data
stopDate = datenum(2014, 12, 31); % stop date of training catalog (up to which rates are estimated)
endDate = datenum(2015, 06, 30);  % end date of test catalog

%Critical bayesFactor
critBF = 0.001;

% Critical number of events required to detect change
critEQ = 2; % If a region has less than critEQ events, no processing is done

%% Add paths to required folders
currentPath = fileparts(mfilename('fullpath'));
addpath(fullfile(currentPath, 'Declustering'));
addpath(fullfile(currentPath, 'ChangePoint'));
addpath(fullfile(currentPath, 'CatalogFunctions'));
addpath(fullfile(currentPath, 'Mapping'));

%% Read input file into structure
catalog = importOKCatalogfile(fileIn);

%% Process inputs in function input format
paramsCP = [k1, theta1; % Parameters for change point model
        k2, theta2];    
paramsNC = [k0, theta0];    %Parameters for no change model used in Bayes Factor calculation

%% Create grid of points
latVect = latRange(1):gridStep:latRange(2);
longVect = longRange(1):gridStep:longRange(2);
[latGrid, longGrid] = meshgrid(latVect, longVect);

%% Determine which points are within the state boundary
[stateLats, stateLongs] = getStatePolygon(stateName);
inpoints = inpolygon(latGrid, longGrid, stateLats, stateLongs);

%% Filter catalog for dates and magnitude ranges
catalogDC = declusterCatalog(catalog); % Decluster complete catalog
catalogDC = catalogFilter(catalogDC, startDate, endDate, magMin, magMax);

% Decluster training catalog separately, since we are assuming that we have
% not seen any events post the training catalog
catalogTraining = catalogFilter(catalog, startDate, stopDate, -2, 10);  % Avoid magnitude filtering at this point
catalogTraining = declusterCatalog(catalogTraining);    % Decluster training catalog
catalogTraining = catalogFilter(catalogTraining, startDate, stopDate, magMin, magMax);  % magnitude filtering
totTrainEvents = length(catalogTraining);
numTrainDays = daysact(startDate, stopDate);

% Calculate number of events in test catalog
[catalogTest, lastDate] = catalogFilter(catalogDC, stopDate + 1, endDate, magMin, magMax);
totTestEvents = length(catalogTest);
numTestDays = daysact(stopDate + 1, lastDate);

%% Run change point detection on all grid points
parobj = parpool(4);    % Parallel processing on Matlab
changeDates = NaN(size(latGrid));   % NaN changes to date if change is detected
trainingRates = changeDates;    % NaN changes to rate if gridpoint is within the state
numTrainEvents = changeDates;
numTestEvents = changeDates;
parfor i = 1:size(latGrid, 1)
    % Additional variables required for parfor
    changeDatesVect = changeDates(i, :);
    trainingChangeRateVect = changeDatesVect;
    probObservedVect = changeDatesVect;
    numTrainEventsVect = changeDatesVect;
    numTestEventsVect = changeDatesVect;
    
    for j = 1:size(latGrid, 2)
        % Continue only if this point is within state boundary
        if (inpoints(i, j) == 0), continue; end     % pass control to next iteration
        trainingRate = 0;   %training rate for gridpoint(i, j)
        % Get earthquake catalog within circular region
        localVect = [latGrid(i, j), longGrid(i, j), radkm];   % Circular region
        catalogTrainingLocal = catalogFilter(catalogTraining, startDate, stopDate, magMin, magMax, localVect);
        
        % Obtain timevect for inter-event times
        if ~isempty(catalogTrainingLocal)
            [dataVect, startDateCP] = getTimeBetweenEvents(catalogTrainingLocal);
            % Change data vector to include information that no events occured between
            % startdate and startdateCP (date of first event in local calaog)
            numDaysNoEvent = daysact(startDate, startDateCP);
            if numDaysNoEvent > 0
                dataVect(2:(end + 1)) = dataVect;
                dataVect(1) = numDaysNoEvent;
            end
                
            % Perform change-point analysis only if events in local region exceed
            % critEQ
            if length(catalogTrainingLocal) >= critEQ
                % Calculate bayes factor and change point for stop date data set
                bFCP = changePointBayesFactorGeneralized(dataVect, paramsCP, paramsNC, startDate);
                if bFCP <= critBF % Implement change point model
                    [datesVect, probVect] = changePointProbability(dataVect, paramsCP, startDate);
                    [maxProb, maxIndx] = max(probVect);
                    dateCP = datesVect(maxIndx);

                    % Calculate probability of post-change point event rate
                    % till stop date
                    trainingRateVect =  logspace(-10, 0, 201);
                    trainingRateProbVect = changePointRateProbability(dataVect,...
                        trainingRateVect, paramsCP, startDate, 'post');
                    [~, maxIndx] = max(trainingRateProbVect);
                    trainingRate = trainingRateVect(maxIndx);

                    changeDatesVect(j) = dateCP;
                end
            end

            if trainingRate == 0 % if no change is detected, compute rate for no change case
                trainingRateVect =  logspace(-10, 0, 201);
                trainingRateProbVect = noChangeRateProbability(dataVect,...
                    trainingRateVect, paramsNC);
                [~, maxIndx] = max(trainingRateProbVect);
                trainingRate = trainingRateVect(maxIndx);   
            end
            
            % Normalize rate to per sq. km (it is already per day)
            trainingRate = trainingRate/(pi*radkm*radkm); 
        end 

        trainingChangeRateVect(j) = trainingRate;
                        
        % Calculate number of events in test catalog
        localRect = [latGrid(i, j), longGrid(i, j), gridStep/2, gridStep/2];
        [catalogTestLocal] = catalogFilterRect(catalogDC, stopDate + 1, endDate, magMin, magMax, localRect);
        regionArea = calculateMapAreaRect(localRect);
        numTestEventsVect(j) = length(catalogTestLocal);
        numTrainEventsVect(j) = length(catalogTrainingLocal);
    end
    changeDates(i, :) = changeDatesVect;
    trainingRates(i, :) = trainingChangeRateVect;
    numTrainEvents(i, :) = numTrainEventsVect;
    numTestEvents(i, :) = numTestEventsVect;
end

delete(parobj);

%% Generate probability gain wrt uniform rate model
modelRates = trainingRates(:);
modelRates = modelRates.*(100*numTestDays); % 100 sq. km is the approximate area of each grid cell
unifModelRates = modelRates;
numCells = sum(~isnan(unifModelRates));
unifModelRates(~isnan(unifModelRates)) = (totTrainEvents/numCells)*(numTestDays/numTrainDays);

% Calculate log-likelihoods
testEventsVect = numTestEvents(:);
LLmodel = nansum(testEventsVect.*log(modelRates)) - nansum(modelRates);
LLunifModel = nansum(testEventsVect.*log(unifModelRates)) - nansum(unifModelRates);
probGain = exp((LLmodel - LLunifModel)/totTestEvents);

%% Plot dates of change on a map
% Draw state 
figH = figure('position', [1600/4, 400, 1200, 600]);
axHandle = drawStates(stateName, figH);   % Draws state

% Draw surface of change dates on state
axes(axHandle);
colormap(jet);
surfH = surfm(latGrid, longGrid, changeDates);

hcb = colorbar;
datetick(hcb, 'y');


%% Plot estimated rates on a map
% Draw state 
figH = figure('position', [1600/4, 400, 1200, 600]);
axHandle = drawStates(stateName, figH);   % Draws state

% Draw surface of change dates on state
axes(axHandle);
colormap(jet);
surfH = surfm(latGrid, longGrid, trainingRates);    % Rates are per sq. km per day

hcb = colorbar;
