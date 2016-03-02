function [catalog, lastDate] = catalogFilterRect(catalog,...
                                    startDate, endDate, magMin, magMax, localVect)
% This function filters the catalog for events between startDate and
% endDate, magnitudes between magMin and MagMax, and in a local region
% defined by localVect - [siteLat, siteLong, deltaLat, deltaLong]. localVect is optional.
% deltaLat and deltaLong define the +- ranges over which data is gathered.
% The function also returns the last date recorded in the catalog, prior to
% filtering.

if nargin <= 5 || isempty(localVect)
    localFlag = 0;  % function is same as catalogFilter
else
    localFlag = 1;
    siteLat = localVect(1);
    siteLong = localVect(2);
    dLat = localVect(3);
    dLong = localVect(4);
end

%% Filter the catalog over dates and magnitude range
[catalog, lastDate] = catalogFilter(catalog, startDate, endDate, magMin, magMax);

%% Determine the local region over which to count events
if localFlag
    latVect = [catalog.latitude];
    lonVect = [catalog.longitude];
    indx = (abs(latVect - siteLat) <= dLat) & (abs(lonVect - siteLong) <= dLong);
    catalog = catalog(indx);
end