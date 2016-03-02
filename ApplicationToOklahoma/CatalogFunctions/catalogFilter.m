function [catalog, lastDate] = catalogFilter(catalog,...
                                    startDate, endDate, magMin, magMax, localVect)
% This function filters the catalog for events between startDate and
% endDate, magnitudes between magMin and MagMax, and in a local region
% defined by localVect - [siteLat, siteLong, radkm]. localVect is optional.
% The function also returns the last date recorded in the catalog, prior to
% filtering.
% Catalog is a struct obtained from importOKCatalogfile.m

if nargin <= 5 || isempty(localVect)
    localFlag = 0;
else
    localFlag = 1;
    siteLat = localVect(1);
    siteLong = localVect(2);
    radkm = localVect(3);
end

if length(catalog) < 1
    lastDate = NaN;
    return;
end

%% Determine the time range over which to count events
catalogDates = datenum([catalog.originYear], [catalog.originMonth], [catalog.originDay]);
lastDate = catalogDates(end);
indx = catalogDates >= startDate & catalogDates <= endDate;
catalog = catalog(indx);    %Reduce catalog to catalog of interest

%% Determine the magnitude range over which to count events
catalogMags = [catalog.prefMag];
indx = catalogMags >= magMin & catalogMags <= magMax;
catalog = catalog(indx);

%% Determine the local region over which to count events
if localFlag
    latVect = [catalog.latitude];
    lonVect = [catalog.longitude];
    siteLatVect = siteLat.*ones(size(latVect));
    siteLonVect = siteLong.*ones(size(lonVect));
    distanceVect = deg2km(distance(siteLatVect, siteLonVect, latVect, lonVect));
    indx = distanceVect <= radkm;
    catalog = catalog(indx);
end