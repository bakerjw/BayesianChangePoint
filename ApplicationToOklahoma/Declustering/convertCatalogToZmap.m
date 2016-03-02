function zmap = convertCatalogToZmap(catalog)
%This function converts the catalog as prepared for Oklahoma events to a
%format that is required to be used in Reasenberg Declustering algorithm.

zmap = [catalog.longitude];
zmap(2, :) = [catalog.latitude];
% Decimal year
zmap(3, :) = decyear([catalog.originYear], [catalog.originMonth], [catalog.originDay],...
    [catalog.originHour], [catalog.originMinute], [catalog.originSec]);
zmap(4, :) = [catalog.originMonth];
zmap(5, :) = [catalog.originDay];
zmap(6, :) = [catalog.prefMag];
zmap(7, :) = [catalog.depth];
zmap(8, :) = [catalog.originHour];
zmap(9, :) = [catalog.originMinute];
zmap(10, :) = [catalog.originSec];

zmap = zmap';