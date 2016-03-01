function [numEventsVect, totEvents] = getNumEventsVect(dataVect, totDays)
%getDatesVect returns a number of events vector containing number of events
%that have occured till every day in the time range
%specified in dataVect and starting on startDate to be used in changePoint
%calculations

% timeVect = cumsum(dataVect);
% timeVect(end + 1) = 0;
% timeVect(2:end) = timeVect(1:(end - 1));
% timeVect(1) = 0;

totEvents = length(dataVect) + 1;   % 1 is added to account for first event on day 0

numEventsVect = ones(totDays - 1, 1);
% numEvents = 1;  %Number of events till time tau. First event occured on day 0

%% This code block modified per Jack Baker's streamlined code
numEventsVect(1:dataVect(1)-1) = 1; % first event is actally at day = 0, so end index has a '-1'
startDay = dataVect(1); % first day for the next set of assignments

for j = 2:length(dataVect)
    numEventsVect(startDay:startDay+dataVect(j)) = j; 
    startDay = startDay+dataVect(j);
end

numEventsVect(end) = numEventsVect(end) + 1; % one event happened on the last day

%% Legacy Code
% for tau = 1:(totDays - 1)   %tau is the change point starting the day after startDate
%     %Find number of events till time tau
%     ind1 = find(timeVect < tau, 1, 'last');
%     if isempty(ind1), ind1 =  1; end
%     ind2 = find(timeVect > tau, 1, 'first');
%     if isempty(ind2), ind2 = length(timeVect); end
%     tauRange = timeVect(ind1:ind2);
%     numEvents = numEvents + length(find(tauRange == tau));
%     numEventsVect(tau) = numEvents;
% end