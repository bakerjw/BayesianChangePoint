function [datesVect, totDays] = getDatesVect(dataVect, startDate)
%getDatesVect returns a dates vector containing every day in the time range
%specified in dataVect and starting on startDate to be used in changePoint
%calculations

timeVect = cumsum(dataVect);

totDays = timeVect(end) + 1;    % 1 is added to account for startDate

datesVect = zeros(totDays - 1, 1);
for tau = 1:(totDays - 1)   %tau is the change point starting the day after startDate
    datesVect(tau) = addtodate(startDate, tau, 'day');
end