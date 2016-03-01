function eventsDate = getEventDates(dataVect, startDate)
% This function returns the dates when events happened based startDate and
% dataVect (which contains the number of days to an event after previous
% event)

timeFromStart = cumsum(dataVect);

eventsDate = zeros(size(dataVect));
for i = 1:length(dataVect)
    eventsDate(i) = addtodate(startDate, timeFromStart(i), 'day');
end