function [datesVect, probVect] = changePointProbability(dataVect, params, startDate) 
%changePointProbability calculates the probability distribution of change
%point using Raftery Akman (1986) approach.
%   ChangePointProbability calculates the probability distribtuion of the
%   change point over the time range specified in dataVect. dataVect is the
%   vector with the time between each successive event. params is a 2x2 
%   matrix where each row defines the shape and scale of the event rate
%   before and after the change point in a gamma distribution. The event
%   rates are assumed with gamma distribtuion priors and the change point
%   is assumed with a continuous uniform distribtuion over the total
%   duration of data.


k1 = params(1, 1);
theta1 = params(1, 2);
k2 = params(2, 1);
theta2 = params(2,2);

%% Calculate probability of change point given data
[datesVect, totDays] = getDatesVect(dataVect, startDate);
[numEventsVect, totEvents] = getNumEventsVect(dataVect, totDays);

probVect = zeros(totDays - 1, 1);
for i = 1:length(datesVect)
    tau = daysact(startDate, datesVect(i));
    numEvents = numEventsVect(i);
    r1 = numEvents + k1;
    r2 = totEvents - numEvents + k2;
    s1 = tau + 1/theta1;
    s2 = totDays - tau + 1/theta2;
    probVect(i) = -log(totDays) + gammaln(r1) + gammaln(r2) - r1*log(s1) - r2*log(s2);
end

[probVect, scale] = scaleExponentials(probVect);   % Exponentiate probability vector to account for very large or very small values.
probVect = probVect./sum(probVect);  %Normalize to obtain pdf (this is indeed the pdf since days are equally spaced @ 1 day)