function [bayesFactor] = changePointBayesFactorGeneralized(dataVect, paramsCM, paramsNC, startDate) 
%changePointBayesFactor calculates the bayes factor for no change against a
%change point using Raftery Akman (1986) approach.
%   ChangePointBayesFactor calculates the Bayes Factor of the
%   change point over the time range specified in dataVect. dataVect is the
%   vector with the time between each successive event. paramsCM is a 2x2 
%   matrix for parameters of the change model prior rates 
%   where each row defines the shape and scale of the event rate
%   before and after the change point in a gamma distribution. paramsNC is
%   1x2 matrix for parameters of no change model. The event
%   rates are assumed with gamma distribtuion priors and the change point
%   is assumed with a continuous uniform distribtuion over the total
%   duration of data. Smaller than 1 Bayes Factor imply a shift towards
%   change model.


k1 = paramsCM(1, 1);
theta1 = paramsCM(1, 2);
k2 = paramsCM(2, 1);
theta2 = paramsCM(2,2);

k0 = paramsNC(1);
theta0 = paramsNC(2);


%% Obtain date and event vectors from data
[datesVect, totDays] = getDatesVect(dataVect, startDate);
[numEventsVect, totEvents] = getNumEventsVect(dataVect, totDays);


%% Calculate numerater of Bayes Factor
if isinf(theta0)
    numLog = gammaln(totEvents + k0) - gammaln(k0) - (totEvents + k0)*log(1/theta0 + totDays);
else
    numLog = gammaln(totEvents + k0) - gammaln(k0) - k0*log(theta0) - (totEvents + k0)*log(1/theta0 + totDays);
end

%% Calculate denominator of Bayes Factor
tau = zeros(totDays - 1, 1);
for i = 1:length(datesVect)
    tau(i) = daysact(startDate, datesVect(i));
end
numEvents = numEventsVect;
r1 = numEvents + k1;
r2 = totEvents - numEvents + k2;
s1 = tau + 1/theta1;
s2 = totDays - tau + 1/theta2;
denVect = -log(totDays) + gammaln(r1) + gammaln(r2) - r1.*log(s1) - r2.*log(s2);

[denVect, scale] = scaleExponentials(denVect);   % Exponentiate probability vector to account for very large or very small values.
if isinf(theta1)
    if isinf(theta2)
        denLog = log(sum(denVect)) + scale - gammaln(k1) - gammaln(k2);
    else
        denLog = log(sum(denVect)) + scale - k2*log(theta2) - gammaln(k1) - gammaln(k2);
    end
elseif isinf(theta2)
    denLog = log(sum(denVect)) + scale - k1*log(theta1) - gammaln(k1) - gammaln(k2);
else
    denLog = log(sum(denVect)) + scale - k1*log(theta1) - k2*log(theta2) - gammaln(k1) - gammaln(k2);
end
    
bayesFactor = exp(numLog - denLog);

% Use Eq 3.2 in Raftery Akman when certian conditions are met
% if isinf(theta0) && isinf(theta1) && isinf(theta2) && k0 == 0.5 && k1 == 0.5 && k2 == 0.5
%     bayesFactor = bayesFactor*((4*pi)^0.5)*(totDays^-0.5);
% else
%     warning(['The constant in Bayes factor calculation is evaluated only for', char(10),...
%         'the case of gamma priors of 0.5 and infinity. For other prior values,' char (10),...
%         'Bayes factor may only be used for comparison between datasets.']);
% end

% Calculate the proportioanlity factor based on description by Raftery
% Akman
if k1 ~= k2
    warning(['The constant in Bayes factor calculation requires the', ...
        ' parameters k1 and k2 to be equal.', char(10), 'Currently, k1 = ', num2str(k1), ...
        ' and k2 = ', num2str(k2), char(10), 'This Bayes factor may only be used for ', ...
        'comparison between datasets with equal days of observation.']);
else
    % Calculate appropriate constant factor by equating the Bayes factor
    % for the boundary condition to 1. The boundary condition is one event
    % happening at the middle of the time range.
    totEvents = 1;
    if isinf(theta0)
        numLogBoundary = gammaln(totEvents + k0) - gammaln(k0) - (totEvents + k0)*log(1/theta0 + totDays);
    else
        numLogBoundary = gammaln(totEvents + k0) - gammaln(k0) - k0*log(theta0) - (totEvents + k0)*log(1/theta0 + totDays);
    end
    tau = (1:(totDays - 1))';
    numEvents = zeros(size(tau));
    numEvents(ceil(totDays/2):end) = 1;
    r1 = numEvents + k1;
    r2 = totEvents - numEvents + k2;
    s1 = tau + 1/theta1;
    s2 = totDays - tau + 1/theta2;
    denVectBoundary = -log(totDays) + gammaln(r1) + gammaln(r2) - r1.*log(s1) - r2.*log(s2);
    [denVectBoundary, scale] = scaleExponentials(denVectBoundary);   % Exponentiate probability vector to account for very large or very small values.
    if isinf(theta1)
        if isinf(theta2)
            denLogBoundary = log(sum(denVectBoundary)) + scale - gammaln(k1) - gammaln(k2);
        else
            denLogBoundary = log(sum(denVectBoundary)) + scale - k2*log(theta2) - gammaln(k1) - gammaln(k2);
        end
    elseif isinf(theta2)
        denLogBoundary = log(sum(denVectBoundary)) + scale - k1*log(theta1) - gammaln(k1) - gammaln(k2);
    else
        denLogBoundary = log(sum(denVectBoundary)) + scale - k1*log(theta1) - k2*log(theta2) - gammaln(k1) - gammaln(k2);
    end
    bayesFactorBoundary = exp(numLogBoundary - denLogBoundary);
    constFactor = 1/bayesFactorBoundary;
    bayesFactor = constFactor*bayesFactor;
end

