function [probVect] = changePointRateProbability(dataVect, rateVect, params, startDate, flag) 
%changePointRateProbability calculates the probability distribution of the
%event rate before or after the change point occured using Raftery Akman (1986) approach.
%Pre or Post rate is defined using flag. flag = 'pre' for pre rate and
%'post' for post-rate, 'prePostRatio' or 'postPreRatio'.
%rateVect is the vector of rates over which probabilities will
%be caculated. The ratios calculation was found to have an error in the
%original paper. The equation used here is described in 
%Gupta, A., and Baker, J. W. (2015). “A Bayesian change point model to detect 
%changes in event occurrence rates, with application to induced seismicity.” 
%12th International Conference on Applications of Statistics and Probability in Civil Engineering (ICASP12).

%   changePointRateProbability calculates the probability distribtuion of the
%   pre/post-change point event rate over the time range specified in dataVect. 
%   dataVect is the vector with the time between each successive event. 
%   params is a 2x2 matrix where each row defines the shape and scale of the event rate
%   before and after the change point in a gamma distribution. The event
%   rates are assumed with gamma distribtuion priors and the change point
%   is assumed with a continuous uniform distribtuion over the total
%   duration of data.

if strcmpi(flag, 'pre')
    flag = 0;
elseif strcmpi(flag, 'post')
    flag = 1;
elseif strcmpi(flag, 'prePostRatio')
    flag = 2;
elseif strcmpi(flag, 'postPreRatio')
    flag = 3;
else
    error('flag is invalid');
end

k1 = params(1, 1);
theta1 = params(1, 2);
k2 = params(2, 1);
theta2 = params(2,2);

%% Calculate probability of event rate before the change point occured
[datesVect, totDays] = getDatesVect(dataVect, startDate);
[numEventsVect, totEvents] = getNumEventsVect(dataVect, totDays);

probVect = zeros(size(rateVect));
probVectUnscaled = zeros(size(rateVect));
logScaleVect = zeros(size(rateVect));   % This vector stores the mean of log values at each rate
for i = 1:length(rateVect)
    xVal = rateVect(i);
    tau = zeros(size(datesVect));
    for j = 1:length(datesVect)
        tau(j) = daysact(startDate, datesVect(j));
    end
    numEvents = numEventsVect;
    r1 = numEvents + k1;
    r2 = totEvents - numEvents + k2;
    s1 = tau + 1/theta1;
    s2 = totDays - tau + 1/theta2;
    switch flag
        case 0 %Pre-rate
            probLogVect = -log(totDays) + (r1 - 1).*log(xVal) - xVal.*s1 + gammaln(r2) - r2.*log(s2);
        case 1 %Post-rate
            probLogVect = -log(totDays) + (r2 - 1).*log(xVal) - xVal.*s2 + gammaln(r1) - r1.*log(s1);
        case 2 %Pre-Post ratio
%                 probLogVect(j) = -log(totDays) + gammaln(r1) + gammaln(r2) + (r1 - 1)*log(xVal)...
%                     + r1*log(s2/s1) - (r1 + r2)*log(s2 + xVal*s1);
            probLogVect = -log(totDays) + (r1 - 1).*log(xVal)...
                - (r1 + r2).*log(s2 + xVal.*s1);
        case 3  %Post-pre ratio
%                 probLogVect(j) = -log(totDays) + gammaln(r1) + gammaln(r2) + (r2 - 1)*log(xVal)...
%                     + r2*log(s1/s2) - (r1 + r2)*log(s1 + xVal*s2);
            probLogVect = -log(totDays) + (r2 - 1).*log(xVal)...
                - (r1 + r2).*log(s1 + xVal.*s2);
    end
    [probExpVect, logScaleVect(i)] = scaleExponentials(probLogVect);    % Exponentiate probability vector to account for very large or very small values.
    probVectUnscaled(i) = sum(probExpVect);
end

%% Scale unscaled probability vector such that each value is scaled to the
% same amount, and can be used to calculate probabilities

scale = mean(logScaleVect);
scaleRatio = 10;    %Ratio by which scale is reduced if sum of scaled probabilities is infinity
while true
    probVect = probVectUnscaled.*exp(logScaleVect - scale);
    if ~isinf(trapz(rateVect, probVect))
        break;
    else
        scale = scale + log(scaleRatio);
    end
end

%% Scale probVect to obtain pdf
% The probVect obtained up to this point is a represenatation of the pdf of
% the posterior distribution which is unscaled. To scale it, we make an
% approximation by integrating this pdf at the discrete values provided in
% rateVect, and making the resulting integrand equal to 1 by proper
% scaling.
integrand = trapz(rateVect, probVect);
probVect = probVect./integrand;  %Normalize probability vector
