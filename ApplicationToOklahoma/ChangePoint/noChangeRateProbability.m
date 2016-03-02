function [probVect] = noChangeRateProbability(dataVect, rateVect, params) 
%noChangeRateProbability calculates the probability distribution of the
%event rate given the datavect assuming that no change in rates has
%occurred. 
%   noChangeRateProbability calculates the probability distribtuion of the
%   event rate simply based on gamma distribution as the conjugate
%   prior. In this methodology, the parameters of the gamma posterior are
%   (k+n, 1/(1/theta + sum(x)). dataVect is the
%   vector with the time between each successive event. params is a 1x2 
%   matrix where the row defines the shape and scale of the event rate
%   in a gamma distribution. 

kPrior = params(1, 1);
thetaPrior = params(1, 2);

%% Calculate probability of event rate 
kPosterior = kPrior + length(dataVect);
thetaPosterior = 1/(1/thetaPrior + sum(dataVect));

probVect = gampdf(rateVect, kPosterior, thetaPosterior);

