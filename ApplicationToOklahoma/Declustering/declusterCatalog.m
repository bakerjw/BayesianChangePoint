function [catalog, indx] = declusterCatalog(catalog)
% This function declusters the input catalog based on Reasenberg (1985) and
% using the algorithm implemented by A. Allman.

currentPath = fileparts(mfilename('fullpath'));
addpath(fullfile(currentPath, 'ReasenbergZMAPversion'));

% Convert catalog to a form that can be read by declustering algortihm
% function
zmapCat = convertCatalogToZmap(catalog);

% Define parameters for declustering algorithm
% rfact  is factor for interaction radius for dependent events (default 10)
% xmeff  is "effective" lower magnitude cutoff for catalog,it is raised 
%         by a factor xk*cmag1 during clusters (default 1.5)
% xk     is the factor used in xmeff    (default .5)
% taumin is look ahead time for not clustered events (default one day)
% taumax is maximum look ahead time for clustered events (default 10 days)
% P      to be P confident that you are observing the next event in 
%        the sequence (default is 0.95)
rfact = 10; 
xmeff = 3;  % Minimum magnitude of completion based on Oklahoma data
xk = 0.5;
taumin = 1;
taumax = 10;
P = 0.95;
err = 0;
derr = 0;

% Get index of mainshocks in catalog using Resenberg algorithm
[declus, indx] = ReasenbergDeclus(taumin, taumax, xk, xmeff, P, rfact, err, derr, zmapCat);

% Apparently, the algorithm returns indx that may not be monotonically
% increasing, hence sort index
indx = sort(indx);
% Return declustered catalog containing only the mainshocks
catalog = catalog(indx);