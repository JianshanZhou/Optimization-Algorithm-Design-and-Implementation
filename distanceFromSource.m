function d = distanceFromSource(v,numSources,origins)
% distanceFromSource Distance function for opticalInterferenceDemo.

%   Copyright 2009 The MathWorks, Inc.  

d = zeros(numSources,1);
for k = 1:numSources
    d(k) = norm(origins(k,:) - [v 0]);
end

