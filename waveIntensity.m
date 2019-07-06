function p = waveIntensity(x,amp,phase,relFreqConst,numSources,origins)
% WaveIntensity Intensity function for opticalInterferenceDemo.

%   Copyright 2009 The MathWorks, Inc.  

d = distanceFromSource(x,numSources,origins);
ampVec = [sum(amp./d .* cos(phase - d*relFreqConst));
    sum(amp./d .* sin(phase - d*relFreqConst))];

% Intensity is ||AmpVec||^2
p = ampVec'*ampVec;