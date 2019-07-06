function [c,ceq] = apertureConstraint(x,xcoords,ycoords)
% apertureConstraint Aperture constraint function for opticalInterferenceDemo.

%   Copyright 2009 The MathWorks, Inc.  

ceq = []; 
c = (x(1) - xcoords).^2 + (x(2) - ycoords).^2 - 9;