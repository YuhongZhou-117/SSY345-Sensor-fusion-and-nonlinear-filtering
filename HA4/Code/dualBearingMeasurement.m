function [hx, Hx] = dualBearingMeasurement(x, s1, s2)
%DUOBEARINGMEASUREMENT calculates the bearings from two sensors, located in 
%s1 and s2, to the position given by the state vector x. Also returns the
%Jacobian of the model at x.
%
%Input:
%   x           [n x 1] State vector, the two first element are 2D position
%   s1          [2 x 1] Sensor position (2D) for sensor 1
%   s2          [2 x 1] Sensor position (2D) for sensor 2
%
%Output:
%   hx          [2 x 1] measurement vector
%   Hx          [2 x n] measurement model Jacobian
%
% NOTE: the measurement model assumes that in the state vector x, the first
% two states are X-position and Y-position.

% Get the 2D position
px = x(1);
py = x(2);

% Calculate the measurements
hx1 = atan2(py - s1(2),px - s1(1));
hx2 = atan2(py - s2(2),px - s2(1));

hx = [hx1;hx2];

% Calculate the Jacobian matrix
Hx = zeros(2,length(x));

dhx1_dx = -(py - s1(2)) / ((px - s1(1))^2 + (py - s1(2))^2);
dhx1_dy = (px - s1(1)) / ((px - s1(1))^2 + (py - s1(2))^2);
    
dhx2_dx = -(py - s2(2)) / ((px - s2(1))^2 + (py - s2(2))^2);
dhx2_dy = (px - s2(1)) / ((px - s2(1))^2 + (py - s2(2))^2);

Hx(1,1) = dhx1_dx;
Hx(1,2) = dhx1_dy;

Hx(2,1) = dhx2_dx;
Hx(2,2) = dhx2_dy;


end