function [x, P,vk] = linearUpdateVk(x, P, y, H, R)
%LINEARPREDICTION calculates mean and covariance of predicted state
%   density using a linear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   y           [m x 1] Measurement
%   H           [m x n] Measurement model matrix
%   R           [m x m] Measurement noise covariance
%
%Output:
%   x           [n x 1] updated state mean
%   P           [n x n] updated state covariance
%

% Your code here
% Kalman filter gain
K = P*H'*inv(H*P*H' + R);

% innovation
vk= y - H*x;

% Update
x = x + K*vk;
P = P - K*H*P;
end