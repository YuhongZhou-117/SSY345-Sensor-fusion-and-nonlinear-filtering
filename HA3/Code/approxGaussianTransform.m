function [mu_y, Sigma_y, y_s] = approxGaussianTransform(mu_x, Sigma_x, f, N)
%approxGaussianTransform takes a Gaussian density and a transformation 
%function and calculates the mean and covariance of the transformed density.
%
%Inputs
%   MU_X        [m x 1] Expected value of x.
%   SIGMA_X     [m x m] Covariance of x.
%   F           [Function handle] Function which maps a [m x 1] dimensional
%               vector into another vector of size [n x 1].
%   N           Number of samples to draw. Default = 5000.
%
%Output
%   MU_Y        [n x 1] Approximated mean of y.
%   SIGMA_Y     [n x n] Approximated covariance of y.
%   ys          [n x N] Samples propagated through f


    if nargin < 4
        N = 5000;
    end

%Your code here
    % draw N samples from the Gaussian density with parameters given in the input parameters 
    x_s = mvnrnd(mu_x, Sigma_x, N)';
    
    % Samples propagated through f
    y_s = f(x_s);
    
    % Approximated mean of y.
    mu_y = mean(y_s,2);
    
    % Approximated covariance of y.
    Sigma_y = (1/(N-1))*(y_s - mu_y)*(y_s - mu_y)';
    
end
