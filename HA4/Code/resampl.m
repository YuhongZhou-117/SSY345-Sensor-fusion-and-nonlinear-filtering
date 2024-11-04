function [Xr, Wr, j] = resampl(X, W)
%RESAMPLE Resample particles and output new particles and weights.
% resampled particles. 
%
%   if old particle vector is x, new particles x_new is computed as x(:,j)
%
% Input:
%   X   [n x N] Particles, each column is a particle.
%   W   [1 x N] Weights, corresponding to the samples
%
% Output:
%   Xr  [n x N] Resampled particles, each corresponding to some particle 
%               from old weights.
%   Wr  [1 x N] New weights for the resampled particles.
%   j   [1 x N] vector of indices refering to vector of old particles

% Your code here!
    N = size(X,2);
    n = size(X,1);
    
    % Getting samples from the uniform distribution, unif[0, 1],
    samples = rand([1 N]);
    
    % Compute cumulative sum of weights
    sum_cumlative = cumsum(W);
    
    for i = 1:N
        
        % Find the random number in the cumulative sum
        index = find(sum_cumlative >= samples(i), 1);
    
        % Resample particle and weight
        Xr(:,i) = X(:,index);
        Wr(i) = 1/N;
        j(i) = index;
    end
end