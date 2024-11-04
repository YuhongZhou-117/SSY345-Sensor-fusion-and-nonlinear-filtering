function Y = genNonLinearMeasurementSequence(X, h, R)
%GENNONLINEARMEASUREMENTSEQUENCE generates ovservations of the states 
% sequence X using a non-linear measurement model.
%
%Input:
%   X           [n x N+1] State vector sequence
%   h           Measurement model function handle
%   h           Measurement model function handle
%               [hx,Hx]=h(x) 
%               Takes as input x (state) 
%               Returns hx and Hx, measurement model and Jacobian evaluated at x
%   R           [m x m] Measurement noise covariance
%
%Output:
%   Y           [m x N] Measurement sequence
%

% Dimensions
n = size(X, 1);
N = size(X, 2)-1;
m = length(R);

% Initial setting
Y = zeros(m, N);

for k = 1:N
    % for time instant k 
    [hx,~]=h(X(:,k+1));
    
    Y(:, k) = hx + mvnrnd(zeros(m, 1), R)';
end

end