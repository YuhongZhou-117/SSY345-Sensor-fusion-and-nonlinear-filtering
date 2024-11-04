function [x, P] = mu_m(x, P, mag, m0, Rm)
% the EKF update using ykm
% under the assumption that fkm = 0.
% mag       is shorthand for ykm, the measurement
% m0        is m0
% Rm        the measurement noise covariance matrix

% measurement function
h = Qq(x)'*m0;

% The derivative of Q(q) wrt qi, i={0,1,2,3}
[Q0, Q1, Q2, Q3] = dQqdq(x);

H = [Q0'*m0, Q1'*m0, Q2'*m0, Q3'*m0];

% the update step
% Innovation
v = mag - h;

% Innovation Covariance
S = H*P*H' + Rm;

% Kalman Gain 
K = P*H'*inv(S);

% update 
x = x + K*v;

P = P - K*S*K';

% normalizes the quaternion
[x, P] = mu_normalizeQ(x, P);
end