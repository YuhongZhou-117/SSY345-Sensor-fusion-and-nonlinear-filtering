function [x, P] = mu_g(x, P, yacc, Ra, g0)
% the EKF update using yka
% under the assumption that fak = 0  
% yacc      yacc is shorthand for yak, the measurement
% Ra        the measurement noise covariance matrix

% measurement function
h = Qq(x)'*g0;

% The derivative of Q(q) wrt qi, i={0,1,2,3}
[Q0, Q1, Q2, Q3] = dQqdq(x);

H = [Q0'*g0, Q1'*g0, Q2'*g0, Q3'*g0];

% the update step
% Innovation
v = yacc - h;

% Innovation Covariance
S = H*P*H' + Ra;

% Kalman Gain 
K = P*H'*inv(S);

% update 
x = x + K*v;

P = P - K*S*K';

% normalizes the quaternion
[x, P] = mu_normalizeQ(x, P);

end