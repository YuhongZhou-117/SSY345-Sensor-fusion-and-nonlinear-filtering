function [x, P] = tu_qw(x, P, omega, T, Rw)
% the time update function
% omega    the measured angular rate
% T        the time since the last measurement
% Rw       the process noise covariance matrix

% expressions for F and G
F = eye(size(x, 1)) + (T/2)* Somega(omega);
G = (T/2)*Sq(x);

% update
x = F*x;
P = F*P*F' + G*Rw*G';

% normalizes the quaternion
[x, P] = mu_normalizeQ(x, P);
end
