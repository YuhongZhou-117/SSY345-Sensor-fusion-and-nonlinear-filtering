close all
clear
clc

%% Problem 1
% Plot colors
color1 = [0.4660 0.6740 0.1880];
color2 = [0.9290 0.6940 0.1250];
color3 = [0.8500 0.3250 0.0980];
color4 = [0 0.4470 0.7410];
color5 = [0.4940 0.1840 0.5560];
color6 = [0.6350 0.0780 0.1840];

% % True track
% Sampling period
T = 0.1;
% Length of time sequence
K = 600;
% Allocate memory
omega = zeros(1,K+1);
% Turn rate
omega(150:450) = -pi/301/T;
% Initial state
x0 = [0 0 20 0 omega(1)]';
% Allocate memory
X = zeros(length(x0),K+1);
X(:,1) = x0;
% Create true track
for i=2:K+1
    % Simulate
    X(:,i) = coordinatedTurnMotion(X(:,i-1), T);
    % Set turn−rate
    X(5,i) = omega(i);
end

% The initial prior
x_0 = [0 0 0 0 0]';
P_0 = diag([10^2 10^2 10^2 (5*pi/180)^2 (pi/180)^2]);

% The sensor position
s1 = [300 -100]';
s2 = [300 -300]';

% The measurement noise standard deviations
Sigma_phi1 = pi/180;
Sigma_phi2 = pi/180;
R = diag([Sigma_phi1^2, Sigma_phi2^2]);

%% 1a)
% Non-linear Kalman filter

% function
f = @(x) coordinatedTurnMotion(x, T);
h = @(x) dualBearingMeasurement(x, s1, s2);

Sigma_v = 1*0.1;
Sigma_w = pi/180*0.16;
Q = diag([0 0 Sigma_v 0 Sigma_w]).^2;

% measurement sequences
Y = genNonLinearMeasurementSequence(X, h, R);

% Get the position from the measurements
pos = zeros(600,2); 
for i = 1:K
    [x_pos, y_pos] = getPosFromMeasurement(Y(1,i), Y(2,i), s1, s2);
    pos(i,:) = [x_pos, y_pos];
end

% filter
[xf_f,Pf_f,~,~] = nonLinearKalmanFilter(Y, x_0, P_0, f, Q, h, R, 'CKF');

% Non-linear smoother
genSigmaPoints = @sigmaPoints;
S = zeros(2,length(Y));

% functions for non-linear smoother
fs = @(x,T) coordinatedTurnMotion(x, T);
hs = @(x,T) dualBearingMeasurement(x,s1,s2);

% smoother
[xs, Ps, ~, ~, ~, ~] = nonLinRTSsmoother(Y, x_0, P_0, fs, T, Q, S, hs, R, genSigmaPoints, 'CKF');


% filter
figure()
grid on 
hold on
% Sensor position
scatter(s1(1),s1(2),150,'o','LineWidth',2, 'DisplayName', 'Sensor1');
scatter(s2(1),s2(2),150,'+','LineWidth',2, 'DisplayName', 'Sensor2');
% True positions
plot(X(1,:), X(2,:), 'color', color1, 'LineWidth', 1.5, 'DisplayName', 'True Position');
% estimated positions
plot(xf_f(1,:), xf_f(2,:), 'color', color2, 'LineWidth', 1.5, 'DisplayName', 'Estimated positions');
% measurement positions
plot(pos(1:end,1), pos(1:end,2), '--', 'color', color3, 'LineWidth', 1.5, 'DisplayName', 'Measurements');
% 3σ-contours at every 5th estimate
for i= 1:5:length(xf_f)
        xy = sigmaEllipse2D(xf_f(1:2,i), Pf_f(1:2,1:2,i), 3, 50);
        plot(xy(1,:),xy(2,:),'color',color4,'HandleVisibility', 'off');
end
% Add legend for 3σ-contours at every 5th estimate
plot(NaN,NaN,'color',color4,'DisplayName','$3\sigma$-contours at 5th');

legend('show', 'Interpreter', 'latex');

% Smoother
figure()
grid on 
hold on
% Sensor position
scatter(s1(1),s1(2),150,'o','LineWidth',2, 'DisplayName', 'Sensor1');
scatter(s2(1),s2(2),150,'+','LineWidth',2, 'DisplayName', 'Sensor2');
% True positions
plot(X(1,:), X(2,:), 'color', color1, 'LineWidth', 1.5, 'DisplayName', 'True Position');
% smooth positions
plot(xs(1,:), xs(2,:), 'color', color5, 'LineWidth', 1.5, 'DisplayName', 'Smoothed positions');
% measurement positions
plot(pos(1:end,1), pos(1:end,2), '--', 'color', color3, 'LineWidth', 1.5, 'DisplayName', 'Measurements');
% 3σ-contours at every 5th estimate
for i= 1:5:length(xs)
        xy = sigmaEllipse2D(xs(1:2,i), Ps(1:2,1:2,i), 3, 50);
        plot(xy(1,:),xy(2,:),'color',color6,'HandleVisibility', 'off');
end
% Add legend for 3σ-contours at every 5th estimate
plot(NaN,NaN,'color',color6,'DisplayName','$3\sigma$-contours at 5th');

legend('show', 'Interpreter', 'latex');

% Zoom in
figure()
grid on 
hold on

% Sensor position
scatter(s1(1),s1(2),150,'o','LineWidth',2, 'DisplayName', 'Sensor1');
scatter(s2(1),s2(2),150,'+','LineWidth',2, 'DisplayName', 'Sensor2');

% True positions
plot(X(1,:), X(2,:), 'color', color1, 'LineWidth', 1.5, 'DisplayName', 'True Position');

% Measurement positions
plot(pos(1:end,1), pos(1:end,2), '--', 'color', color3, 'LineWidth', 1.5, 'DisplayName', 'Measurements');

% Estimated positions
plot(xf_f(1,:), xf_f(2,:), 'color', color2, 'LineWidth', 1.5, 'DisplayName', 'Estimated positions');

% 3σ-contours at every 5th estimate
for i= 1:5:length(xf_f)
    xy = sigmaEllipse2D(xf_f(1:2,i), Pf_f(1:2,1:2,i), 3, 50);
    plot(xy(1,:), xy(2,:), 'color', color4,'HandleVisibility', 'off');
end

% Smooth positions
plot(xs(1,:), xs(2,:), 'color', color5, 'LineWidth', 1.5, 'DisplayName', 'Smoothed positions');

% 3σ-contours at every 5th estimate
for i= 1:5:length(xs)
    xy = sigmaEllipse2D(xs(1:2,i), Ps(1:2,1:2,i), 3, 50);
    plot(xy(1,:), xy(2,:), 'color', color6,'HandleVisibility', 'off');
end
% Add legend for 3σ-contours at every 5th estimate
plot(NaN,NaN,'color',color4,'DisplayName','$3\sigma$-contours at 5th');
plot(NaN,NaN,'color',color6,'DisplayName','$3\sigma$-contours at 5th');

legend('show', 'Interpreter', 'latex');


%% b)
% Modify a measurement at one time instance k = 300
Y_1b = Y;
Y_1b(:,300) = Y_1b(:,300) + 20*rand();

% Get the position from the measurements
pos_1b = zeros(600,2); 
for i = 1:K
    [x_pos, y_pos] = getPosFromMeasurement(Y_1b(1,i), Y_1b(2,i), s1, s2);
    pos_1b(i,:) = [x_pos, y_pos];
end

% filter
[xf_f_1b,Pf_f_1b,~,~] = nonLinearKalmanFilter(Y_1b, x_0, P_0, f, Q, h, R, 'CKF');

% smoother
[xs_1b, Ps_1b, ~, ~, ~, ~] = nonLinRTSsmoother(Y, x_0, P_0, fs, T, Q, S, hs, R, genSigmaPoints, 'CKF');

% Zoom in
figure()
grid on 
hold on

% Sensor position
scatter(s1(1),s1(2),150,'o','LineWidth',2, 'DisplayName', 'Sensor1');
scatter(s2(1),s2(2),150,'+','LineWidth',2, 'DisplayName', 'Sensor2');

% True positions
plot(X(1,:), X(2,:), 'color', color1, 'LineWidth', 1.5, 'DisplayName', 'True Position');

% Measurement positions
plot(pos_1b(1:end,1), pos_1b(1:end,2), '--', 'color', color3, 'LineWidth', 1.5, 'DisplayName', 'Measurements');

% Estimated positions
plot(xf_f_1b(1,:), xf_f_1b(2,:), 'color', color2, 'LineWidth', 1.5, 'DisplayName', 'Estimated positions');

% 3σ-contours at every 5th estimate
for i= 1:5:length(xf_f_1b)
    xy = sigmaEllipse2D(xf_f_1b(1:2,i), Pf_f_1b(1:2,1:2,i), 3, 50);
    plot(xy(1,:), xy(2,:), 'color', color4,'HandleVisibility', 'off');
end

% Smooth positions
plot(xs_1b(1,:), xs_1b(2,:), 'color', color5, 'LineWidth', 1.5, 'DisplayName', 'Smoothed positions');

% 3σ-contours at every 5th estimate
for i= 1:5:length(xs_1b)
    xy = sigmaEllipse2D(xs_1b(1:2,i), Ps_1b(1:2,1:2,i), 3, 50);
    plot(xy(1,:), xy(2,:), 'color', color6,'HandleVisibility', 'off');
end
% Add legend for 3σ-contours at every 5th estimate
plot(NaN,NaN,'color',color4,'DisplayName','$3\sigma$-contours at 5th');
plot(NaN,NaN,'color',color6,'DisplayName','$3\sigma$-contours at 5th');

legend('show', 'Interpreter', 'latex');