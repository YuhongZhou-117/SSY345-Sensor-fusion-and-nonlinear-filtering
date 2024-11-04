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

%% Problem 2
clear
clc

% Plot colors
color1 = [0.4660 0.6740 0.1880];
color2 = [0.9290 0.6940 0.1250];
color3 = [0.8500 0.3250 0.0980];
color4 = [0 0.4470 0.7410];
color5 = [0.4940 0.1840 0.5560];
color6 = [0.6350 0.0780 0.1840];

% Parameters
Q = 1.5;
R = 3;

% The initial prior
x_0 = 2; 
P_0 = 8;

%% a)
% time steps
k = 30;

% Number of Particles
N = 100;

% functions
A = 1;
H = 1;
f = @(x) A*x;
h = @(x) H*x;

% state sequence and measurement sequence
X = genLinearStateSequence(x_0, P_0, A, Q, k);
Y = genLinearMeasurementSequence(X, H, R);

% Kalman filter
[xf, Pf] = kalmanFilter(Y, x_0, P_0, A, Q, H, R);

% Particle filters without resampling
bResample = 1;
[xfp, Pfp, Xp, Wp] = pfFilter(x_0, P_0, Y, f, Q, h, R, N, ~bResample, [], 0);

% Particle filters with resampling
[xfp_re, Pfp_re, Xp_re, Wp_re] = pfFilter(x_0, P_0, Y, f, Q, h, R, N, bResample, [], 0);

% mean square error (MSE)
disp('KF mse:')
immse(xf,X(:,2:end))
disp('PF without resampling')
immse(xfp,X(:,2:end))
disp('PF with resampling')
immse(xfp_re,X(:,2:end))

% Plot 
figure()
grid on
hold on
plot(0:k, X, 'color', color5, 'LineWidth',2, 'DisplayName','True Trajectory');
plot(1:k, Y, '--', 'color', color4, 'LineWidth',2, 'DisplayName','Measurement sequence');
plot(1:k, xf, 'color', color3, 'LineWidth',2, 'DisplayName','KF estimate');
errorbar(1:k,xf,reshape(Pf,[1,length(Pf)]),'Color',color3, 'LineWidth',1.5, 'HandleVisibility', 'off');
plot(1:k, xfp, 'color', color2, 'LineWidth',2, 'DisplayName','PF without Resample');
errorbar(1:k,xfp,reshape(Pfp,[1,length(Pfp)]),'Color',color2, 'LineWidth',1.5 ,'HandleVisibility', 'off');
plot(1:k, xfp_re, 'color', color1, 'LineWidth',2, 'DisplayName','PF with Resample');
errorbar(1:k,xfp_re,reshape(Pfp_re,[1,length(Pfp_re)]),'Color',color1, 'LineWidth',1.5 ,'HandleVisibility', 'off');
legend('show');
% print('2a.eps','-depsc') ;

% the approximations to the posterior densities
sigma = 1;
k1 = 1;
k2 = 14;
k3 = 28;

% Particle filters without resampling 
figure()
plotPostPdf(k1, Xp(:,:,k1), Wp(:,k1)', xf, Pf, ~bResample, sigma, []); % k = 1

figure()
plotPostPdf(k2, Xp(:,:,k2), Wp(:,k2)', xf, Pf, ~bResample, sigma, []); % k = 14

figure()
plotPostPdf(k3, Xp(:,:,k3), Wp(:,k3)', xf, Pf, ~bResample, sigma, []); % k = 28

% Particle filters with resampling
figure()
plotPostPdf(k1, Xp_re(:,:,k1), Wp_re(:,k1)', xf, Pf, bResample, sigma, []); % k = 1

figure()
plotPostPdf(k2, Xp_re(:,:,k2), Wp_re(:,k2)', xf, Pf, bResample, sigma, []); % k = 14

figure()
plotPostPdf(k3, Xp_re(:,:,k3), Wp_re(:,k3)', xf, Pf, bResample, sigma, []); % k = 28

%% b)
% The initial prior
x_0_b = -20; 
P_0_b = 2;

% Kalman filter
[xf, Pf] = kalmanFilter(Y, x_0_b, P_0_b, A, Q, H, R);

% Particle filters without resampling
bResample = 1;
% plot = @(k, Xk, Xkmin1, j) plotPartTrajs(k, Xk, Xkmin1,[], j);
[xfp, Pfp, Xp, Wp] = pfFilter(x_0_b, P_0_b, Y, f, Q, h, R, N, ~bResample, [], 0);

% Particle filters with resampling
[xfp_re, Pfp_re, Xp_re, Wp_re] = pfFilter(x_0_b, P_0_b, Y, f, Q, h, R, N, bResample, [], 0);

% Plot 
figure()
grid on
hold on
plot(0:k, X, 'color', color5, 'LineWidth',2, 'DisplayName','True Trajectory');
plot(1:k, Y, '--', 'color', color4, 'LineWidth',2, 'DisplayName','Measurement sequence');
plot(1:k, xf, 'color', color3, 'LineWidth',2, 'DisplayName','KF estimate');
errorbar(1:k,xf,reshape(Pf,[1,length(Pf)]),'Color',color3, 'LineWidth',1.5, 'HandleVisibility', 'off');
plot(1:k, xfp, 'color', color2, 'LineWidth',2, 'DisplayName','PF without Resample');
errorbar(1:k,xfp,reshape(Pfp,[1,length(Pfp)]),'Color',color2, 'LineWidth',1.5 ,'HandleVisibility', 'off');
plot(1:k, xfp_re, 'color', color1, 'LineWidth',2, 'DisplayName','PF with Resample');
errorbar(1:k,xfp_re,reshape(Pfp_re,[1,length(Pfp_re)]),'Color',color1, 'LineWidth',1.5 ,'HandleVisibility', 'off');
legend('show');
print('2b.eps','-depsc') ;

%% c) PF without resampling
close all
color = [0.3010 0.7450 0.9330];
plotfunc = @(i, Xk, Xkmin1, Wk, j)plotPartTrajs(i, Xk, Xkmin1, Wk, j);

figure()
hold on
[xfp_c, Pfp_c, ~, ~] = pfFilter(x_0, P_0, Y, f, Q, h, R, N, ~bResample, plotfunc, 1);
plot(0:k,X, 'Color', color1, 'LineWidth', 2, 'DisplayName','True Trajectory');
plot(1:k,xfp_c, 'Color', color3, 'LineWidth', 2, 'DisplayName','Estimation');
errorbar(1:k,xfp_c, reshape(Pfp_c,[1,length(Pfp_c)]),'Color',color3, 'LineWidth',1.5 ,'HandleVisibility', 'off');
plot(NaN,NaN,'color',color,'DisplayName','Particle trajectories');
xlabel('Timesteps')
legend('show');
% print('2c.eps','-depsc') ;

%% d) PF with resampling
close all

figure()
hold on
[xfp_c, Pfp_c, ~, ~] = pfFilter(x_0, P_0, Y, f, Q, h, R, N, bResample, plotfunc, 1);
plot(0:k,X, 'Color', color1, 'LineWidth', 2, 'DisplayName','True Trajectory');
plot(1:k,xfp_c, 'Color', color3, 'LineWidth', 2, 'DisplayName','Estimation');
errorbar(1:k,xfp_c, reshape(Pfp_c,[1,length(Pfp_c)]),'Color',color3, 'LineWidth',1.5 ,'HandleVisibility', 'off');
plot(NaN,NaN,'color',color,'DisplayName','Particle trajectories');
xlabel('Timesteps')
legend('show');
% print('2d.eps','-depsc') ;



%% Problem 3
clear
clc

color1 = [0.4660 0.6740 0.1880];
color2 = [0.9290 0.6940 0.1250];
color3 = [0.8500 0.3250 0.0980];
color4 = [0 0.4470 0.7410];
color5 = [0.4940 0.1840 0.5560];
color6 = [0.6350 0.0780 0.1840];

%% a)
% Down in the MapProblemGetPoint.m function

%% b)
load('Xk.mat');
X = Xk(1,:);
Y = Xk(2,:);

% assume the x_0 by our self
X_0 = mvnrnd([11.5,9.5],0.7*eye(2))';

% X_k-1
X_previous = [X_0 Xk(:,1:end-1)];

% X_k - X_k-1
V_k = Xk - X_previous;

% Get the measurement
H = eye(2);
sigma_r = pi/180;
R = diag([sigma_r^2 sigma_r^2]);
V_measurement = genLinearMeasurementSequence(V_k, H, R);

% Plot
figure()
grid on
hold on
plot(V_measurement(1,:));
plot(V_measurement(2,:));
legend('Velocity of X','Velocity of Y')

%% c)
% Motivated in the report

%% d)
% Using CV model
T = 1;
N = 100;
bResample = 1;
Q = diag([0 0 0.1 0.1]);
sigma_r = 0.05;
R = diag([sigma_r^2 sigma_r^2]);

% The know initial prior
x_0 = [X(1) Y(1) V_k(1,2) V_k(2,2)]';
P_0 = diag([0 0 0 0]);

% functions
A = [1 0 T 0;
     0 1 0 T;
     0 0 1 0;
     0 0 0 1];

H = [0 0 1 0;
     0 0 0 1];

f = @(x) A*x;
h = @(x) H*x;


[xfp_c, ~, ~, ~] = pfFilter(x_0, P_0, V_measurement, f, Q, h, R, N, bResample, [], 0);

% Plot 
% position
figure()
hold on
plotfunctrack(X, Y, xfp_c)

% Velocity of X
figure()
grid on
hold on
plot(V_k(1,:),'Color',color1,'DisplayName','Velocity of X');
plot(xfp_c(3,:),'Color',color3,'DisplayName','Estimated velocity of X');
legend('show')

% Velocity of Y
figure()
grid on
hold on
plot(V_k(2,:),'Color',color1,'DisplayName','Velocity of Y');
plot(xfp_c(4,:),'Color',color3,'DisplayName','Estimated velocity of Y');
legend('show')

%% e)
% generate the random points on the road
% number of random points
n = 25000;

% Generate random points in the map
x_r = 0.8 + (11.2 - 0.8) * rand(n, 1);
y_r = 0.8 + (9.2 - 0.8) * rand(n, 1);

% check which points are on the road
u = isOnRoad(x_r, y_r);

% get the points that are on the road

points_on = [0 0]';
for i = 1:n
    if u(i) == 1
        points_on = [points_on,[x_r(i) y_r(i)]'];
    end
end

figure()
hold on
scatter(points_on(1,:), points_on(2,:),5,'filled');
axis([0.8 11.2 0.8 9.2])
title('Random Points on the Road');

% extra one point from the random points on the road to be the initial
% position
index = randi(length(points_on));
x_0e = [points_on(1,index), points_on(2,index),0,0]';
P_0e = diag([0 0 10 10]);

% PF filter
[xfp_e, ~, ~, ~] = pfFilter(x_0e, P_0e, V_measurement, f, Q, h, R, N, bResample, [], 0);

% Plot 
% position
figure()
hold on
plotfunctrack(X, Y, xfp_e)

% Velocity of X
figure()
grid on
hold on
plot(V_k(1,:),'Color',color1,'DisplayName','Velocity of X');
plot(xfp_e(3,:),'Color',color3,'DisplayName','Estimated velocity of X');
legend('show')

% Velocity of Y
figure()
grid on
hold on
plot(V_k(2,:),'Color',color1,'DisplayName','Velocity of Y');
plot(xfp_e(4,:),'Color',color3,'DisplayName','Estimated velocity of Y');
legend('show')

%% How to export your source code as .txt file. 
filename = fullfile('HA4.m');  % You should change "mainExample.m" with the name of your source file!
copyfile(filename,'HA4.txt','f') % Here, 'mainExample.txt' is the output. You should upload the 'main.txt' (or whatever you name it).

%% functions
function plotfunctrack(X, Y, xfp)
plot([1+i 1+9*i 5+9*i],'HandleVisibility', 'off')
plot([7+9*i 11+9*i 11+i 7+i],'HandleVisibility', 'off');
plot([5+i 1+i],'HandleVisibility', 'off')
plot([2+5.2*i 2+8.3*i 4+8.3*i 4+5.2*i 2+5.2*i],'HandleVisibility', 'off')%House 1
plot([2+3.7*i 2+4.4*i 4+4.4*i 4+3.7*i 2+3.7*i],'HandleVisibility', 'off')%House 2
plot([2+2*i 2+3.2*i 4+3.2*i 4+2*i 2+2*i],'HandleVisibility', 'off')%House 3
plot([5+i 5+2.2*i 7+2.2*i 7+i],'HandleVisibility', 'off')%House 4
plot([5+2.8*i 5+5.5*i 7+5.5*i 7+2.8*i 5+2.8*i],'HandleVisibility', 'off')%House 5
plot([5+6.2*i 5+9*i],'HandleVisibility', 'off');
plot([7+9*i 7+6.2*i 5+6.2*i],'HandleVisibility', 'off')%House 6
plot([8+4.6*i 8+8.4*i 10+8.4*i 10+4.6*i 8+4.6*i],'HandleVisibility', 'off')%House 7
plot([8+2.4*i 8+4*i 10+4*i 10+2.4*i 8+2.4*i],'HandleVisibility', 'off')%House 8
plot([8+1.7*i 8+1.8*i 10+1.8*i 10+1.7*i 8+1.7*i],'HandleVisibility', 'off')%House 9

% axis([0.8 11.2 0.8 9.2])
title('A map of the village','FontSize',20)
color1 = [0.4660 0.6740 0.1880];
color2 = [0.8500 0.3250 0.0980];

plot([X+Y*1i],'-*','Color',color1,'LineWidth', 1,'DisplayName','True Trajectory');
plot([xfp(1,:)+xfp(2,:)*1i],'-*','Color',color2,'LineWidth', 1, 'DisplayName','PF Estimation');
legend('show')
end