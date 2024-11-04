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
N = 1000;

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
plotPostPdf(k1, Xp(:,:,k1), Wp(:,k1)', xf, Pf, ~bResample, 1, []); % k = 1
% print('2a1.eps','-depsc')

figure()
plotPostPdf(k2, Xp(:,:,k2), Wp(:,k2)', xf, Pf, ~bResample, 1, []); % k = 14
% print('2a2.eps','-depsc')

figure()
plotPostPdf(k3, Xp(:,:,k3), Wp(:,k3)', xf, Pf, ~bResample, 1, []); % k = 28
% print('2a3.eps','-depsc')

% Particle filters with resampling
figure()
plotPostPdf(k1, Xp_re(:,:,k1), Wp_re(:,k1)', xf, Pf, bResample, 1.5, []); % k = 1
% print('2a4.eps','-depsc')

figure()
plotPostPdf(k2, Xp_re(:,:,k2), Wp_re(:,k2)', xf, Pf, bResample, 1.8, []); % k = 14
% print('2a5.eps','-depsc')

figure()
plotPostPdf(k3, Xp_re(:,:,k3), Wp_re(:,k3)', xf, Pf, bResample, 1.5, []); % k = 28
% print('2a6.eps','-depsc')
