close all
clear
clc

%% Problem 1
%% a)
N = 10000;
s1 = [0; 100];
s2 = [100; 0];
Sigma_phi = diag([(0.1*pi/180)^2 (0.1*pi/180)^2]);
h = @(x) dualBearingMeasurement(x, s1, s2);

% three different state densities
mean1 = [125; 125];
mean2 = [-25; 125];
mean3 = [60; 60];

Var = [10^2   0;
          0 5^2];

% Get N sample states
x1 = mvnrnd(mean1,Var,N)';
x2 = mvnrnd(mean2,Var,N)';
x3 = mvnrnd(mean3,Var,N)';

% generate samples of y = h(x) + r
y1 = genNonLinearMeasurementSequence(x1, h, Sigma_phi);
y2 = genNonLinearMeasurementSequence(x2, h, Sigma_phi);
y3 = genNonLinearMeasurementSequence(x3, h, Sigma_phi);

% approximate the mean and covariance
% y1
mu_y1 = mean(y1,2);
Sigma_y1 = (1/(N-1))*(y1 - mu_y1)*(y1 - mu_y1)';

fprintf('The approximated mean of y1 is\n%.4f\n%.4f\n',mu_y1)
fprintf('The approximated covariance of y1 is\n%.4f  %.4f\n%.4f  %.4f\n',Sigma_y1)

% y2
mu_y2 = mean(y2,2);
Sigma_y2 = (1/(N-1))*(y2 - mu_y2)*(y2 - mu_y2)';

fprintf('The approximated mean of y2 is\n%.4f\n%.4f\n',mu_y2)
fprintf('The approximated covariance of y2 is\n%.4f  %.4f\n%.4f  %.4f\n',Sigma_y2)

% y3
mu_y3 = mean(y3,2);
Sigma_y3 = (1/(N-1))*(y3 - mu_y3)*(y3 - mu_y3)';

fprintf('The approximated mean of y3 is\n%.4f\n%.4f\n',mu_y3)
fprintf('The approximated covariance of y3 is\n%.4f  %.4f\n%.4f  %.4f\n',Sigma_y3)

%% b)
% EKF
[mean1_EKF, cov1_EKF] = nonLinKFprediction(mean1, Var, h, Sigma_phi, 'EKF');
[mean2_EKF, cov2_EKF] = nonLinKFprediction(mean2, Var, h, Sigma_phi, 'EKF');
[mean3_EKF, cov3_EKF] = nonLinKFprediction(mean3, Var, h, Sigma_phi, 'EKF');

% UKF
[mean1_UKF, cov1_UKF] = nonLinKFprediction(mean1, Var, h, Sigma_phi, 'UKF');
[mean2_UKF, cov2_UKF] = nonLinKFprediction(mean2, Var, h, Sigma_phi, 'UKF');
[mean3_UKF, cov3_UKF] = nonLinKFprediction(mean3, Var, h, Sigma_phi, 'UKF');

% CKF
[mean1_CKF, cov1_CKF] = nonLinKFprediction(mean1, Var, h, Sigma_phi, 'CKF');
[mean2_CKF, cov2_CKF] = nonLinKFprediction(mean2, Var, h, Sigma_phi, 'CKF');
[mean3_CKF, cov3_CKF] = nonLinKFprediction(mean3, Var, h, Sigma_phi, 'CKF');


%% c)
% samples
xy1 = sigmaEllipse2D(mu_y1, Sigma_y1, 3, 100);
xy2 = sigmaEllipse2D(mu_y2, Sigma_y2, 3, 100);
xy3 = sigmaEllipse2D(mu_y3, Sigma_y3, 3, 100);

% EKF
xy1_EKF = sigmaEllipse2D(mean1_EKF, cov1_EKF, 3, 100);
xy2_EKF = sigmaEllipse2D(mean2_EKF, cov2_EKF, 3, 100);
xy3_EKF = sigmaEllipse2D(mean3_EKF, cov3_EKF, 3, 100);

% UKF
xy1_UKF = sigmaEllipse2D(mean1_UKF, cov1_UKF, 3, 100);
xy2_UKF = sigmaEllipse2D(mean2_UKF, cov2_UKF, 3, 100);
xy3_UKF = sigmaEllipse2D(mean3_UKF, cov3_UKF, 3, 100);

% CKF
xy1_CKF = sigmaEllipse2D(mean1_CKF, cov1_CKF, 3, 100);
xy2_CKF = sigmaEllipse2D(mean2_CKF, cov2_CKF, 3, 100);
xy3_CKF = sigmaEllipse2D(mean3_CKF, cov3_CKF, 3, 100);

% setting
alpha = 0.3;
color1 = [0.4660 0.6740 0.1880];
color2 = [0.9290 0.6940 0.1250];
color3 = [0.8500 0.3250 0.0980];
color4 = [0 0.4470 0.7410];

% Scenario 1
figure('Position',[300 300 800 500])
grid on
hold on
scatter(y1(1,:),y1(2,:),4,'filled', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha);
scatter(mu_y1(1,:),mu_y1(2,:),150,'o','MarkerEdgeColor',color1,'LineWidth',1.5);
plot(xy1(1,:),xy1(2,:),'-','color',color1,'Linewidth',2)

scatter(mean1_EKF(1,:),mean1_EKF(2,:),150,'+','MarkerEdgeColor',color2,'LineWidth',1.5);
plot(xy1_EKF(1,:),xy1_EKF(2,:),'--','color',color2,'Linewidth',2)

scatter(mean1_UKF(1,:),mean1_UKF(2,:),150,'d','MarkerEdgeColor',color3,'LineWidth',1.5);
plot(xy1_UKF(1,:),xy1_UKF(2,:),'-.','color',color3,'Linewidth',2)

scatter(mean1_CKF(1,:),mean1_CKF(2,:),150,'x','MarkerEdgeColor',color4,'LineWidth',1.5);
plot(xy1_CKF(1,:),xy1_CKF(2,:),'--','color',color4,'Linewidth',2)

legend('Samples','Sample mean','Sample covariance','EFK-mean', ...
    'EKF-covariance','UFK-mean', 'UKF-covariance','CKF-mean', 'CKF-covariance')
% print('1c1.eps','-depsc') ;

% Scenario 2
figure('Position',[300 300 800 500])
grid on
hold on
scatter(y2(1,:),y2(2,:),4,'filled', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha);
scatter(mu_y2(1,:),mu_y2(2,:),150,'o','MarkerEdgeColor',color1,'LineWidth',1.5);
plot(xy2(1,:),xy2(2,:),'-','color',color1,'Linewidth',2)

scatter(mean2_EKF(1,:),mean2_EKF(2,:),150,'+','MarkerEdgeColor',color2,'LineWidth',1.5);
plot(xy2_EKF(1,:),xy2_EKF(2,:),'--','color',color2,'Linewidth',2)

scatter(mean2_UKF(1,:),mean2_UKF(2,:),150,'d','MarkerEdgeColor',color3,'LineWidth',1.5);
plot(xy2_UKF(1,:),xy2_UKF(2,:),'-.','color',color3,'Linewidth',2)

scatter(mean2_CKF(1,:),mean2_CKF(2,:),150,'x','MarkerEdgeColor',color4,'LineWidth',1.5);
plot(xy2_CKF(1,:),xy2_CKF(2,:),'--','color',color4,'Linewidth',2)

legend('Samples','Sample mean','Sample covariance','EFK-mean', ...
    'EKF-covariance','UFK-mean', 'UKF-covariance','CKF-mean', 'CKF-covariance')
% print('1c2.eps','-depsc') ;

% Scenario 3
figure('Position',[300 300 800 500])
grid on
hold on
scatter(y3(1,:),y3(2,:),4,'filled', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha);
scatter(mu_y3(1,:),mu_y3(2,:),150,'o','MarkerEdgeColor',color1,'LineWidth',1.5);
plot(xy3(1,:),xy3(2,:),'-','color',color1,'Linewidth',2)

scatter(mean3_EKF(1,:),mean3_EKF(2,:),150,'+','MarkerEdgeColor',color2,'LineWidth',1.5);
plot(xy3_EKF(1,:),xy3_EKF(2,:),'--','color',color2,'Linewidth',2)

scatter(mean3_UKF(1,:),mean3_UKF(2,:),150,'d','MarkerEdgeColor',color3,'LineWidth',1.5);
plot(xy3_UKF(1,:),xy3_UKF(2,:),'-.','color',color3,'Linewidth',2)

scatter(mean3_CKF(1,:),mean3_CKF(2,:),150,'x','MarkerEdgeColor',color4,'LineWidth',1.5);
plot(xy3_CKF(1,:),xy3_CKF(2,:),'--','color',color4,'Linewidth',2)

legend('Samples','Sample mean','Sample covariance','EFK-mean', ...
    'EKF-covariance','UFK-mean', 'UKF-covariance','CKF-mean', 'CKF-covariance')
% print('1c3.eps','-depsc') ;



%% Problem 2
% clear
% close all
% clc
% color1 = [0.4660 0.6740 0.1880];
% color2 = [0.9290 0.6940 0.1250];
% color3 = [0.8500 0.3250 0.0980];
% color4 = [0 0.4470 0.7410];

% Time steps
N = 100;

% Sampling time
T = 1;

% Sensor location
s1 = [-200 100]';
s2 = [-200 -100]';

% The initial prior
x_0 = [0 0 20 0 5*pi/180]';
P_0 = diag([10^2 10^2 2^2 (pi/180)^2 (pi/180)^2]);

% Noise parameters
Sigma_v = 1;
Sigma_w = pi/180;

% Process noise
Q = diag([0 0 Sigma_v 0 Sigma_w]).^2;

% Function
f = @(x)coordinatedTurnMotion(x, T);
h = @(x) dualBearingMeasurement(x, s1, s2);

% Case 1 2 3
Case = 1;

switch Case
    case 1
        Sigma_phi1 = 2*pi/180;
        Sigma_phi2 = 2*pi/180;
    case 2
        Sigma_phi1 = 2*pi/180;
        Sigma_phi2 = 0.1*pi/180;
    case 3
        Sigma_phi1 = 0.1*pi/180;
        Sigma_phi2 = 0.1*pi/180;
    otherwise
        error('Incorrect Case')
end

%% a) b)
% Measurement noise
R = diag([Sigma_phi1^2, Sigma_phi2^2]);

% Generate the state sequence
X = genNonLinearStateSequence(x_0, P_0, f, Q, N);

% Generate the measurement sequence
Y = genNonLinearMeasurementSequence(X, h, R);

% Get the position from the measurements
pos = zeros(100,2); 
for i = 1:N
    [x_pos, y_pos] = getPosFromMeasurement(Y(1,i), Y(2,i), s1, s2);
    pos(i,:) = [x_pos, y_pos];
end

% The three different non-linear Kalman filters
% EKF
[xf_EKF, Pf_EKF, ~, ~] = nonLinearKalmanFilter(Y, x_0, P_0, f, Q, h, R, 'EKF');

% UKF
[xf_UKF, Pf_UKF, ~, ~] = nonLinearKalmanFilter(Y, x_0, P_0, f, Q, h, R, 'UKF');

% CKF
[xf_CKF, Pf_CKF, ~, ~] = nonLinearKalmanFilter(Y, x_0, P_0, f, Q, h, R, 'CKF');

% plot the result
% EKF
figure()
grid on 
hold on
% Sensor position
scatter(s1(1),s1(2),150,'o','LineWidth',2);
scatter(s2(1),s2(2),150,'+','LineWidth',2);
% rue positions
plot(X(1,:), X(2,:),'color',color1,'LineWidth', 1.5);
% estimated position
plot(xf_EKF(1,:), xf_EKF(2,:),'color',color2,'LineWidth', 1.5);
% measurements
plot(pos(1:end,1),pos(1:end,2),'--','color',color3,'LineWidth', 1.5);
% 3σ-contours at every 5th estimate
for i= 1:5:length(xf_EKF)
        xy = sigmaEllipse2D(xf_EKF(1:2,i), Pf_EKF(1:2,1:2,i), 3, 50);
        plot(xy(1,:),xy(2,:),'color',color4);
        % fill(xy(1,:),xy(2,:),color4,'FaceAlpha', 0.3);
end
legend('Sensor1','Sensor2','True Position','Estimated position','Measurements','3\sigma-contours at 5th')
subtitle('EKF')
% print('2aEKF.eps','-depsc') ;
% print('2a2EKF.eps','-depsc') ;
% print('2a3EKF.eps','-depsc') ;

% UKF
figure()
grid on 
hold on
% Sensor position
scatter(s1(1),s1(2),150,'o','LineWidth',2);
scatter(s2(1),s2(2),150,'+','LineWidth',2);
% rue positions
plot(X(1,:), X(2,:),'color',color1,'LineWidth', 1.5);
% estimated position
plot(xf_UKF(1,:), xf_UKF(2,:),'color',color2,'LineWidth', 1.5);
% measurements
plot(pos(1:end,1),pos(1:end,2),'--','color',color3,'LineWidth', 1.5);
% 3σ-contours at every 5th estimate
for i= 1:5:length(xf_UKF)
        xy = sigmaEllipse2D(xf_UKF(1:2,i), Pf_UKF(1:2,1:2,i), 3, 50);
        plot(xy(1,:),xy(2,:),'color',color4);
        % fill(xy(1,:),xy(2,:),color4,'FaceAlpha', 0.3);
end
legend('Sensor1','Sensor2','True Position','Estimated position','Measurements','3\sigma-contours at 5th')
subtitle('UKF')
% print('2aUKF.eps','-depsc') ;
% print('2a2UKF.eps','-depsc') ;
% print('2a3UKF.eps','-depsc') ;

% CKF
figure()
grid on 
hold on
% Sensor position
scatter(s1(1),s1(2),150,'o','LineWidth',2);
scatter(s2(1),s2(2),150,'+','LineWidth',2);
% rue positions
plot(X(1,:), X(2,:),'color',color1,'LineWidth', 1.5);
% estimated position
plot(xf_CKF(1,:), xf_CKF(2,:),'color',color2,'LineWidth', 1.5);
% measurements
plot(pos(1:end,1),pos(1:end,2),'--','color',color3,'LineWidth', 1.5);
% 3σ-contours at every 5th estimate
for i= 1:5:length(xf_UKF)
        xy = sigmaEllipse2D(xf_CKF(1:2,i), Pf_CKF(1:2,1:2,i), 3, 50);
        plot(xy(1,:),xy(2,:),'color',color4);
        % fill(xy(1,:),xy(2,:),color4,'FaceAlpha', 0.3);
end
legend('Sensor1','Sensor2','True Position','Estimated position','Measurements','3\sigma-contours at 5th')
subtitle('CKF')
% print('2aCKF.eps','-depsc') ;
% print('2a2CKF.eps','-depsc') ;
% print('2a3CKF.eps','-depsc') ;

%% c)
MC = 150;
error_x_EKF = cell(1, MC);
error_y_EKF = cell(1, MC);

error_x_UKF = cell(1, MC);
error_y_UKF = cell(1, MC);

error_x_CKF = cell(1, MC);
error_y_CKF = cell(1, MC);

for imc = 1:MC
    % Simulate state sequence
    X = genNonLinearStateSequence(x_0, P_0, f, Q, N);

    % Simulate measurements
    Y = genNonLinearMeasurementSequence(X, h, R);

    % Run Kalman filter (you need to run all three, for comparison)
    % EKF
    [xf_EKF,Pf_EKF,~,~] = nonLinearKalmanFilter(Y,x_0,P_0,f,Q,h,R,'EKF');
    
    % UKF
    [xf_UKF,Pf_UKF,~,~] = nonLinearKalmanFilter(Y,x_0,P_0,f,Q,h,R,'UKF');

    % CKF
    [xf_CKF,Pf_CKF,~,~] = nonLinearKalmanFilter(Y,x_0,P_0,f,Q,h,R,'CKF');

    % Save the estimation errors and the prediction errors!
    error_x_EKF{imc} = X(1,2:end) - xf_EKF(1,:);
    error_y_EKF{imc} = X(2,2:end) - xf_EKF(2,:);

    error_x_UKF{imc} = X(1,2:end) - xf_UKF(1,:);
    error_y_UKF{imc} = X(2,2:end) - xf_UKF(2,:);
    
    error_x_CKF{imc} = X(1,2:end) - xf_CKF(1,:);
    error_y_CKF{imc} = X(2,2:end) - xf_CKF(2,:);
end

error_x_EKF = cell2mat(error_x_EKF);
error_y_EKF = cell2mat(error_y_EKF);

error_x_UKF = cell2mat(error_x_UKF);
error_y_UKF = cell2mat(error_y_UKF);

error_x_CKF = cell2mat(error_x_CKF);
error_y_CKF = cell2mat(error_y_CKF);

% plot the result
figure('Position',[300 300 1200 800])
grid on
subplot(2, 3, 1)
hold on
plotc(error_x_EKF,N)
title('X-position-EKF')
subplot(2, 3, 2)
hold on
plotc(error_x_UKF,N)
title('X-position-UKF')
subplot(2, 3, 3)
hold on
plotc(error_x_CKF,N)
title('X-position-CKF')
subplot(2, 3, 4)
hold on
plotc(error_y_EKF,N)
title('Y-position-EKF')
subplot(2, 3, 5)
hold on
plotc(error_y_UKF,N)
title('Y-position-UKF')
subplot(2, 3, 6)
hold on
plotc(error_y_CKF,N)
title('Y-position-CKF')
% print('2c1.eps','-depsc') ;
% print('2c2.eps','-depsc') ;
% print('2c3.eps','-depsc') ;

%% Problem 3
% close all
%% True track
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

% sensor position 
s1 = [300 -100]';
s2 = [300 -300]';

% The measurement noise standard deviations
Sigma_phi1 = pi/180;
Sigma_phi2 = pi/180;
R = diag([Sigma_phi1^2, Sigma_phi2^2]);

% function
f = @(x) coordinatedTurnMotion(x, T);
h = @(x) dualBearingMeasurement(x, s1, s2);

%% a) b) c)
Sigma_v = 1*0.1;
Sigma_w = pi/180*0.16;
Q = diag([0 0 Sigma_v 0 Sigma_w]).^2;
Y = genNonLinearMeasurementSequence(X, h, R);

% Get the position from the measurements
pos3 = zeros(600,2); 
for i = 1:K
    [x_pos, y_pos] = getPosFromMeasurement(Y(1,i), Y(2,i), s1, s2);
    pos3(i,:) = [x_pos, y_pos];
end

[xf,Pf,xp,Pp] = nonLinearKalmanFilter(Y,x_0,P_0,f,Q,h,R,'CKF');

% plot the result
figure()
grid on
hold on
% Sensor position
scatter(s1(1),s1(2),150,'o','LineWidth',2);
scatter(s2(1),s2(2),150,'+','LineWidth',2);
% True positions
plot(X(1,:), X(2,:),'color',color1,'LineWidth', 1.5);
% estimated positions
plot(xf(1,:), xf(2,:),'color',color2, 'LineWidth', 1.5);
% measurement positions
plot(pos3(1:end,1),pos3(1:end,2),'--','color',color3,'LineWidth', 1.5);
% 3σ-contours at every 5th estimate
for i= 1:5:length(xf)
        xy = sigmaEllipse2D(xf(1:2,i), Pf(1:2,1:2,i), 3, 50);
        plot(xy(1,:),xy(2,:),'color',color4);
end
legend('Sensor1','Sensor2','True Position','Estimated positions', ...
    'Measurements','$3\sigma-contours at 5th$','Interpreter', 'latex')


% error
error  = vecnorm(X(1:2,2:end)-xf(1:2,:));
figure('Position',[300 300 600 200])
plot(0:599,error)
legend('$||P_k - \hat{P}_{k|k}||_2$','Interpreter', 'latex')

%% How to export your source code as .txt file. 
% filename = fullfile('HA3.m');  % You should change "mainExample.m" with the name of your source file!
% copyfile(filename,'HA3.txt','f') % Here, 'mainExample.txt' is the output. You should upload the 'main.txt' (or whatever you name it).

%% FUNCTIONS
function plotc(error,N)
    histogram(error, 'Normalization', 'pdf')
    mu = mean(error);
    sigma = std(error)^2;
    x = linspace(mu-3*sqrt(sigma), mu+3*sqrt(sigma), N);
    y = normpdf(x, mu, sqrt(sigma));
    plot(x,y,'LineWidth',1.5)
    lgd = legend('Histograms','Gaussian distribution');
    title(lgd,['$\mu$:',num2str(mu),'$\sigma$:',num2str(sigma)],'Interpreter', 'latex')
end