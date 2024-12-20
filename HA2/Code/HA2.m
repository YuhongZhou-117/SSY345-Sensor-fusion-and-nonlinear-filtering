close all
clear
clc
%% Scenario 1 – A first Kalman filter and its properties
%% a)
Q = 1.5;
R = 3;

% model matrix
A = 1;
H = 1;

% initial prior
x_0 = 2;
P_0 = 8;

% length of the sequence
N = 35;

% Generate a state sequence
X_1a = genLinearStateSequence(x_0, P_0, A, Q, N);

% Generate a measurement sequence
Y_1a = genLinearMeasurementSequence(X_1a, H, R);

% Plot the results
figure()
hold on
grid on
plot(0:N, X_1a,'Linewidth', 2)
plot(1:N, Y_1a,'linewidth', 2)
legend('State sequence X','Measurement sequence Y')
xlabel('Sequence length N');
% print('1a.eps','-depsc') ;

%% 1b)
% kalman filter
[X_est, P_est] = kalmanFilter(Y_1a, x_0, P_0, A, Q, H, R);

% x estimate + - 2 sigma
X_sigma1 = X_est + 3*sqrt(P_est);
X_sigma2 = X_est - 3*sqrt(P_est);

% Plot the results
figure()
hold on
grid on
plot(0:N, X_1a,'Linewidth', 1.5)
plot(1:N, Y_1a,'Linewidth', 1.5)
plot(1:N, X_est,'Linewidth', 1.5)
plot(1:N, X_sigma1(1,1:N),'--','color',[0.4660 0.6740 0.1880],'Linewidth', 1.5)
plot(1:N, X_sigma2(1,1:N),'--','color',[0.4660 0.6740 0.1880],'Linewidth', 1.5)
legend('Correct states','Measurements','$\hat{X_k}$','$\hat{X_k}+3\sigma$','$\hat{X_k}-3\sigma$','Interpreter', 'latex')
xlabel('Sequence length N');
% print('1b1.eps','-depsc') ;

% plot the error density
k = [1; 2; 4; 30];
range = [-10:0.1:10];

figure()
hold on
grid on
for i = 1:length(k)
    error_pdf = normpdf(range, 0, sqrt(P_est(i)));
    plot(range,error_pdf,'LineWidth',1.5)
end
legend('K=1', 'K=2','K=4','K=30')
xlabel('error')
ylabel('pdf')
% print('1b2.eps','-depsc') ;

%% 1c)
% kalman filter with incorrect prior x_0 = 12
[X_est_in, ~] = kalmanFilter(Y_1a, 12, P_0, A, Q, H, R);

% Plot the results
figure()
hold on
grid on
plot(0:N, X_1a,'Linewidth', 1.5)
plot(1:N, X_est,'--','color',[0.4660 0.6740 0.1880],'linewidth', 1.5)
plot(1:N, X_est_in,'--r','linewidth', 1.5)
legend('Correct states','Correct Kalman filter','Incorrect Kalman filter')
xlabel('Sequence length N');
% print('1c.eps','-depsc') ;


%% 1d)
%  choice of k
k_1d = 20;
x_range = [X_est(:,k_1d)-10:0.1:X_est(:,k_1d)+10];


% x_k-1 given y1_k-1
d1 = normpdf(x_range,X_est(k_1d-1),sqrt(P_est(k_1d-1)));

% x_k given y1_k-1
[x_pre,P_pre]=linearPrediction(X_est(:,k_1d-1),P_est(:,:,k_1d-1),A,Q); 
d2 = normpdf(x_range,x_pre,sqrt(P_pre));

% x_k given y1_k 
d3 = normpdf(x_range,X_est(k_1d),sqrt(P_est(k_1d)));

% plot the result
figure()
hold on
grid on
plot(x_range,d1,'LineWidth',1.5)
plot(x_range,d2,'LineWidth',1.5)
plot(x_range,d3,'LineWidth',1.5)
xline(Y_1a(:,k_1d),'LineWidth',1.5);
legend('$p(x_{k-1}|y_{1:k-1})$','$p(x_{k}|y_{1:k-1})$','$p(x_{k}|y_{1:k})$','$y_k$','Interpreter', 'latex')
% print('1d.eps','-depsc') ;

%% 1e)
Ne = 10000;

% Generate a state sequence
X_1e = genLinearStateSequence(x_0, P_0, A, Q, Ne);

% Generate a measurement sequence
Y_1e = genLinearMeasurementSequence(X_1e, H, R);

% Filter the measurement sequence using the Kalman filter
[X_est_1e, P_est_1e] = kalmanFilter(Y_1e, x_0, P_0, A, Q, H, R);

% estimation error
error = X_1e(2:end) - X_est_1e;

% the pdf of N(x;0,P_N|N)
range_1e = [-5*sqrt(P_est_1e(:,:,end)):0.1:5*sqrt(P_est_1e(:,:,end))];
pdf_1e = normpdf(range_1e, 0, sqrt(P_est_1e(:,:,end)));

% plot the result
figure()
hold on
grid on
histogram(error,range_1e,'Normalization','pdf')
plot(range_1e,pdf_1e,'LineWidth',2)
legend('$x_k - \hat{x}_{k|k}$','$\mathcal{N}(x;0,P_{N|N})$','Interpreter', 'latex')
% print('1e1.eps','-depsc') ;

% innovation v_k = y_k - H*x_k|k-1
[~, P_v,vk,Sk] = kalmanFilterVk(Y_1e, x_0, P_0, A, Q, H, R);

% mean
mean_vk = mean(vk);

% plot the result
figure()
hold on
grid on
plot(Sk^(-1/2)*vk)
yline(mean_vk,'color',[1,77/255,0],'LineWidth',2)
yline(3*sqrt(P_v(length(Y_1e))),'--r')
yline(-3*sqrt(P_v(length(Y_1e))),'--r')
legend('$S_k^{-1/2}v_k$',['mean:',num2str(mean_vk)],'$3\sigma$','$-3\sigma$','Interpreter', 'latex')
xlabel('k')
ylabel('$S_k^{-1/2}v_k$','Interpreter', 'latex')
xlim([0,10000])
% print('1e2.eps','-depsc') ;


figure()
hold on
grid on
autocorr(vk)
xlabel('Lag')
ylabel('Sample Autocorrelation')
% print('1e3.eps','-depsc') ;


%% Scenario 2 - Tuning a Kalman filter
%% 2a)
load SensorMeasurements.mat;

% y_k = C(v_k + r_k)
% stationary 
v0 = CalibrationSequenceVelocity_v0;

% velocity of the train was constantly 10 m/s
v10 = CalibrationSequenceVelocity_v10;

% velocity of the train was constantly 20 m/s 
v20 = CalibrationSequenceVelocity_v20;

% estimate c
C10 = (mean(v10)-mean(v0))/10;
C20 = (mean(v20)-mean(v0))/20;
C = mean([C10,C20]);

fprintf('The scaling constant C is %f\n', C)

% variance var[y_k] = C^2 var[r_k]
var0 = var(v0)/(C^2);
var10 = var(v10)/(C^2);
var20 = var(v20)/(C^2);
var_rk = mean([var0,var10,var20]);

fprintf('The variance of the measurement noise is %f', var_rk)


%% 2b)
Y_2b = Generate_y_seq();

n_2b = size(Y_2b,2);
temple = [];
for i = 1:n_2b
    if ~isnan(Y_2b(1,i))
        temple = [temple,Y_2b(:,i)];
    end
end
Y_2b = temple;


% function [X, P,vk] = kalmanFilter_update(Y, x_0, P_0, A, Q, H, R)
% %KALMANFILTER Filters measurements sequence Y using a Kalman filter. 
% %
% %Input:
% %   Y           [m x N] Measurement sequence
% %   x_0         [n x 1] Prior mean
% %   P_0         [n x n] Prior covariance
% %   A           [n x n] State transition matrix
% %   Q           [n x n] Process noise covariance
% %   H           [m x n] Measurement model matrix
% %   R           [m x m] Measurement noise covariance
% %
% %Output:
% %   x           [n x N] Estimated state vector sequence
% %   P           [n x n x N] Filter error convariance
% %
% 
% %% Parameters
% N = size(Y,2);
% 
% n = length(x_0);
% m = size(Y,1);
% 
% %% remove the data with NaN to form the new Y
% temple = [];
% for k = 1:N
%     if isnan(Y(:,k)) ~= 1
%         temple = [temple,Y(:,k)];
%     end
% end
% 
% Y = temple;
% N = size(Y,2);
% 
% %% Data allocation
% x = zeros(n,N);
% P = zeros(n,n,N);
% vk= zeros(m,N);
% %% Initial 
% x(:,1) = x_0;
% P(:,:,1) = P_0;
% 
% for i =2:N+1
%     % prediction
%     [x_pred, P_pred] = linearPrediction(x(:,i-1), P(:,:,i-1), A, Q);
% 
%     % update
%     [x(:,i), P(:,:,i),v_k] = linearUpdateVk(x_pred, P_pred, Y(:,i-1), H, R);    
% 
%     % innovation
%     vk(:,i)=v_k;
% end
% 
% X = x(:,2:end);
% P = P(:,:,2:end);
% 
% end



%% 2c)
% Constant Velocity Model
Y_2c = Generate_y_seq();
T = 0.2;
R_2c = 3;

A_cv = [1, T;
        0, 1];

H_cv = [1, 0;
        0, C];

R_cv = [1,        0;
        0, C^2*R_2c];

% Constant Acceleratiion Model
A_ca = [1, T, T^2/2;
        0, 1,     T;
        0, 0,     1];

H_ca = [1, 0, 0;
        0, C, 0];

R_ca = [1,        0;
        0, C^2*R_2c];

% initially stationary
% CV Model
x1_0 = [0; 0];
P1_0 = eye(2);

% CA Model
x2_0 = [0;0;0];
P2_0 = eye(3);

% test1
Q_1 = 0.0001;
Q_cv = [0,   0;
        0, Q_1];

Q_ca = [0, 0,   0;
        0, 0,   0;
        0, 0, Q_1];

% CV Model
[X_cv1, ~,~,Y_2cn] = kalmanFilter_update(Y_2c, x1_0, P1_0, A_cv, Q_cv, H_cv, R_cv);

% CA Model
[X_ca1, ~,~] = kalmanFilter_update(Y_2c, x2_0, P2_0, A_ca, Q_ca, H_ca, R_ca);

% plot the result
figure()
subplot(2,1,1)
hold on 
grid on
plot(1:length(X_cv1),X_cv1(1,:),'linewidth',1.5)
plot(1:length(X_cv1),X_ca1(1,:),'--','color',[0.4660 0.6740 0.1880],'linewidth',1.5)
plot(1:length(X_cv1),Y_2cn(1,:),'-.','linewidth',1)
legend('CV Model','CA Model','Measurment')
title('Position')
subplot(2,1,2)
hold on 
grid on
plot(1:length(X_cv1),X_cv1(2,:),'linewidth',1.5)
plot(1:length(X_cv1),X_ca1(2,:),'--','color',[0.4660 0.6740 0.1880],'linewidth',1.5)
plot(1:length(X_cv1),Y_2cn(2,:)./C,'linewidth',1)
legend('CV Model','CA Model','Measurment')
title('Velocity')
% print('2c1.eps','-depsc') ;

% test2
Q_2 = 1;
Q_cv = [0,   0;
        0, Q_2];

Q_ca = [0, 0,   0;
        0, 0,   0;
        0, 0, Q_2];

% CV Model
[X_cv1, ~,~,Y_2cn] = kalmanFilter_update(Y_2c, x1_0, P1_0, A_cv, Q_cv, H_cv, R_cv);

% CA Model
[X_ca1, ~,~,~] = kalmanFilter_update(Y_2c, x2_0, P2_0, A_ca, Q_ca, H_ca, R_ca);

% plot the result
figure()
subplot(2,1,1)
hold on 
grid on
plot(1:length(X_cv1),X_cv1(1,:),'linewidth',1.5)
plot(1:length(X_cv1),X_ca1(1,:),'--','color',[0.4660 0.6740 0.1880],'linewidth',1.5)
plot(1:length(X_cv1),Y_2cn(1,:),'-.','linewidth',1)
legend('CV Model','CA Model','Measurment')
title('Position')
subplot(2,1,2)
hold on 
grid on
plot(1:length(X_cv1),X_cv1(2,:),'linewidth',1.5)
plot(1:length(X_cv1),X_ca1(2,:),'--','color',[0.4660 0.6740 0.1880],'linewidth',1.5)
plot(1:length(X_cv1),Y_2cn(2,:)./C,'linewidth',1)
legend('CV Model','CA Model','Measurment')
title('Velocity')
% print('2c2.eps','-depsc') ;


% test3
Q_3 = 10000;
Q_cv = [0,   0;
        0, Q_3];

Q_ca = [0, 0,   0;
        0, 0,   0;
        0, 0, Q_3];

% CV Model
[X_cv1,~,~,Y_2cn] = kalmanFilter_update(Y_2c, x1_0, P1_0, A_cv, Q_cv, H_cv, R_cv);

% CA Model
[X_ca1,~,~] = kalmanFilter_update(Y_2c, x2_0, P2_0, A_ca, Q_ca, H_ca, R_ca);

% plot the result
figure()
subplot(2,1,1)
hold on 
grid on
plot(1:length(X_cv1),X_cv1(1,:),'linewidth',1.5)
plot(1:length(X_cv1),X_ca1(1,:),'--','color',[0.4660 0.6740 0.1880],'linewidth',1.5)
plot(1:length(X_cv1),Y_2cn(1,:),'-.','linewidth',1)
legend('CV Model','CA Model','Measurment')
title('Position')
subplot(2,1,2)
hold on 
grid on
plot(1:length(X_cv1),X_cv1(2,:),'linewidth',1.5)
plot(1:length(X_cv1),X_ca1(2,:),'--','color',[0.4660 0.6740 0.1880],'linewidth',1.5)
plot(1:length(X_cv1),Y_2cn(2,:)./C,'linewidth',1)
legend('CV Model','CA Model','Measurment')
title('Velocity')
% print('2c3.eps','-depsc') ;


%% How to export your source code as .txt file. 
% filename = fullfile('HA2.m');  % You should change "mainExample.m" with the name of your source file!
% copyfile(filename,'HA2.txt','f') % Here, 'mainExample.txt' is the output. You should upload the 'main.txt' (or whatever you name it).



%% function
function [x, P,vk] = linearUpdateVk(x, P, y, H, R)
%LINEARPREDICTION calculates mean and covariance of predicted state
%   density using a linear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   y           [m x 1] Measurement
%   H           [m x n] Measurement model matrix
%   R           [m x m] Measurement noise covariance
%
%Output:
%   x           [n x 1] updated state mean
%   P           [n x n] updated state covariance
%

% Your code here
% Kalman filter gain
K = P*H'*inv(H*P*H' + R);

% innovation
vk= y - H*x;

% Update
x = x + K*vk;
P = P - K*H*P;
end

function [X, P,vk,Sk] = kalmanFilterVk(Y, x_0, P_0, A, Q, H, R)
%KALMANFILTER Filters measurements sequence Y using a Kalman filter. 
%
%Input:
%   Y           [m x N] Measurement sequence
%   x_0         [n x 1] Prior mean
%   P_0         [n x n] Prior covariance
%   A           [n x n] State transition matrix
%   Q           [n x n] Process noise covariance
%   H           [m x n] Measurement model matrix
%   R           [m x m] Measurement noise covariance
%
%Output:
%   x           [n x N] Estimated state vector sequence
%   P           [n x n x N] Filter error convariance
%

%% Parameters
N = size(Y,2);

n = length(x_0);
m = size(Y,1);

%% Data allocation
x  = zeros(n,N);
P  = zeros(n,n,N);
vk = zeros(m,N);
%% Initial 
x(:,1) = x_0;
P(:,:,1) = P_0;

for i =2:N+1
    % prediction
    [x_pred, P_pred] = linearPrediction(x(:,i-1), P(:,:,i-1), A, Q);


    Sk = H*P_pred*H' + R;
    % update
    [x(:,i), P(:,:,i),v_k] = linearUpdateVk(x_pred, P_pred, Y(:,i-1), H, R);    

    % innovation
    vk(:,i)=v_k;
end

X = x(:,2:end);
P = P(:,:,2:end);

end

function [X, P,vk,Y] = kalmanFilter_update(Y, x_0, P_0, A, Q, H, R)
%KALMANFILTER Filters measurements sequence Y using a Kalman filter. 
%
%Input:
%   Y           [m x N] Measurement sequence
%   x_0         [n x 1] Prior mean
%   P_0         [n x n] Prior covariance
%   A           [n x n] State transition matrix
%   Q           [n x n] Process noise covariance
%   H           [m x n] Measurement model matrix
%   R           [m x m] Measurement noise covariance
%
%Output:
%   x           [n x N] Estimated state vector sequence
%   P           [n x n x N] Filter error convariance
%

%% Parameters
N = size(Y,2);

n = length(x_0);
m = size(Y,1);

%% remove the data with NaN to form the new Y
temple = [];
for k = 1:N
    if ~isnan(Y(1,k))
        temple = [temple,Y(:,k)];
    end
end

Y = temple;
N = size(Y,2);
m = size(Y,1);

%% Data allocation
x = zeros(n,N);
P = zeros(n,n,N);
vk= zeros(m,N);
%% Initial 
x(:,1) = x_0;
P(:,:,1) = P_0;

for i =2:N+1
    % prediction
    [x_pred, P_pred] = linearPrediction(x(:,i-1), P(:,:,i-1), A, Q);
    
    % update
    [x(:,i), P(:,:,i),v_k] = linearUpdateVk(x_pred, P_pred, Y(:,i-1), H, R);    

    % innovation
    vk(:,i)=v_k;
end

X = x(:,2:end);
P = P(:,:,2:end);

end

