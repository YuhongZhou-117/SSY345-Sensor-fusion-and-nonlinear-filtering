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

% Pf filter
[xfp_d, ~, ~, ~] = pfFilter(x_0, P_0, V_measurement, f, Q, h, R, N, bResample, [], 0);
plotfunc = @(i, Xk, Xkmin1, Wk, j)plotPartTrajs(i, Xk, Xkmin1, Wk, j);

% Plot 
% position
figure()
hold on
% [xfp_d, ~, ~, ~] = pfFilter(x_0, P_0, V_measurement, f, Q, h, R, N, bResample, plotfunc, 1)
plotfunctrack(X, Y, xfp_d)

% Velocity of X
figure()
grid on
hold on
% [xfp_d, ~, ~, ~] = pfFilter(x_0, P_0, V_measurement, f, Q, h, R, N, bResample, plotfunc, 1);
plot(V_k(1,:),'Color',color1,'DisplayName','Velocity of X');
plot(xfp_d(3,:),'Color',color3,'DisplayName','Estimated velocity of X');
legend('show')

% Velocity of Y
figure()
grid on
hold on
plot(V_k(2,:),'Color',color1,'DisplayName','Velocity of Y');
plot(xfp_d(4,:),'Color',color3,'DisplayName','Estimated velocity of Y');
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