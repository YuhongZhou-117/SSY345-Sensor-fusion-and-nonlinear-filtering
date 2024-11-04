% clear
clc
%% Task 5
% [xhat, meas] = filterTemplateTask5();

%% Task 7
% [xhat, meas] = filterTemplateTask7();

%% Task 8
% [xhat, meas] = filterTemplateTask8();

%% Task 10
% [xhat, meas] = filterTemplateTask10();

%% Task 11
% [xhat, meas] = filterTemplateTask11();

%% Task 12
% [xhat, meas] = filterTemplateTask12();


% Removal of NaN values
% index =  ~any(isnan(meas.orient));

% t_meas = meas.t(index);
% orient_meas = meas.orient(:, ~any(isnan(meas.orient), 1));
% Euler_meas = q2euler(orient_meas);
% 
% t_xhat = xhat.t;
% 
% 
% Euler_xhat = q2euler(xhat.x);

[t_own, angler_own] = Eulerangler(xhat.t, xhat.x);
[t_phone, angler_phone] = Eulerangler(meas.t, meas.orient);


% Plot the Euler angles of both own orientation filter and the built in filter in the phone.
color1 = [0 0.4470 0.7410];
color2 = [0.8500 0.3250 0.0980];


figure

subplot(3,1,1)
hold on
plot(t_own, angler_own(3,:), 'Color', color1,'LineWidth',2, 'DisplayName', 'own filter')
plot(t_phone,angler_phone(3,:), 'Color', color2,'LineStyle','--','LineWidth', 2, 'DisplayName', 'phone filter')
legend('show')
title('X-axis')

subplot(3,1,2)
hold on
plot(t_own, angler_own(2,:), 'Color', color1,'LineWidth',2, 'DisplayName', 'own filter')
plot(t_phone,angler_phone(2,:), 'Color', color2,'LineStyle','--','LineWidth', 2, 'DisplayName', 'phone filter')
legend('show')
title('Y-axis')

subplot(3,1,3)
hold on
plot(t_own, angler_own(1,:), 'Color', color1,'LineWidth',2, 'DisplayName', 'own filter')
plot(t_phone,angler_phone(1,:), 'Color', color2,'LineStyle','--','LineWidth', 2, 'DisplayName', 'phone filter')
legend('show')
title('Z-axis')



function [te, angler] = Eulerangler(t, data)
    index =  ~any(isnan(data));
    te = t(index);
    datae = data(:, ~any(isnan(data), 1));
    angler = rad2deg(q2euler(datae));
end