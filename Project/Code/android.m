close all
% clear
clc

% color
color1 = '#0072BD';
color2 = '#D95319';
color3 = '#EDB120';
color4 = '#7E2F8E';
color5 = '#77AC30';

% [xhat, meas] = filterTemplate();

%% Task 2
t = meas.t;

% Plot measurement
plotmeasurement(meas.acc,t,'Accelerometer');

plotmeasurement(meas.gyr,t,'Gyroscope');

plotmeasurement(meas.mag,t,'Magnetometer');

% plotmeasurement(meas.orient,t,'Orientation');

% mean
measure_acc = meas.acc(:, ~any(isnan(meas.acc), 1));
measure_gyr = meas.gyr(:, ~any(isnan(meas.gyr), 1));
measure_mag = meas.mag(:, ~any(isnan(meas.mag), 1));

mean_acc = mean(measure_acc, 2);
printm(mean_acc,'mean of accelerometer')
mean_gyr = mean(measure_gyr, 2);
printm(mean_gyr,'mean of gyroscope')
mean_mag = mean(measure_mag, 2);
printm(mean_mag,'mean of magnetometer')

% covariance
% Accelerometer
cov_accX = cov(measure_acc(1,:));
fprintf('covariance of accelerometer in X: %e\n',cov_accX);
cov_accY = cov(measure_acc(2,:));
fprintf('covariance of accelerometer in Y: %e\n',cov_accY);
cov_accZ = cov(measure_acc(3,:));
fprintf('covariance of accelerometer in Z: %e\n',cov_accZ);

cov_acc = [cov_accX;cov_accY;cov_accZ];

% Gyroscope
cov_gyrX = cov(measure_gyr(1,:));
fprintf('covariance of Gyroscope in X: %e\n',cov_gyrX);
cov_gyrY = cov(measure_gyr(2,:));
fprintf('covariance of Gyroscope in Y: %e\n',cov_gyrY);
cov_gyrZ = cov(measure_gyr(3,:));
fprintf('covariance of Gyroscope in Z: %e\n',cov_gyrZ);

covz_gyr = [cov_gyrX;cov_gyrY;cov_gyrZ];

% Magnetometer
cov_magX = cov(measure_mag(1,:));
fprintf('covariance of Magnetometer in X: %e\n',cov_magX);
cov_magY = cov(measure_mag(2,:));
fprintf('covariance of Magnetometer in Y: %e\n',cov_magY);
cov_magZ = cov(measure_mag(3,:));
fprintf('covariance of Magnetometer in Z: %e\n',cov_magZ);

cov_mag = [cov_magX;cov_magY;cov_magZ];

% Plot the histograms
plothistogram(meas.acc,'Accelerometer')
plothistogram(meas.gyr,'Gyroscope')
plothistogram(meas.mag,'Magnetometer')


%% Task 5
% [xhat, meas] = filterTemplateTask5();

%% functions
function plotmeasurement(measure,t,name)
    % color
    color1 = '#0072BD';
    color2 = '#D95319';
    color3 = '#7E2F8E';

    k = length(t);
    figure()
    grid on
    hold on
    subplot(3,1,1)
    plot(t(:,1:k),measure(1,1:k),"Color",color1);
    legend('X')
    subplot(3,1,2)
    plot(t(:,1:k),measure(2,1:k),"Color",color2);
    legend('Y')
    subplot(3,1,3)
    plot(t(:,1:k),measure(1,1:k),"Color",color3);
    legend('Z')
    sgtitle(sprintf('%s', name));
end

function printm(value,name)
    fprintf('%s = \n',name);
    for i = 1:length(value)
        fprintf('%.4f\n', value(i));
    end
end

function plothistogram(measure,name)
    color1 = '#0072BD';
    color2 = '#D95319';
    color3 = '#7E2F8E';

    figure
    grid on
    hold on
    subplot(3,1,1)

    histogram(measure(1,:), 'Normalization', 'pdf','FaceColor',color1);
    title(sprintf('%s - X',name))
    subplot(3,1,2)
    histogram(measure(2,:), 'Normalization', 'pdf','FaceColor',color2);
    title(sprintf('%s - Y',name))
    subplot(3,1,3)
    histogram(measure(3,:), 'Normalization', 'pdf','FaceColor',color3);
    title(sprintf('%s - Z',name))
end