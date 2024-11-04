% HA1
clear
clc
close all
%% 1c
% Parameters
mu = [0; 10];

Sigma = [0.3, 0; 
           0, 8];

level = 3;

npoints = 150;

A = [1, 0.5;
     0,   1];

% for q
[ xy_q ] = sigmaEllipse2D(mu, Sigma, level, npoints);

% for z
[mu_z, Sigma_z] = affineGaussianTransform(mu, Sigma, A, 0);

[ xy_z ] = sigmaEllipse2D(mu_z, Sigma_z, level, npoints);

% Plot
figure();
h1 = plot(xy_q(1,:), xy_q(2,:));
axis equal
hold on
grid on
plot(mu(1), mu(2), '*', 'color', h1.Color);

h2 = plot(xy_z(1,:), xy_z(2,:));
axis equal
hold on
grid on
plot(mu_z(1), mu_z(2), '*', 'color', h2.Color);
legend('3\sigma ellipse for q', 'mean q','3\sigma ellipse for z','mean z')
xlabel('x_1'); 
ylabel('x_2');
% print('1c.eps','-depsc')

%% 2a
mu_x = 0;
Sigma_x = 2;
N = 10000;
f = @(x) 3*x;
x = [-20:0.5:20];

%  analytically
[mu_z1, Sigma_z1] = affineGaussianTransform(mu_x, Sigma_x, 3, 0);
anal_pdf = normpdf(x,mu_z1, sqrt(Sigma_z1));

% numerical Gaussian approximation
[mu_z2, Sigma_z2, y_s] = approxGaussianTransform(mu_x, Sigma_x, f, N);
approx_pdf = normpdf(x,mu_z2, sqrt(Sigma_z2));

% PLOT
figure()
histogram(y_s,x,'Normalization','pdf')
hold on
grid on
plot(x,anal_pdf,'LineWidth',2)
legend('z', ' analytical-pdf')
% print('2a1.eps','-depsc')

figure()
histogram(y_s,x,'Normalization','pdf')
hold on
grid on
plot(x,approx_pdf,'LineWidth',2)
legend('z', 'approximation-pdf')
% print('2a2.eps','-depsc')

%% 2b 
f_2b = @(x) x.^3;
x_2b = [-40:1:40];

%  analytically
p_2b = @(x) x.^6 .* exp(-x.^2 / 4) / sqrt(4 * pi);
Var = integral(p_2b, -Inf, Inf);
fprintf('Cov[z] = %.6f\n',Var);
anal_pdf_2b = normpdf(x_2b,0, sqrt(Var));

% numerical Gaussian approximation
[mu_z_2b, Sigma_z_2b, y_s_2b] = approxGaussianTransform(mu_x, Sigma_x, f_2b, N);
approx_pdf_2b = normpdf(x_2b,mu_z_2b, sqrt(Sigma_z_2b));

% plot
figure()
histogram(y_s_2b,x_2b,'Normalization','pdf')
hold on
grid on
plot(x_2b,anal_pdf_2b,'LineWidth',1.5)
legend('z', ' analytical-pdf')
% print('2b1.eps','-depsc')

figure()
histogram(y_s_2b,x_2b,'Normalization','pdf')
hold on
grid on
plot(x_2b,approx_pdf_2b,'LineWidth',1.5)
legend('z', 'approximation-pdf')
% print('2b2.eps','-depsc')

%% 3e
% linear function
f_linear =  @(x) 6.*x;

% nonlinear function 
f_nonlinear =  @(x) x.^3;

%number of sample
N_3e = 100000;

% Gaussian distribution for r, N(0,9)
r_sample =  mvnrnd(0, 9, N_3e)';

% Gaussian distribution for x, N(1,9)
gaussian_sample =  mvnrnd(1, 4, N_3e)';

% uniform distribution for x, U(2,8)
uniform_sample = (8 - 2) * rand(1, N_3e) +2;

%% linear transforamtion with uniform distribution
y_u_linear= f_linear(uniform_sample) + r_sample;

y_1 = [0:0.5:60];
% analytically
anal_pdf_1=[];
for i= 1:length(y_1)
    % define the function
    fun1 = @(x)((1/6)*(1/sqrt(18*pi))*exp(-(y_1(:,i)-f_linear(x)).^2/18));
    int1 =integral(fun1,2,8);
    anal_pdf_1=[anal_pdf_1,int1];
end

% Plot 
figure()
histogram(y_u_linear,y_1,'Normalization','pdf')
hold on
grid on
plot(y_1,anal_pdf_1,'LineWidth',2)
legend('histogram of y', 'analytical-pdf')
title('linear transforamtion with uniform distribution')
% print('3e1.eps','-depsc')

%% nonlinear transforamtion with uniform distribution
y_u_nonlinear= f_nonlinear(uniform_sample) + r_sample;

y_2 = [0:5:600];
% analytically
anal_pdf_2=[];
for i= 1:length(y_2)
    % define the function 
    fun2 = @(x)((1/6)*(1/sqrt(18*pi))*exp(-(y_2(:,i)-f_nonlinear(x)).^2/18));
    int2 =integral(fun2,2,8);
    anal_pdf_2=[anal_pdf_2,int2];
end

% Plot 
figure()
histogram(y_u_nonlinear,y_2,'Normalization','pdf')
hold on
grid on
plot(y_2,anal_pdf_2,'LineWidth',2)
legend('histogram of y', 'analytical-pdf')
title('nonlinear transforamtion with uniform distribution')
% print('3e2.eps','-depsc')
%% linear transforamtion with Gaussian distribution
y_n_linear= f_linear(gaussian_sample) + r_sample;

y_3 = [-60:1:80];

% analytically
mu_y3 = 6;
Sigma_y3 = 36*4+9;
anal_pdf_3 = normpdf(y_3,mu_y3,sqrt(Sigma_y3));

% Plot 
figure()
histogram(y_n_linear,y_3,'Normalization','pdf')
hold on
grid on
plot(y_3,anal_pdf_3,'LineWidth',2)
legend('histogram of y', 'analytical-pdf')
title('linear transforamtion with Gaussian distribution')
% print('3e3.eps','-depsc')

%% nonlinear transforamtion with Gaussian distribution
y_n_nonlinear= f_nonlinear(gaussian_sample) + r_sample;

y_4 = [-200:5:300];

% analytically
anal_pdf_4=[];
for i= 1:length(y_4)
    %define the function that will be integralled
    fun4 = @(x)((1/sqrt(8*pi)).*exp(-(x-1).^2/8).*(1/sqrt(18*pi)).*exp(-(y_4(:,i)-f_nonlinear(x)).^2/18));
    %get the function of p(y) and the value of p(y with)
    int4 =integral(fun4,-inf,+inf);
    anal_pdf_4=[anal_pdf_4,int4];
end

% Plot 
figure()
histogram(y_n_nonlinear,y_4,'Normalization','pdf')
hold on
grid on
plot(y_4,anal_pdf_4,'LineWidth',2)
legend('histogram of y', 'analytical-pdf')
title('nonlinear transforamtion with Gaussian distribution')
% print('3e4.eps','-depsc')

%% p(y|x) linear function with Gaussian distribution
gaussian_sample_g = mvnrnd(1, 4, 1);
y_ng_linear= f_linear(gaussian_sample_g) + r_sample;

y_5 = [-100:1:100];
% analytically
Sigma_y5 = 9;
anal_pdf_5 = normpdf(y_5, f_linear(gaussian_sample_g),sqrt(Sigma_y5));

% Plot 
figure()
histogram(y_ng_linear,'Normalization','pdf')
hold on
grid on
plot(y_5,anal_pdf_5,'LineWidth',2)
legend('histogram of y', 'analytical-pdf')
title('p(y|x) for linear function with Gaussian distribution')
% print('3e5.eps','-depsc')

%% p(y|x) nonlinear function with Gaussian distribution
y_ng_nonlinear= f_nonlinear(gaussian_sample_g) + r_sample;

y_6 = [-100:1:100];
% analytically
Sigma_y6 = 9;
anal_pdf_6 = normpdf(y_6,f_nonlinear(gaussian_sample_g),sqrt(Sigma_y6));

% Plot 
figure()
histogram(y_ng_nonlinear,y_6,'Normalization','pdf')
hold on
grid on
plot(y_6,anal_pdf_6,'LineWidth',2)
legend('histogram of y', 'analytical-pdf')
title('p(y|x) for nonlinear function with Gaussian distribution')
% print('3e6.eps','-depsc')

%% p(y|x) linear function with uniform distribution
uniform_sample_g = (8 - 2) * rand +2;
y_ug_linear= f_linear(uniform_sample_g) + r_sample;

y_7 = [-100:1:100];
% analytically
Sigma_y7 = 9;
anal_pdf_7 = normpdf(y_7,f_linear(uniform_sample_g),sqrt(Sigma_y7));

% Plot 
figure()
histogram(y_ug_linear,y_7,'Normalization','pdf')
hold on
grid on
plot(y_7,anal_pdf_7,'LineWidth',2)
legend('histogram of y', 'analytical-pdf')
title('p(y|x) for linear function with uniform distribution')
% print('3e7.eps','-depsc')

%% %% p(y|x) nonlinear function with uniform distribution
uniform_sample_g = (8 - 2) * rand +2;
y_ug_nonlinear= f_nonlinear(uniform_sample_g) + r_sample;

y_8 = [-100:1:100];
% analytically
Sigma_y8 = 9;
anal_pdf_8 = normpdf(y_8,f_nonlinear(uniform_sample_g),sqrt(Sigma_y8));

% Plot 
figure()
histogram(y_ug_nonlinear,y_8,'Normalization','pdf')
hold on
grid on
plot(y_8,anal_pdf_8,'LineWidth',2)
legend('histogram of y', 'analytical-pdf')
title('p(y|x) for nonlinear function with uniform distribution')
% print('3e8.eps','-depsc')

%% 4a
% Drawing samples
Sigma_4a = 0.5^2;

% number of samples
n = 100000;
y = zeros(1,n);

% sample w
w_sample = mvnrnd(0,Sigma_4a,n);

% sample theta
p = [0.5, 0.5];
    
theta_sample = randsample([-1, 1], n, true, p);

% add them up
for i = 1:n
    y(i) = theta_sample(i) + w_sample(i);
end

x_4a = [-4:0.05:4];

% analytical-pdf
anal_pdf_4a = 0.5*normpdf(x_4a, -1 , 0.5) + 0.5*normpdf(x_4a, 1 , 0.5);

% Plot
figure()
histogram(y,x_4a,'Normalization','pdf')
hold on
grid on
plot(x_4a,anal_pdf_4a,'LineWidth',2)
legend('histogram of y', 'analytical-pdf')
% print('4a.eps','-depsc')

%% 4b
1/(sqrt(2*pi)*0.5)*exp(-(0.7-1)^2/(2*0.5^2))

1/(sqrt(2*pi)*0.5)*exp(-(0.7+1)^2/(2*0.5^2))

tanh(2.8)

%% How to export your source code as .txt file. 
filename = fullfile('HA1.m');  % You should change "mainExample.m" with the name of your source file!
copyfile(filename,'HA1.txt','f') % Here, 'mainExample.txt' is the output. You should upload the 'main.txt' (or whatever you name it).