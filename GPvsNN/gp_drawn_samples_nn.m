
% Demo for drawing samples from a one-hidden neural network
% 
% Copyright Zhidi Lin

% close all;
clc; clear; close all;
addpath('functions')  

x = linspace(-15,15,1000)';  
tsize =size(x,1);

I = 1;                                       % dimensionality of input space 
H = 300;                                   %  # of Hidden unit 
mean_para = 0;
sigma_a = 1;                                % variance of bias a
sigma_u = 1;                                 % variance of weight parameters
sigma_b = 1;                                 % variance of bias b
omega_v = 1;

% variance of weight parameters v
sigma_v = omega_v * 1/(sqrt(H));            

% define one-layer hidden neural network with one-dimensional input data x
% and one-dimensional output 
a = normrnd(mean_para,sigma_a,1,H);
u = normrnd(mean_para,sigma_u,1,H);
b = normrnd(mean_para,sigma_b,1,1);
v = normrnd(mean_para,sigma_v,1,H);

actv = a + u .* x;
h = step_func(actv);   % hidden units

%% covariance out
mean = zeros(size(x));
C = covariance(h');
cov = sigma_b.^2 + omega_v.^2 .* C; 
K = chol(cov + 1e-9 * eye(size(cov)),'lower');   % Cholesky decomposition of covariance matrix

figure 
N = 2;                                             % # of samples want to draw
Color = linspecer(N);                                  % set color
for ii = 1:N
    uids = randn(size(mean));                     %  generate independent random numbers from standard normal distribution.
%     uids = sort(uids);
    f = mean +  K * uids;                           %  sample from existing
    plot(x,f,'color',Color(ii,:),'linewidth',2); hold on;
   % scatter(x,f,[],C(ii,:),'filled');hold on;
    if ii == N
        hold off;
    end
end
axis([-5 5 -6 6])
grid on
xlabel('input, x')
ylabel('output, fx')