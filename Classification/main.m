% demo for Laplace approximation on a 2-d classification task using GP
%
% Author: Zhidi Lin
% Reference: GPML code
%
% 2020-10-31

clc; clear all; close all;
%% data preparation: 10 data for each class
% One can use datageneration function to generate simulated data

load('data.mat');

n = length(r)/2;
y = [ones(n,1); -1 * ones(n,1)];
scatter(r(1:n,1),r(1:n,2), 100, 'filled', 'o');
hold on;
scatter(r(1+n:end,1),r(1+n:end,2), 100, 'filled','o'); 
grid on; 
hold on;

%% Kernel for demonstration
l = 0.15;              % lengthscale 
sigma_f = 9;           % sigma_f

sqdist = sq_dist(r',r');         % pairwise squared distances distance 
K = sigma_f * exp(- sqdist ./ (2 * l.^2));    % the covariance function 

%% Newton's iterations
tol = 1e-6;                   % tolerance for when to stop the Newton iterations

Psi_old = -Inf;               % make sure while loop starts
Psi_new = 0;
alpha = 0;
f = zeros(2*n,1);

% [lp,dlp,d2lp,d3lp] = logistic(y, f, 'deriv');
[lp,dlp,d2lp,d3lp] = cumGauss(y, f, 'deriv');

W = -d2lp;

while Psi_new - Psi_old > tol            % begin Newton's iterations
  Psi_old = Psi_new; 
  alpha_old = alpha; 
  
  sW = sqrt(W);                     
  L = chol(eye(2*n)+sW * sW'.* K);       % L'*L=B=I+sW*K*sW
  b = W.*f+dlp;
  temp = L\(L'\(sW.*(K*b)));
  alpha = b - sW .* temp;
  f = K * alpha;
%   [lp,dlp,d2lp,d3lp] = logistic(y, f, 'deriv');
  [lp,dlp,d2lp,d3lp] = cumGauss(y, f, 'deriv');

  W=-d2lp;

  Psi_new = -alpha'*f/2 + lp;
  
end  % end Newton's iterations

sW = sqrt(W);                                  % recalculate L
L  = chol(eye(2 * n)+sW * sW'.*K);             % L'*L=B=I+sW*K*sW
nlZ = alpha'*f/2 - lp + sum(log(diag(L)));     % approx neg. log marg likelihood

%% prediction
[t1, t2] = meshgrid(0:0.02:1, 0:0.02:1);

% data for testing
t = [t1(:) t2(:)]; 

sqdist_ss = sq_dist(t',t');                        % pairwise squared distances distance 
kstarstar = sigma_f * exp(- sqdist_ss ./ (2 * l.^2));    % the covariance function 

sqdist_s = sq_dist(r',t');   
kstar = sigma_f * exp(- sqdist_s ./ (2 * l.^2));

mu = kstar' * alpha;                                % predictive means
v  = L'\(repmat(sW, 1, length(t1(:))) .*  kstar);
s2 = kstarstar - sum(v.*v,1)';                      % predictive variances

% p = logistic([], mu, diag(s2));                   % predictive probabilities
p = cumGauss([], mu, diag(s2));

%% plot the results
[cs, h] = contour(t1, t2, reshape(p, size(t1)), [0:0.3:1,0.5], 'linewidth',1);
clabel(cs, h, 'FontSize', 15);
legend('class +1','class -1')
set(gca,'FontSize',22)
