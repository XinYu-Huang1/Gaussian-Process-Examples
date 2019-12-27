% Prediction with noise free observations
% Gaussian Process with SE kernel
% K = gamma * exp(- |x1-x2|/(2*l^2))
% Copyright Zhidi Lin 

clc;clear;

l = 3;
sig_f = 1;

x_test = linspace(-10,10,4000)';
m_f = zeros(size(x_test));

% 4 observations (training points)
x_train = [-8 -3 -2 2 7];
fx_train = [5 1 0 2 -2];


% Calculate the partitions of the joint covariance matrix
[covXXInd1,covXXInd2] = meshgrid(x_train,x_train);
[covXXsInd1,covXXsInd2] = meshgrid(x_train,x_test);
[covXsXsInd1,covXsXsInd2] = meshgrid(x_test,x_test);
covXsXs = sig_f * exp(-(covXsXsInd1-covXsXsInd2).^2 ./ l.^2);
covXX = sig_f * exp(-(covXXInd1-covXXInd2).^2 ./ l.^2);
covXXs = sig_f * exp(-(covXXsInd1-covXXsInd2).^2 ./ l.^2);

% Cholesky decomposition of covariance matrix
chol_covXX = chol(covXX + 1e-9);


%posterior_mean = covXXs/covXX * fx_train';
posterior_mean = (covXXs/chol_covXX)/chol_covXX' * fx_train';
%posterior_cov = covXsXs - covXXs/covXX * covXXs';
posterior_cov = covXsXs - (covXXs/chol_covXX)/chol_covXX' * covXXs';

s2 = diag(posterior_cov);
f = [posterior_mean+2*sqrt(s2); flip(posterior_mean-2*sqrt(s2),1)];

% plot 
figure;
fill([x_test; flip(x_test,1)], f, [7 7 7]/8);
hold on; 
plot(x_train, fx_train, 'r+', 'MarkerSize', 12); hold on;  plot(x_train, fx_train, 'ro', 'MarkerSize', 12); hold on;
h1 = plot(x_test, posterior_mean,'k-', 'LineWidth', 2); 
axis([-10 10 -6 6])
grid on
xlabel('input, x')
ylabel('output, fx')
hold on;

% sampling from posterior 
K = chol(posterior_cov + 1e-10 * eye(size(posterior_cov)),'lower');   % Cholesky decomposition of covariance matrix

N = 3;                                                      % # of samples want to draw
C = linspecer(N);                                            % set color

for ii = 1:N
    uids = randn(size(posterior_mean));                                %  generate independent random numbers from standard normal distribution.
    f = posterior_mean' + transpose(uids) * K;                                     %  sample from existing
    h(ii) = plot(x_test',f,'color',C(ii,:),'linewidth',3); hold on;
    %scatter(x_test',f,[],C(ii,:),'filled');hold on;
    if ii == N
        hold off;
    end
end
legend([h1,h],'mean function','samples1','sample2','sample3')
title(' Prediction with noise free observations')
