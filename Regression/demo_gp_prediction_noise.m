% Prediction with noisy observations
% Gaussian Process with SE kernel
% K = gamma * exp(- |x1-x2|/(2*l^2))
% Copyright Zhidi Lin 

clc;clear;

l = 1.5;
sig_f = 1;

sig_n = 0.1;

x_test = (-10:0.05:10)';
m_f = zeros(size(x_test));

% 4 observations (training points)
x_train = [-3 -2 2 7];
y_train = [1 0 2 -2];

% Calculate the partitions of the joint covariance matrix
[covXXInd1,covXXInd2] = meshgrid(x_train,x_train);
[covXXsInd1,covXXsInd2] = meshgrid(x_train,x_test);
[covXsXsInd1,covXsXsInd2] = meshgrid(x_test,x_test);
K_test = sig_f * exp(-(covXsXsInd1-covXsXsInd2).^2 ./ l.^2);
K_cross = sig_f * exp(-(covXXsInd1-covXXsInd2).^2 ./ l.^2);
K_train = sig_f * exp(-(covXXInd1-covXXInd2).^2 ./ l.^2);
C = K_train + sig_n^2 *eye(size(K_train));


L = chol(C);
%posterior_mean = covXXs/covXX_noisy * fx_train';
posterior_mean = (K_cross/L)/L' * y_train';
%posterior_cov = covXsXs - covXXs/covXX_noisy * covXXs';
posterior_cov = K_test - (K_cross/L)/L' * K_cross';

s2 = diag(posterior_cov);
f = [posterior_mean+2*sqrt(s2); flip(posterior_mean-2*sqrt(s2),1)];

% plot 
figure;
fill([x_test; flip(x_test,1)], f, [7 7 7]/8);
hold on; 
plot(x_train, y_train, 'r+', 'MarkerSize', 12); hold on;  plot(x_train, y_train, 'ro', 'MarkerSize', 12);hold on;
h1 = plot(x_test, posterior_mean,'k-', 'LineWidth', 2); 
axis([-10 10 -6 6])
grid on
xlabel('input, x')
ylabel('output, fx')
hold on;

% sampling from posterior 
K = chol(posterior_cov + 1e-10 * eye(size(posterior_cov)),'lower');   % Cholesky decomposition of covariance matrix

N = 2;                                                                % # of samples want to draw
C = linspecer(N);                                                     % set color

for ii = 1:N
    uids = randn(size(posterior_mean));                                %  generate independent random numbers from standard normal distribution.
    f = posterior_mean' + transpose(uids) * K;                                     %  sample from existing
    h(ii) = plot(x_test',f,'color',C(ii,:),'linewidth',3); hold on;
%     scatter(x_test',f,[],C(ii,:),'filled');hold on;
    if ii == N
        hold off;
    end
end
legend([h1,h],'mean function','smaple1','sample2')
title(strcat([' Prediction with noisy observations ', 'l = '],num2str(l)))
