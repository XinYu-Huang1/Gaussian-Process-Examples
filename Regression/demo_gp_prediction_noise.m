% Prediction with noisy observations
% Gaussian Process with SE kernel
% K = gamma * exp(- |x1-x2|/(2*l^2))
% Copyright Zhidi Lin 

clc;clear;
% close all;

l = 2.50;
gamma = 1;

sig_n = 0;

% 4 observations (training points)
x_train = [-6 -2 2 5];
y_train = [1 0 2 -1];

% test data
x_test = linspace(-10,10,100)';

% Calculate the partitions of the joint covariance matrix
[XX1,XX2] = meshgrid(x_train,x_train);
[XXs1,XXs2] = meshgrid(x_train,x_test);
[XsXs1,XsXs2] = meshgrid(x_test,x_test);

K_XsXs = gamma * exp(-(XsXs1-XsXs2).^2 ./ l.^2);
K_XXs = gamma * exp(-(XXs1-XXs2).^2 ./ l.^2);
K_XX = gamma * exp(-(XX1-XX2).^2 ./ l.^2);
K_XX_noisy = K_XX + sig_n^2 *eye(size(K_XX));


chol_KXX_noisy = chol(K_XX_noisy);
%posterior_mean = covXXs/covXX_noisy * fx_train';

posterior_mean = (K_XXs/chol_KXX_noisy)/chol_KXX_noisy' * y_train';
%posterior_cov = covXsXs - covXXs/covXX_noisy * covXXs';

posterior_cov = K_XsXs - (K_XXs/chol_KXX_noisy)/chol_KXX_noisy' * K_XXs';

s2 = diag(posterior_cov);
f = [posterior_mean+2*sqrt(s2); flip(posterior_mean-2*sqrt(s2),1)];

% plot 
figure;
h11 = fill([x_test; flip(x_test,1)], f, [7 7 7]/8);
hold on; 
h0 = plot(x_train, y_train, 'ok', 'linewidth',4, 'MarkerSize', 20, 'MarkerFaceColor', 'w');hold on;
% h0 = plot(x_train, y_train, 'xk', 'linewidth',3, 'MarkerSize', 20, 'MarkerFaceColor', 'k'); hold on;  
h1 = plot(x_test, posterior_mean,'k-', 'LineWidth', 4); 
axis([-10 10 -3 3])
grid on
xlabel('x')
ylabel('f(x)')
hold on;

% sampling from posterior 
K = chol(posterior_cov + 1e-12 * eye(size(posterior_cov)),'lower');   % Cholesky decomposition of covariance matrix

N = 3;                                                                % # of samples want to draw
% C = linspecer(N);                                                     % set color
C = ["-.k","-.k","-.k"];
% C = ["-ok","-sk","-^k"];
for ii = 1:N
    uids = randn(size(posterior_mean));                                %  generate independent random numbers from standard normal distribution.
    f = posterior_mean + K * uids;                                     %  sample from existing
%     h(ii) = plot(x_test',f,'--','color',C(ii,:),'linewidth',1); hold on;
    h(ii) = plot(x_test, f, C(ii),'linewidth',3, 'MarkerSize', 10, 'MarkerFaceColor', 'w'); hold on;
%     scatter(x_test',f,[],C(ii,:),'filled');hold on;
    if ii == N
        hold off;
    end
end
legend1 = legend([h0,h1,h(1),h11],'observations','mean function','samples','95% CR');
set(gca,'FontSize',45)
% title(strcat([' Prediction with noisy observations ', 'l = '],num2str(l)))
