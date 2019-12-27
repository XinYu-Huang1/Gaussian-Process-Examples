% Demo for drawing samples from a Gaussian Process with RBF
% K = gamma * exp(- |x1-x2|/(2*l^2))
% Copyright Zhidi Lin 

clc;clear;

l = 5;                                            % lengthscale 
gamma = 2;                                        %  gamma = 1
x = linspace(-10,10,2000)';                       % points where the function will be evaluated 
mean = zeros(size(x));                            % Mean of GP (0 here)

[covd1,covd2] = meshgrid(x,x);                    % Generate all the possible pairs of points
cov = gamma * exp(-(covd1 - covd2).^2 ./ l.^2);   % compute the covariance function 
L = chol(cov + 1e-10 * eye(size(cov)),'lower');   % Cholesky decomposition of covariance matrix

N = 5;                                             % # of samples want to draw
C = linspecer(N);                                  % set color

% plot confidence zone
figure;
s1 = diag(cov);
f = [mean+3*sqrt(s1); flip(mean-3*sqrt(s1),1)]; 
fill([x; flip(x,1)], f, [7 7 7]/8); hold on;

for ii = 1:N
    uids = randn(size(mean));                     %  generate independent random numbers from standard normal distribution.
    f = mean +  L * uids;                           %  sample from existing
    plot(x,f,'color',C(ii,:),'linewidth',3); hold on;
   % scatter(x,f,[],C(ii,:),'filled');hold on;
    if ii == N
        hold off;
    end
end
axis([-10 10 -6 6])
grid on
xlabel('input, x')
ylabel('output, fx')
%title(strcat(['Sampling from prior ', 'l = '],num2str(l)));
set(gca,'fontsize',28,'linewidth',2);
legend('95% CR','samples1','samples2','samples3','samples4','samples5','samples6')
hh = findobj('tag','legend');%|
set(hh,'fontsize',15) %| ÉèÖÃlegend×ÖºÅ´óÐ¡
