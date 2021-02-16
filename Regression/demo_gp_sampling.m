% Demo for drawing samples from a Gaussian Process with RBF
% 
% K = gamma * exp(- |x1-x2|/(l^2))
% Copyright Zhidi Lin 

clc; clear; close all

% lengthscale 
l = 2.5;                                         

% gamma = 1
gamma = 1;                                    

% points where the function will be evaluated 
x = linspace(-10,10,100)';                        

% Mean of GP (0 here)
m = zeros(size(x));                           

% Generate all the possible pairs of points
[covd1,covd2] = meshgrid(x,x);                   

% compute the covariance function 
cov = gamma * exp(-(covd1 - covd2).^2 ./ l.^2);   

N = 3;                                             % # of samples want to draw
% C = linspecer(N);                                % set color
% C = ["-ok","-sk","-^k"];
C = ["-.k","-.k","-.k", "-.k", "-.k"];
% mc = {'w','k',[90 90 90]/255};

%% plot confidence zone
figure;
s1 = diag(cov);
f = [m+2*sqrt(s1); flip(m-2*sqrt(s1),1)]; 
h0 = fill([x; flip(x,1)], f, [7 7 7]/8); hold on;

% Cholesky decomposition of covariance matrix
K = chol(cov + 1e-12 * eye(size(cov)), 'lower');  

for ii = 1:N
    %  generate independent random numbers from standard normal distribution.
    uids = randn(size(m));                   
    
    %  sample from existing
    f = m +  K * uids;        
    
    h(ii) = plot(x,f, C(ii),'linewidth',4, 'MarkerSize', 10, 'MarkerFaceColor', 'k'); hold on;
    if ii == N
        hold off;
    end
end

axis([-10 10 -3 3])
xlabel('x')
ylabel('f(x)')
legend1 = legend([h(1), h0], 'samples', '95% CR'); %,'sample 2','sample 3');
set(gca,'FontSize',45)
grid on;
