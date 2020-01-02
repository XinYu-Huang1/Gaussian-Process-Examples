% Demo for drawing samples from a one-hidden neural network
% 
% Copyright Zhidi Lin

close all;
clc; clear;
addpath('functions')  

s = 35; randn('seed',s)                         % making sure the same random seed 

x = linspace(-5,5,1000)';  

tsize =size(x,1);

N = 2;                                           % # of samples want to draw
Color = linspecer(N);                            % set color

eta = 3;                                   % parameter for signz activiation function 
I = 1;                                       % dimensionality of input space
H = 10000;                                   %  # of Hidden unit 
mean = 0;

omega_a = 1;                                 % omega of bias a
omega_u = 1;                                 % omega of weight parameters
sigma_b = 1;                                 % variance of bias b
omega_v = 1;    
% variance of weight parameters v
sigma_v = omega_v * 1/(sqrt(H));     

% gaussian distribution of a & u
x_par = 0.01:0.1:20;
A_j = pdfrnd(x_par, invgamdis(x_par,eta), H);
for ii = 1:H
    sigma_a(ii) = A_j(ii) * omega_a;
    a(ii) = normrnd(mean,sigma_a(ii),1,1);
    
    sigma_u(ii) = A_j(ii) * omega_u;
    u(ii) = normrnd(mean,sigma_u(ii),1,1);
end

figure

for jj = 1:N  
    % define one-layer hidden neural network with one-dimensional input data x
    % and one-dimensional output 

    b = normrnd(mean,sigma_b,1,1);

    v = normrnd(mean,sigma_v,1,H); 

    % construct a neural network

    for ii = 1:tsize
        % Input to hidden layer 
        % activiation function 
        actv{ii} = a + u .* x(ii);

        h{ii} = tanh_func(actv{ii});   % hidden units

        % Hidden to output 

        y_temp{ii} = v * h{ii}';

        y{jj,ii} = b + y_temp{ii};
    end
    aa = cell2mat(y);
    
    plot(x',aa(jj,:),'color',Color(jj,:),'linewidth',2); hold on;
end

grid on;

xlabel('x');ylabel('f(x)');

title(strcat(['Functions drawn from fractional Brownian priors for NN (H = 10000 tanh units) with ','\eta = '],num2str(eta)));

saveas(gcf,strcat('figs/tanh_ita=',num2str(eta),'.jpg'))

hold off;