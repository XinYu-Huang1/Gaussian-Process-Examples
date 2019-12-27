% Demo for drawing samples from a one-hidden neural network
% 
% Copyright Zhidi Lin

close all;
clc; clear;

addpath('functions')  

% s = 35; randn('seed',s)                         % making sure the same random seed 

x = linspace(-5,5,1000)';  

tsize =size(x,1);

N = 2;                                           % # of samples want to draw

Color = linspecer(N);                            % set color

% parameters setup for simulations
I = 1;                                       % dimensionality of input space
H = 10000;                                   %  # of Hidden unit 
mean = 0;
sigma_a = 5;                                % variance of bias a
sigma_u = sigma_a * 1;                       % variance of weight parameters
sigma_b = 1;                                 % variance of bias b
omega_v = 1;    
sigma_v = omega_v * 1/(sqrt(H));             % variance of weight parameters v
figure

for jj = 1:N         

    % define one-layer hidden neural network with one-dimensional input data x
    % and one-dimensional output 

    a = normrnd(mean,sigma_a,1,H);

    u = normrnd(mean,sigma_u,1,H);

    b = normrnd(mean,sigma_b,1,1);

    v = normrnd(mean,sigma_v,1,H); 

    % construct a neural network

    for ii = 1:tsize
        % Input to hidden layer 
        % activiation function 
        actv{ii} = a + u .* x(ii);

        h{ii} = sigmoid_func(actv{ii});   % hidden units

        % Hidden to output 

        y_temp{ii} = v * h{ii}';

        y{jj,ii} = b + y_temp{ii};
    end
    aa = cell2mat(y);
    
    plot(x',aa(jj,:),'color',Color(jj,:),'linewidth',2); hold on;
end

grid on;

xlabel('x');ylabel('f(x)');

title(strcat(['Functions drawn from soomth priors for NN (sigmoid hidden units) with ','\sigma_u = '],num2str(sigma_u)));

saveas(gcf,strcat('figs/tanh_siamg_u=',num2str(sigma_u),'.jpg'))

hold off;


% define sigmoid function 
% 

function y = sigmoid_func(x)
 y = 1./(1+exp(-x));
end
