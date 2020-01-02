% Demo for drawing samples from a one-hidden neural network
% 
% Copyright Zhidi Lin

close all;
clc; clear;
addpath('functions')  

% s = 32; randn('seed',s)                         % making sure the same random seed 

x = linspace(-5,5,1000)';  
tsize =size(x,1);

N = 2;                                           % # of samples want to draw
Color = linspecer(N);                            % set color

% parameters setup for simulations
I = 1;                                       % dimensionality of input space
H = 10000;                                   %  # of Hidden unit 
mean = 0;
sigma_a = 1;                                % variance of bias a
sigma_u = 1;                                 % variance of weight parameters
sigma_b = 1;                                 % variance of bias b
omega_v = 1;    
% variance of weight parameters v
sigma_v = omega_v * 1/(sqrt(H));   


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

        h{ii} = step_func(actv{ii});   % hidden units

        % Hidden to output 

        y_temp{ii} = v * h{ii}';

        y{jj,ii} = b + y_temp{ii};
    end
    aa = cell2mat(y);
    
    plot(x',aa(jj,:),'color',Color(jj,:),'linewidth',2); hold on;
end

grid on;

xlabel('x');ylabel('f(x)');

title(strcat(['Functions drawn from Brownian priors for NN (step hidden units) with ','H = '],num2str(H)));

saveas(gcf,strcat('figs/setp_H=',num2str(H),'.jpg'))

hold off;

set(gca,'fontsize',28,'linewidth',2);
hh = findobj('tag','legend');%|
set(hh,'fontsize',15) %| …Ë÷√legend◊÷∫≈¥Û–°