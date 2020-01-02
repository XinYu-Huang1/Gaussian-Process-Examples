% define tanh function 
% 

function y = tanh_func(x)
x1 = exp(x) - exp(-x);
x2 = exp(x) + exp(-x);
y = x1./x2;
end
