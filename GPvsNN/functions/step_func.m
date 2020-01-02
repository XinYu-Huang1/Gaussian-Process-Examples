% define step function 
% 

function y = step_func(x)
 x(x>=0) = 1;
 x(x<0) = -1;
 y = x;
end