% define relu function 
% 
function y = relu_func(x,lamda)
 x(x>=0) = lamda .* (x(x>=0)).^2;
 x(x<0) = 0;
 y = x;
end