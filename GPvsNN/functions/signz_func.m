% define signz function 
% 

function y = signz_func(z,ita)
temp = (ita-1)/2;
y = sign(z).*abs(z).^temp;
end