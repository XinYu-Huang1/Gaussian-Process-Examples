%%invgamdis(A,eta): function returning a inverse gamma pdf sampled at x
%%coords, gamma given by a/k and b/theta, depending on your notation.

function p = invgamdis(x,eta)
% return the invgam distribution with parameters x,ita
% p = invgamdis(x,eta);

if nargin<2
   error('Requires three input arguments.'); 
end

p = (x.^(eta)).*exp(-(eta -1)./ (2 .* x.^2));
