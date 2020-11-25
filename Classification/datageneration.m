function r = datageneration(n)

if nargin == 0
    n = 10;
end

mu1 = [0.2 0.5]; Sigma1 = [.09 .04; .04 .03];
r1 = mvnrnd(mu1, Sigma1, n); 
scatter(r1(:,1),r1(:,2), 50, 'filled', 'o');
hold on;
mu2 = [0.7 0.7]; Sigma2 = [.09 .04; .04 .03];
r2 = mvnrnd(mu2, Sigma2, n); 
scatter(r2(:,1),r2(:,2), 50, 'o');

r = [r1;r2];
end