function cov = covariance(X)
% 由随机变量样本矩阵计算协方差矩阵
%---- 输入：    
% X:     M*N的样本矩阵，其中一行表示一个随机向量样本
%       共有M个随机样本
%---- 输出：
% cov：   N*N的协方差矩阵，表示各个随机变量的协方差
%---- 计算方法
% 随机变量均值用样本均值统计量估计：X_mean = 1/N*ΣXi;
% 随机变量方差用样本方差统计量估计：S = 1/(N-1)*Σ(Xi-X_mean)^2
% 随机变量协方差可以用如下统计量估计：C = 1/(N-1)*Σ(Xi-X_mean)(Yi-Y_mean)

% 各个样本减去其平均值
% X = bsxfun(@minus,X,mean(X));
%或者：
%X = X-repmat(mean(X),size(X,1),1);
cov = 1/(size(X,1)-1)*X'*X;