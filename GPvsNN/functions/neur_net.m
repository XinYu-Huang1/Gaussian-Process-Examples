% define a neural network 
% train_data should be a row matrix
function y = neur_net(weight_u, weight_v, train_data, bias1, bias2)

n = size(train_data,2);
mu = size(weight_u,2);


%  input-to-hidden layer 
train_data1= repmat(train_data',1,mu);
yy = (train_data1) .* weight_u  + bias1;
Hi = step_func(yy);
% Hi = tanh(yy);

%  hidden-to-output layer 
y =  Hi * weight_v + bias2;
end