clear all; close all;

%% BODE PLOTS
H = tf([-0.1,-2.4,-181,-1950],[1,3.3,990,2600]);
w = logspace(-2,3);
bode(H,w)
grid on

%% MACHINE LEARNING
figure;
%% Linear Regression
% Linear Regression (Generation)
a = 10;
b = -3.5;
x = 10*rand(100,1);
Y = a*x+b+5*randn(100,1);
scatter(x,Y)
% Linear Regression (Learning)
X = [ones(100,1) x];
L = X\Y
new_y = linspace(100,1)*L(2)+L(1)
hold on
plot(x,X*L)
title('Linear Regression')
hold off

%% Polynomial Regression
% Polynomial Regression (Generation)
figure;
n = 50;
a = 3;
b = 2*pi;
x = 1*rand(n,1);
Y = a*sin(b*x)+0.25*(randn(n,1));
subplot(131)
scatter(x,Y)
% Polynomial Regression (Learning)
X = [ones(n,1) bsxfun(@power,x,(1:24))];
L = X\Y
hold on
x = linspace(0,1,1000).';
X = [ones(1000,1), bsxfun(@power,x,(1:24))];
plot(linspace(0,1,1000),X*L)
title('Polynomial Regression')
axis([0 1 -5 5])
hold off

n = 50;
a = 3;
b = 2*pi;
X = 1*rand(n,1);
Y = a*sin(b*X)+0.25*(randn(n,1));
X = [ones(n,1) bsxfun(@power,X,(1:24))];
subplot(133)
scatter(X(:,2),Y)
hold on
%lasso (L1) regression
[B1, FitInfo1] = lasso(X,Y,'Alpha',1,'Lambda',0.0005); 
 %ridge (L2) regression
[B, FitInfo] = lasso(X,Y,'Alpha',1e-20,'Lambda',0.0005);
x = linspace(0,1,n).';
x = [ones(n,1), bsxfun(@power,x,(1:24))];
y = x*B;
plot(x(:,2),y)
title('Ridge (L2) Regression')
axis([0 1 -5 5])
hold off

subplot(132)
scatter(X(:,2),Y)
hold on
y = x*B1;
plot(x(:,2),y)
title('Lasso (L1) Regression')
axis([0 1 -5 5])
hold off

%% Cross Validation 
figure;
% For choosing parameters
n = 50;
a = 3;
b = 2*pi;
x = 1*rand(n,1);
y = a*sin(b*x)+0.25*(randn(n,1));
scatter(x,y)

se = zeros(3,1);

ind = randperm(n);
cv_ind = reshape(ind, 5, 10);
lambda_vec = [1e-5 1e-4 1e-3];
for i = 1:length(lambda_vec)
    for j = 1:10
       training_mtx = cv_ind(1:end,1:9);
       testing_vec = cv_ind(1:end,10);
       y_train = y(training_mtx(:));
       x_train = x(training_mtx(:));
       x_train = [ones(45,1) bsxfun(@power,x_train,(1:24))];
       x_test = x(testing_vec);
       x_test = [ones(5,1) bxsfun(@power,x_test,(1:24))];
       y_test = y(testing_vec);
       
       [B, FitInfo1] = lasso(x_train,y_train,'Alpha',1,'Lambda',lambda_vec(i));
       
       y_estimate = FitInfo.Intercept + x_test*B;
       se(i) = se(i) + sum((y_estimate-y_test).^2);
       cv_ind = circshift(cv_ind,1,2);
    end
end
mse = se/50;
[min_val,min_ind] = min(mse);
lambda_winner = lambda_vec(min_ind)
[B, FitInfo1] = lasso([ones(50,1), x.^(1:24)],y,'Alpha',1,'Lambda',lambda_winner);
X = linspace(0,1,1000);
plot(X,[ones(1000,1), (X.').^(1:24)]*B)
hold on
scatter(x,y)
hold off
