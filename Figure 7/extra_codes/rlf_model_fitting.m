function [data_model, r_square, threshold, spls_interp, data_interp] = rlf_model_fitting(data)
% data is the rlf function
% i = 3
% data = wt.rlf(i,:);
spls = 0:10:70; % sound levels
spls_interp = 0:1:70; % interpolation sound levels
% data_interp = interp1(spls,data,spls_interp, 'cubic');
data_interp = interp1(spls,data,spls_interp);

[~, u] = max(data_interp),

y = @(x,level) (x(1) * exp(-(level - x(2)).^2/(2 * x(3).^2)) + x(4)).*(level<=x(2)) +  ( x(1) * exp(-(level - x(2)).^2/(2 * x(5).^2)) + x(6)).*(level>x(2))
x0 = [1, u, 10, 1, 30, 1];
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
% lb = [0,0,0,0,0,0];
lb = [];
ub =[];
% ub = [100, 60, 100, 100, 100, 100];
x = lsqcurvefit(y,x0,spls_interp,data_interp, lb, ub);
data_model = y(x, spls_interp);
r_square = 1 - sum((y(x, spls_interp) - data_interp).^2)/(sum((data_interp - mean(data_interp)).^2));
threshold = x(2);
% figure; 
% plot(spls, data)
% hold on
% plot(spls_interp, y(x, spls_interp))
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
%% try to use fmincon to fit the rlf
% fun = @(x) rlf_modelCost(x, spls_interp, data_interp);
% % options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
% [~, u] = max(data_interp);
% x0 = [1, u, 10, 1, 30, 1];
% % x = fmincon(fun, x0, [], [], [], [], [0,0,0,0,0,0], [100,70,100,100,100,100], [],options);
% x = fmincon(fun, x0, [], [], [], [], [0,0,0,0,0,0], [100,70,100,100,100,100]);
% 
% y = @(x,level) (x(1) * exp(-(level - x(2)).^2/(2 * x(3).^2)) + x(4)).*(level<=x(2)) +  ( x(1) * exp(-(level - x(2)).^2/(2 * x(5).^2)) + x(6)).*(level>x(2))
% figure; 
% plot(spls_interp, data_interp)
% hold on
% plot(spls_interp, y(x, spls_interp))
