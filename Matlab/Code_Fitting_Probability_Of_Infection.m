clear all;
close all;
clc;

% Data
ViralData = [3.5 4.2 5.3 5.4 6.3 6.4 6.5 7.4 7.9];  % log Scale
Probability_Infection = [0 0  1/17 0 5/30 11/30 5/17 19/26 29/31];

% Bounds for parameters
lb = [1e-5 0];
ub = [1e-2 10];

% Initial guess for parameters
k = [0.001  5.2785];


options = optimset('Display', 'iter', 'TOLX', 1e-12, 'TOLFun', 1e-12, 'MaxIter', 10000);

[k, fval] = fmincon(@(k) err_in_data(k), k, [], [], [], [], lb, ub, [], options);

% Extract optimized parameters
a = k(1);
h = k(2);


% Model prediction
V = 0:0.1:9;
Model_Prob_Infection = 1 - exp(-a*(V-3.5).^h);
% Model_Prob_Infection = h./(1+exp(-a*(V-b)));   %fun 3

% AIC calculation
data_points = length(ViralData);
KK = 2;  % Number of parameters
AIC = data_points*log(fval/data_points) + 2*KK + (2*KK*(KK+1))/(data_points - KK -1);

% Plot the data and model
figure;
plot(ViralData, Probability_Infection, 'r.', 'MarkerSize', 25);
hold on;
plot(V, Model_Prob_Infection, 'b', 'LineWidth', 2);
xlim([2, 9]);
ylim([0, 1]);
set(gca, 'FontSize', 15, 'FontName', 'Arial', 'LineWidth', 1, 'FontWeight', 'bold');
set(gca, 'YGrid', 'on', 'XGrid', 'off', 'LineWidth', 1, 'FontSize', 14);
xlabel('Viral Load (log scale)', 'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold');
ylabel('Probability of Infection', 'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold');

% Error function for optimization
function error_in_data = err_in_data(k)
   ViralData = [3.5 4.2 5.3 5.4 6.3 6.4 6.5 7.4 7.9];  % log Scale
Probability_Infection = [0 0  1/17 0 5/30 11/30 5/17 19/26 29/31];

    a = k(1);
    h = k(2);
    

    % Model prediction
    Model_Prbobability = 1 - exp(-a*(ViralData-3.5).^h);  % fun 2
    

    % Error (sum of squared residuals)
    error_in_data = sum((Model_Prbobability - Probability_Infection).^2);
end
