%within-host fitting
clear all
close all
clc

Viral_load_data =  [1.699	2.652	2.678	5.202	4.248	6.459	4.991	4.153	4.188	2.722	3.014	1.699]'; %log

time_Viral =[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5 , 6]'; % days


 params = [1.00432067935374e-05	6.99890940630983	19.4974816748663	10.5174125251227	3.48621792230982];
 
 lb = [1e-10 0.5  0  0     0];
 ub = [7.4    7  25  1e8  50];

opts_con = optimset('Display','iter', 'TolX', 1e-12);


[params, fval] =  fmincon(@err_in_data,params,[],[],[],[],lb,ub,[],opts_con);

dtau = 0.01;
infection_timeforward = 0:dtau:7;
initial_condition = [4e6 0 0 1];

 within_params = params;

[t_r, y_r] = ode15s(@(t,y)Usutu_inHost(y,within_params),infection_timeforward, initial_condition);

figure;
plot(infection_timeforward, log10(y_r(:,4)), 'b-', 'LineWidth', 2); % model solution Sh
hold on;
plot(time_Viral, Viral_load_data, 'ro', 'MarkerSize', 8, 'LineWidth', 2); % observed data
xlabel('Time (days)');
ylabel('Viral Load');

 [max_val, idx] = max(y_r(:,4))
idx = find(y_r(:,4) < 10)


function error_in_data = err_in_data(k)


    Viral_load_data =  [1.699	2.652	2.678	5.202	4.248	6.459	4.991	4.153	4.188	2.722	3.014	1.699]'; %log

    time_Viral =[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5 , 6]'; % days

    initial_condition = [4e6 0 0 1];
    sol = ode15s(@(t,y)Usutu_inHost(y,k),[0 7], initial_condition);
    V_model = deval(sol, time_Viral, 4)';
    Model_Viral = log10(V_model);

  
   error_in_data = sum((Model_Viral - Viral_load_data).^2) ;
   
   


end




function dy = Usutu_inHost(y,params)

    dy = zeros(4, 1);

    beta = params(1);
    d = params(2);
    delta = params(3);
    pi = params(4);
    c = params(5);

        
        T = y(1);
        E = y(2);
        I = y(3);
        V = y(4);

        dy(1) = -beta * V * T;
        dy(2) = beta * V * T - d * E;
        dy(3) = d * E - delta * I;
        dy(4) = pi * I - c * V;
 end