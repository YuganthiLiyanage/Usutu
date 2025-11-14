clear all
close all
clc


% Data: infected birds & susceptible birds & Viral load

Infected_birds_data = [0.00655021834, 0.03766375546, 0.1359170306, 0.1517467249, 0.2194323144, 0.09497816594, 0.00545851528, 0.01965065502, 0.03220524017, 0]';
Susceptible_birds_data = [0.99, 0.612408962, 0.725485483, 0.617753809, 0.433004761, 0.294464827, 0.289434639, 0.235625616, 0.251140096, 0.187049998]';
% Viral_load_data =  [2.652	2.992977778	5.202	4.636233333	6.459	5.0999875	4.153	4.0209125	2.722	2.138871429	1.698985002]'; %log
Viral_load_data = [1.699	2.652	2.678	5.202	4.248	6.459	4.991	4.153	4.188	2.722	3.014	1.699]'; %new data


time_inf_birds = [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5]'-1; % months
time_susceptible_birds = [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5]'-1; % months
time_Viral =[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5 , 6]'; % days


% params = [beta_h c_v gamma_h mu_v];
params = [4.78842852816869	548.275533316287	0.142860922436804	2.15068360150543	9.40812203003236e-06	6.99950120249231	20.3685300195583	11.6811012156835	3.54021078095144]; %bestfit sim
params = [4.78842852816869	548.275533316287	0.142860922436804	2.15068360150543	9.40812203003236e-06	6.99950120249231	20.3685300195583	11.6811012156835	3.54021078095144];
% params = [7.61744356683431	548.418079761525	0.154556825099877	5.25086908305652	8.67473759280175e-06	6.99955390382402	31.0070335550796	17.8335684786654	3.60087538971448];
% params = [7.61744356683431	548.418079761525	0.154556825099877	5.25086908305652	8.67473759280175e-06	6.99955390382402	31.0070335550796	17.8335684786654	3.60087538971448];
% params = [7.61744356683431	548.418079761525	0.154556825099877	5.25086908305652	8.67473759280175e-06	6.99955390382402	31.0070335550796	17.8335684786654	3.60087538971448];
% params = [6.31335928381058	548.281354220743	0.142877589266024	4.16878758042062	7.10512894355590e-06	9.99999515064854	19.8956082916260	13.1648109314848	3.85634194877104];
% params = [6.31335928381058	548.281354220743	0.142877589266024	4.16878758042062	7.10512894355590e-06	6.99950120249231	19.8956082916260	13.1648109314848	3.85634194877104];
% % params = [4.78842852816869	548.275533316287	0.142860922436804	2.15068360150543	9.40812203003236e-06	6.99950120249231	20.3685300195583	11.6811012156835	3.54021078095144];
params = [4.78842852816869	548.275533316287	0.142860922436804	2.15068360150543	9.40812203003236e-06	6.99950120249231	20.3685300195583	11.6811012156835	3.54021078095144];
% params = [6.31368499980693	548.281355850399	0.142877602700719	4.16852007413467	8.43066031369562e-06	6.99950087821784	19.8885886803839	13.1789571261336	3.84680718784547];
% 
% params = [6.31368499980693	548.281355850399	0.142877602700719	4.16852007413467	8.43066031369562e-06	6.99950087821784	19.8885886803839	13.1789571261336	3.84680718784547];
% params = [6.86310269933821	338.711736899123	0.142864494847676	2.87072521474700  1e-5	 6.99	19.5	10.52	3.49];
% params = [6.86310269933821	338.711736899123	0.142864494847676	2.87072521474700	1.00000000000000e-05	6.99000000000000	19.5000000000000	10.5200000000000	3.49000000000000];


params = [8.46688795601704	548.282664018811	0.143092807616721	4.97386676434899	1.28431083099837e-05	6.99999720268402	23.7225851728925	8.73266351955808	3.22291456150468];
params = [8.46688795601704	548.282664018811	0.143092807616721	4.97386676434899	1.28431083099837e-05	6.99999720268402	23.7225851728925	8.73266351955808	3.22291456150468];
params = [8.46688795602496	548.282664018811	0.143092807616721	4.97386676433347	1.28431083074797e-05	6.99999720268402	23.7225851729151	8.73266351967616	3.22291456151704];
% params = [beta_h c_v gamma_h mu_v beta d delta pi c];

% lb = [0 0 1/7 0.5 exp(-23) 0.5 0 0 0];
% ub = [100 1000 1 6 exp(2) 10 100 1e8 100];

lb = [0 0 1/7 0.5 exp(-23) 0.5 0 0 0];
ub = [100 1000 1 6 exp(2) 7 25 1e8 50]; %new bnds



% opts_con = optimset('Display','iter');
% [params, fval] =  fmincon(@err_in_data,params,[],[],[],[],lb,ub,[],opts_con);

% opts_con = optimset('Display','off', 'TolX', 1e-12, 'TolFun', 1e-12);
% [params, fval] =  fmincon(@err_in_data,params,[],[],[],[],lb,ub,[],opts_con);

options = optimset('Display', 'off', 'TOLX', 1e-12, 'TOLFun', 1e-12);
problem = createOptimProblem('fmincon', 'x0', params, 'objective', @(k)err_in_data(k), 'lb', lb, 'ub', ub, 'options', options);
ms =MultiStart('UseParallel', true,'Display','iter');    %defines a multistart problem
% Run the optimization with 20 initial guesses
[params, fval] = run(ms, problem, 100);




% Compute the model solution with the optimized parameters

[~,~,Sh_model_sol,Ih_model_sol,~] = FD_Between_Host(params);

Time = 6; % final time for between-host model
dtau = 0.01;  %within-host model discretization
kappa = 30;   %conversion for time scales
dt = 0.01/kappa;  % between host discretization
tforward = 0:dt:Time;  
infection_timeforward = 0:dtau:7;
initial_condition = [4e6 0 0 1];

within_params = params(5:9)

[t_r, y_r] = ode15s(@(t,y)Usutu_inHost(y,within_params),infection_timeforward, initial_condition);

% AIC calculation
N = length(Susceptible_birds_data)+length(Infected_birds_data)+length(Viral_load_data) ; %number of data points
KK = length(params); %number of data parameters for fitting
AIC = (N*log(fval/N) + 2*KK)
AICc=AIC + (2*KK*(KK+1))/(N - KK -1)

% Plotting 
figure;
plot(tforward, Ih_model_sol, 'b-', 'LineWidth', 2); % model solution Ih
hold on;
plot(time_inf_birds, Infected_birds_data, 'ro', 'MarkerSize', 8, 'LineWidth', 2); % observed data
xlabel('Time (months)');
ylabel('Infected Birds');

figure;
plot(tforward, Sh_model_sol, 'b-', 'LineWidth', 2); % model solution Sh
hold on;
plot(time_susceptible_birds, Susceptible_birds_data, 'ro', 'MarkerSize', 8, 'LineWidth', 2); % observed data
xlabel('Time (months)');
ylabel('Susceptible Birds');

figure;
plot(infection_timeforward, log10(y_r(:,4)), 'b-', 'LineWidth', 2); % model solution Sh
hold on;
plot(time_Viral, Viral_load_data, 'ro', 'MarkerSize', 8, 'LineWidth', 2); % observed data
xlabel('Time (days)');
ylabel('Viral Load');




function error_in_data = err_in_data(k)

Infected_birds_data =  [0.00655021834, 0.03766375546, 0.1359170306, 0.1517467249, 0.2194323144, 0.09497816594, 0.00545851528, 0.01965065502, 0.03220524017, 0]';
Susceptible_birds_data = [0.99, 0.612408962, 0.725485483, 0.617753809, 0.433004761, 0.294464827, 0.289434639, 0.235625616, 0.251140096, 0.187049998]';
Viral_load_data =  [1.699	2.652	2.678	5.202	4.248	6.459	4.991	4.153	4.188	2.722	3.014	1.699]'; %log

[~,~,Sh,Ih,~, V_model] = FD_Between_Host(k);

   
    tmeasure_inf = 1:1500:15001-1500;
    tmeasure_sus = 1:1500:15001-1500;
   
   
    Model_Ih = Ih(tmeasure_inf(:));
    Model_Sh = Sh(tmeasure_sus(:));
    Model_Viral = log10(V_model);

   % 
   error_in_data = sum((Model_Viral - Viral_load_data).^2)+...
                   sum((Model_Sh - Susceptible_birds_data).^2)+...
                   sum((Model_Ih - Infected_birds_data).^2);

    
end



% between host model finite difference soln
function [Sv,Iv,Sh,Ih,Rh,V_model] = FD_Between_Host(params)

    dtau = 0.01;
    kappa = 30;

    dt = 0.01/kappa; % dtau = kappa*dt;
    Time = 6; %months
    Infection_Age = 7; %days
 

  
    beta_h = params(1);   %bird infectivity
    cv = params(2);            %biting rate
    gamma_h = params(3);%(1/200)*30;  %bird recovery rate
     mu_v = params(4);    %mosquito birth/death rate



    within_params = params(5:9);
    
    tforward = 0:dt:Time;
    infection_timeforward = 0:dtau:Infection_Age;
    initial_condition = [4e6 0 0 1];

    time_Viral =[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5 , 6]'; % days
    sol = ode15s(@(t,y)Usutu_inHost(y,within_params),[0 7], initial_condition);
    V_model = deval(sol, time_Viral, 4)';
    V = deval(sol, infection_timeforward, 4)';%y_sol(:,4);
    V_log = log10(V);
    indc = find(V_log <=3.5);%force all V below 2 to be at 2 (prob 0 in transmission model)
    V_log(indc) = 3.5;

    beta_v = cv*(1-exp(-1.6e-3*(V_log - 3.50).^4.97));   % transmission rate per bite

        Sv = zeros(length(tforward),1); %susceptible Mos:
        Iv = zeros(length(tforward),1); %infected Mos:
        Sh = zeros(length(tforward),1); %susceptible birds
        Rh = zeros(length(tforward),1); %recoverd birds
        ih_old = zeros(length(infection_timeforward),1);
        ih_new = zeros(length(infection_timeforward),1);
        Ih = zeros(length(tforward),1); %infected birds
 
        

        Sv(1) = 0.99; %Mos: ICs
        Iv(1) = 0.01;


        Rh(1) = 0;
        ih_old(:) = 0.000291389974314388/length(ih_old);
        Ih(1) = dtau*trapz(ih_old);  %sum of infected birds
        Sh(1) = 1 - Ih(1); %bird ICs


            for n = 1:length(tforward)-1

                BvInt = dtau*trapz(beta_v.*ih_old);

                Sv(n+1) = (Sv(n) + dt*mu_v)/(1 + dt*BvInt + dt*mu_v);
                Iv(n+1) = (Iv(n) + dt*Sv(n+1)*BvInt)/(1 + dt*mu_v);
                Sh(n+1) = Sh(n)/(1 + dt*beta_h*Iv(n+1));
                ih_new(1) = (beta_h*Sh(n+1)*Iv(n+1))/kappa; %boundary condition
                ih_new(2:end) = ih_old(1:end-1)/( 1 + dt*gamma_h);
                Ih(n+1) = dtau*trapz(ih_new); %sum of infect birds of all ages
                Rh(n+1) = (1 - Sh(n+1) - Ih(n+1));

                ih_old = ih_new;
            end

    

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