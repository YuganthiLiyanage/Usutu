clear all
close all
clc


% Data: infected birds & susceptible birds & Viral load

Infected_birds_data = [0.00655021834, 0.03766375546, 0.1359170306, 0.1517467249, 0.2194323144, 0.09497816594, 0.00545851528, 0.01965065502, 0.03220524017, 0]';
Susceptible_birds_data = [0.99, 0.612408962, 0.725485483, 0.617753809, 0.433004761, 0.294464827, 0.289434639, 0.235625616, 0.251140096, 0.187049998]';
%Viral_load_data =  [2.546	3.363	4.858	5.155	5.918	5.548	3.722	4.277	2.994	2.055	1.699]'; %log

time_inf_birds = [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5]'-1; % months
time_susceptible_birds = [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5]'-1; % months

%time_Viral =[1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5 , 6]'; % days


% params = [beta_h c_v gamma_h mu_v];


params = [4.21468515351522	429.493890865587	0.142941678934300	0.500129336298220];	%best fitted parameters
% params = [4.78842852816869	548.275533316287	0.142860922436804	2.15068360150543];
% params = [4.78842852816869	548.275533316287	0.142860922436804	2.15068360150543];
% params = [4.80711866139034	428.227996244537	0.142860922748490	2.16915607530849];
params = [4.79735745785129	429.498930465263	0.142861115773887	2.17072396126139];
params = [4.79735745785129	429.498930465263	0.142861115773887	2.17072396126139];
params = [6.89249999368509	318.313458526789	0.142863024178789	2.89043666443029];
params = [6.89249999368509	318.313458526789	0.142863024178789	2.89043666443029];  %10
params = [9.37220799823523	215.070600440137	0.142867623104870	2.85378365824097]; %7
params = [6.89265737464509	318.302253475800	0.142864494921824	2.89040378789559];
params = [6.89265737464509	318.302253475800	0.142864494921824	2.89040378789559];
params = [6.86013289731686	427.942992925823	0.142864495036372	2.86857681457153];
params = [6.86013289731686	427.942992925823	0.142864495036372	2.86857681457153];
params = [6.86310269933821	338.711736899123	0.142864494847676	2.87072521474700];
params = [6.86310269933821	338.711736899123	0.142864494847676	2.87072521474700];
params = [6.86310269933821	338.711736899123	0.142864494847676	2.87072521474700];

% params = [beta_h c_v gamma_h mu_v];

lb = [0 0 1/7 0.5 ];
ub = [100 1000 1 6];




opts_con = optimset('Display','iter');
[params, fval] =  fmincon(@err_in_data,params,[],[],[],[],lb,ub,[],opts_con);

% options = optimset('Display', 'off', 'TOLX', 1e-12, 'TOLFun', 1e-12);
% problem = createOptimProblem('fmincon', 'x0', params, 'objective', @(k)err_in_data(k), 'lb', lb, 'ub', ub, 'options', options);
% ms =MultiStart('UseParallel', true,'Display','iter');    %defines a multistart problem
% % Run the optimization with 20 initial guesses
% [params, fval] = run(ms, problem, 1000);




% Compute the model solution with the optimized parameters

[~,~,Sh_model_sol,Ih_model_sol,~] = FD_Between_Host(params);

Time = 6; % final time for between-host model
dtau = 0.01;  %within-host model discretization
kappa = 30;   %conversion for time scales
dt = 0.01/kappa;  % between host discretization
tforward = 0:dt:Time;  
infection_timeforward = 0:dtau:7;
initial_condition = [4e6 0 0 1];

% within_params = [9.0e-6	 3.69	8.38	15.55	9.03];  %new fitting values
within_params = [1e-5	 6.99	19.5	10.52	3.49];  %new fitting values

[t_r, y_r] = ode15s(@(t,y)Usutu_inHost(y,within_params),infection_timeforward, initial_condition);

% AIC calculation
% N = length(Susceptible_birds_data)+length(Infected_birds_data) ; %number of data points
% KK = length(params); %number of data parameters for fitting
% AIC = (N*log(fval/N) + 2*KK)
% AICc=AIC + (2*KK*(KK+1))/(N - KK -1)

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




function error_in_data = err_in_data(k)

Infected_birds_data =  [0.00655021834, 0.03766375546, 0.1359170306, 0.1517467249, 0.2194323144, 0.09497816594, 0.00545851528, 0.01965065502, 0.03220524017, 0]';
Susceptible_birds_data = [0.99, 0.612408962, 0.725485483, 0.617753809, 0.433004761, 0.294464827, 0.289434639, 0.235625616, 0.251140096, 0.187049998]';

% Viral_load_data =  [2.546	3.363	4.858	5.155	5.918	5.548	3.722	4.277	2.994	2.055	1.699]'; %log

[~,~,Sh,Ih,~, V_model] = FD_Between_Host(k);

   
    tmeasure_inf = 1:1500:15001-1500;
    tmeasure_sus = 1:1500:15001-1500;
   
   
    Model_Ih = Ih(tmeasure_inf(:));
    Model_Sh = Sh(tmeasure_sus(:));
    % Model_Viral = log10(V_model);


   error_in_data = sum((Model_Sh - Susceptible_birds_data).^2)+...
                   10*sum((Model_Ih - Infected_birds_data).^2);

     %error_in_data = sum((Model_Ih - Infected_birds_data).^2);
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



    within_params = [1e-5	 6.99	19.5	10.52	3.49];  %new within fitting values
    
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