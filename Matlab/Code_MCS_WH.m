clear all
close all
clc

   
numiter = 1000; 

% params = [beta d delta pi c]

true_params = [1.00432067935374e-05	6.99890940630983	19.4974816748663	10.5174125251227	3.48621792230982];
%% 
         
 
X = zeros(length(true_params),numiter); 

time_Viral =[1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5 , 6]'; % days


initial_cond = [4e6 0 0 1];


sol = ode15s(@(t,y)Usutu_inHost(y,true_params),[0 7], initial_cond);
V_model = deval(sol, time_Viral, 4)';



 Model_Viral = log10(V_model);

 noiselevel = [0.01, 0.05, 0.1, 0.2, 0.5];

 total_ARE =  zeros(length(noiselevel), length(true_params));


 total_ARE_Table = {'beta', 'd', 'delta',  'pi', 'c'};

%%
for noisei = 1:5
    
rng default
noiselev = noiselevel(noisei)

    parfor i = 1:numiter
            i

  ViralData = (noiselev*randn(length(Model_Viral),1)) + Model_Viral % Normally Distributed observations
  
            %params = [beta k delta pi c]
            lb = [exp(-23) 0.5 0 0 0];
            ub = [exp(2)    7 100 1e8 100 ];
             
            options = optimset('TolX', 1e-12);

             [k,fval] = fmincon(@(k)err_in_data(k,ViralData), true_params,[], [], [], [], lb, ub,[], options);
             fval
             X(:,i) = k';
     end
        
        arescore = zeros(1,length(true_params));

    for i = 1:length(true_params)
        arescore(i) = 100*sum(abs(true_params(i) - X(i,:))/abs(true_params(i)))/numiter;
    end
    
    total_ARE(noisei,:) = round(arescore,1);
    total_ARE_Table(noisei+1,:) = num2cell(total_ARE(noisei,:));
end


function error_in_data = err_in_data(k, ViralData) 
 

 initial_cond = [4e6 0 0 1];
time_Viral =[1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5 , 6]'; % days

sol = ode15s(@(t,y)Usutu_inHost(y,k),[0 7], initial_cond);
V_model_f = deval(sol, time_Viral, 4)';

 
 Model_Viral_f = log10(V_model_f);
 
 
 error_in_data = sum((Model_Viral_f - ViralData).^2);
                        

 end

function dy = Usutu_inHost(y,k)

dy = zeros(4,1);

%params = [beta d delta pi c]
beta = k(1);
d = k(2);
delta = k(3);
pi = k(4);
c = k(5);


T = y(1);
E = y(2);
I = y(3);
V = y(4);


dy(1) = - beta*V.*T ;
dy(2) = beta*V.*T  - d*E;
dy(3) = d*E - delta*I;
dy(4) = pi*I - c*V;
 
end