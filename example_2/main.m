%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Example 2: Estimate Global and Local Kendall's Tau between COVID-19 Cumulative Cases and Number of Deaths  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic Settings
options = optimset('Display','off','TolCon',10^-8,'TolFun',10^-8,'TolX',10^-8);
measure_type = 'Kendall';

%% Read Data
u = csvread('data_examples/example_2/data/u_cases.csv',1,0);    % cumulative cases
v = csvread('data_examples/example_2/data/v_deaths.csv',1,0);   % number of deaths
UV = [u v];

%% Parameter(s) of the copula model
copula_type = 'gumbel';
weight1 = 0;                        % the first copula's weight of the mixture copula model
weight2 = 0;                        % the second copula's weight of the mixture copula model
copula_parameter1 = par_gumbel(1);  % for mixture copula model, parameter of the first copula function; for sjc copula, the upper tail dependence coefficient   
copula_parameter2 = 0;              % for mixture copula model, parameter of the second copula function; for sjc copula, the lower tail dependence coefficient 
copula_parameter3 = 0;              % for mixture copula model, parameter of the third copula function

%% Four quantile bounds of the generalized local Kendall's tau
% the following quantiles are calculated by the inverse of marginal distributions of cumulative cases and number of deaths
pl = 0.9610925;    % lower bound of u
pu = 1;            % upper bound of u
ql = 0.9637547;    % lower bound of v
qu = 1;            % upper bound of v

% ld_type = 'll';   % "uu": upper-upper local Kendall's tau; "ll": lower-lower local Kendall's tau; "ul": upper-lower local Kendall's tau; "lu": lower-upper local Kendall's tau
% u_quantile = 0.1;
% v_quantile = 0.1;

%% Initiate values for Gumbel copula
lower = [1.001];
upper = [+Inf];
par0 = [3];

%% Use two estimators to estimate local Kendall's tau respectively
% copula model-based estimator
% copula_based_local_tau = fun_copulald_type_II(copula_type,weight1,weight2,copula_parameter1,copula_parameter2,copula_parameter3,measure_type,ld_type,0.05,0.05);
copula_based_local_tau = fun_copulald_general(copula_type,weight1,weight2,copula_parameter1,copula_parameter2,copula_parameter3,measure_type,pl,pu,ql,qu);

% U-statistic-based estimator
% nonparametric_local_tau = fun_u_statistic_based_ld_type_II(UV,ld_type,u_quantile,v_quantile);
nonparametric_local_tau = fun_u_statistic_based_ld_general(UV,pl,pu,ql,qu);

%% Calculate bootstrap standard error and 95% confidence interval
% start bootstrap: 1000 times
B_times = 1000;
b_nonparametric_local_tau = zeros(B_times,1);
b_copula_par = zeros(B_times,5);
b_copula_based_local_tau = zeros(B_times,1);
b_matrix_index = bootstrp(B_times,@(x)x,[1:1:length(u)]);
  
%% NOTE: the following while loop may consume much time. It depends on your computer's processing speed.
b_success_number = 0;
b_index = 0;
while b_success_number < B_times
    b_index = b_index + 1;
    
    % generate bootstrap sample
    temp_bootsample_UV = UV(b_matrix_index(b_index,:),:);
    
    %% U-statistic-based estimator
    % b_nonparametric_local_tau(b_index,1) = fun_u_statistic_based_ld_type_II(temp_bootsample_UV,ld_type,u_quantile,v_quantile);
    b_nonparametric_local_tau(b_index,1) = fun_u_statistic_based_ld_general(temp_bootsample_UV,pl,pu,ql,qu);
        
    %% copula model-based estimator
    % NOTE: if the 'copula_type' is other copula models, do not forget to change the objective function
    [parameters LL] = fmincon('gumbelCL',par0,[],[],[],[],lower,upper,[],options,[temp_bootsample_UV(:,1),temp_bootsample_UV(:,2)]);
   
    % NOTE: if you use mixture copula models, do not forget to revise the following parameters
    temp_weight1 = 0;                        % for two-component mixture model, temp_weight1 should not be 0
    temp_weight2 = 0;                        % for three-component mixture model, temp_weight2 should not be 0
    temp_copulaparameter1 = parameters(1);
    temp_copulaparameter2 = 0;               % for two-component mixture model, temp_copulaparameter2 should not be 0
    temp_copulaparameter3 = 0;               % for three-component mixture model, temp_copulaparameter3 should not be 0
 
    % b_copula_based_local_tau(b_index,1) = fun_copulald_type_II(copula_type,temp_weight1,temp_weight2,temp_copulaparameter1,temp_copulaparameter2,temp_copulaparameter3,measure_type,ld_type,u_quantile,v_quantile);
    b_temp = fun_copulald_general(copula_type,temp_weight1,temp_weight2,temp_copulaparameter1,temp_copulaparameter2,temp_copulaparameter3,measure_type,pl,pu,ql,qu);
    
    % if b_temp is NaN，fail
    if isnan(b_temp)
        continue;
    else  % success    
        b_success_number = b_success_number + 1;
        b_copula_based_local_tau(b_success_number,1) = b_temp;
    end
    
    left_number = 1000 - b_success_number;
    disp("success number:" + b_success_number + "; left:" + left_number);
end
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Bootstrap standard error and 95% confidence interval          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% for copula model-based estimator
mean_b_copula_based = mean(b_copula_based_local_tau);
% standard error(SE)
se_copula_based = sqrt(sum((b_copula_based_local_tau(:,1) - mean_b_copula_based).^2)./(B_times - 1));
% 95% confidence interval
b_conf_interval_95_copula_based_lower = copula_based_local_tau - 1.96 * se_copula_based;
b_conf_interval_95_copula_based_upper = copula_based_local_tau + 1.96 * se_copula_based;  

% for U-statistic-based estimator
mean_b_nonparametric = mean(b_nonparametric_local_tau);
% standard error(SE)
se_nonparametric = sqrt(sum((b_nonparametric_local_tau(:,1) - mean_b_nonparametric).^2)./(B_times - 1));
% 95% confidence interval
b_conf_interval_95_nonparametric_lower = nonparametric_local_tau - 1.96 * se_nonparametric;
b_conf_interval_95_nonparametric_upper = nonparametric_local_tau + 1.96 * se_nonparametric;


