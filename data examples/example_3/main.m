%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Example 3: estimate global and local Kendall's tau between X and Y   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% basic settings
options = optimset('Display','off','TolCon',10^-8,'TolFun',10^-8,'TolX',10^-8);
measuretype = 'Kendall';

%% read data
u = csvread('data examples/example_3/data/u.csv',1,0);
v = csvread('data examples/example_3/data/v.csv',1,0);
UV = [u v];

%% parameters of the copula model
copula_type = 'mixfg';
weight1 = pars_fg(1);              % the weight of the first copula function in the mixture copula model
weight2 = 0;                       % the weight of the second copula function in the mixture copula model
copulaparameter1 = pars_fg(2);     % for mixture copula model, parameter of the first copula function; for sjc copula, the upper tail dependence coefficient  
copulaparameter2 = pars_fg(3);     % for mixture copula model, parameter of the second copula function; for sjc copula, the lower tail dependence coefficient 
copulaparameter3 = 0;              % for mixture copula model, parameter of the third copula function

%% four quantile parameters of the generalized local Kendall's tau
ul_quantile = 0;      % lower bound of u
uu_quantile = 1;      % upper bound of u
vl_quantile = 0;      % lower bound of v
vu_quantile = 1;      % upper bound of v

%% initiate values for MixFG copula
lower = [0.001 1.001 0.001];
upper = [0.999 +Inf +Inf];
par0 = [0.5 2 2];

%% use two estimators to estimate local Kendall's tau
% copula model-based estimator
% copula_based_local_tau = fun_copulald_type_II(copula_type,weight1,weight2,copulaparameter1,copulaparameter2,copulaparameter3,'kendall',ld_type,0.05,0.05);
copula_based_local_tau = fun_copulald_general(copula_type,weight1,weight2,copulaparameter1,copulaparameter2,copulaparameter3,measuretype,ul_quantile,uu_quantile,vl_quantile,vu_quantile);

% U-statistic-based estimator
% nonparametric_local_tau = fun_u_statistic_based_ld_type_II(UV,ld_type,u_quantile,v_quantile);
nonparametric_local_tau = fun_u_statistic_based_ld_general(UV,ul_quantile,uu_quantile,vl_quantile,vu_quantile);

%% calculate bootstrap standard error and 95% confidence interval
% start bootstrap: 1000 times
B_times = 1000;
b_nonparametric_local_tau = zeros(B_times,1);
b_copula_par = zeros(B_times,5);
b_copula_based_local_tau = zeros(B_times,1);
b_matrix_index = bootstrp(2000,@(x)x,[1:1:length(u)]);

%% NOTE: the following while loop may consume much time. It depends on your computer's processing speed.
b_success_number = 0;
b_index = 0;
while b_success_number < B_times
    b_index = b_index + 1;
    
    % generate bootstrap sample
    temp_bootsample_UV = UV(b_matrix_index(b_index,:),:);
    
    %% U-statistic-based estimator
    % b_nonparametric_local_tau(b_index,1) = fun_u_statistic_based_ld_type_II(temp_bootsample_UV,ld_type,u_quantile,v_quantile);
    b_nonparametric_local_tau(b_index,1) = fun_u_statistic_based_ld_general(temp_bootsample_UV,ul_quantile,uu_quantile,vl_quantile,vu_quantile);
        
    %% copula model-based estimator
    % NOTE: if the 'copula_type' is other copula models, do not forget to change the objective function
    [parameters LL] = fmincon('mixFGCL',par0,[],[],[],[],lower,upper,[],options,[temp_bootsample_UV(:,1),temp_bootsample_UV(:,2)]);
    
    temp_weight1 = parameters(1);             % for two-component mixture model, temp_weight1 should not be 0
    temp_weight2 = 0;                         % for three-component mixture model, temp_weight2 should not be 0
    temp_copulaparameter1 = parameters(2);
    temp_copulaparameter2 = parameters(3);    % for two-component mixture model, temp_copulaparameter2 should not be 0
    temp_copulaparameter3 = 0;                % for three-component mixture model, temp_copulaparameter3 should not be 0
 
    % b_copula_based_local_tau(b_index,1) = fun_copulald_type_II(copula_type,temp_weight1,temp_weight2,temp_copulaparameter1,temp_copulaparameter2,temp_copulaparameter3,measuretype,ld_type,u_quantile,v_quantile);
    b_temp = fun_copulald_general(copula_type,temp_weight1,temp_weight2,temp_copulaparameter1,temp_copulaparameter2,temp_copulaparameter3,measuretype,ul_quantile,uu_quantile,vl_quantile,vu_quantile);

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

%% bootstrap standard error and 95% confidence interval
% for copula model-based estimator
mean_b_copula_based = mean(b_copula_based_local_tau);
% standard error
se_copula_based = sqrt(sum((b_copula_based_local_tau(:,1) - mean_b_copula_based).^2)./(B_times - 1));
% 95% confidence interval
b_conf_interval_95_copula_based_lower = copula_based_local_tau - 1.96 * se_copula_based;
b_conf_interval_95_copula_based_upper = copula_based_local_tau + 1.96 * se_copula_based;  

% for U-statistic-based estimator
mean_b_nonparametric = mean(b_nonparametric_local_tau);
% standard error
se_nonparametric = sqrt(sum((b_nonparametric_local_tau(:,1) - mean_b_nonparametric).^2)./(B_times - 1));
% 95% confidence interval
b_conf_interval_95_nonparametric_lower = nonparametric_local_tau - 1.96 * se_nonparametric;
b_conf_interval_95_nonparametric_upper = nonparametric_local_tau + 1.96 * se_nonparametric;

