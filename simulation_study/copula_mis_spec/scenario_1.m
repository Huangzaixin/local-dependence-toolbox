%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                                            Scenario 1                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Description: The impact of copula model mis-specification on local Kendall's tau estimation.
%
% Scenario 1: the true copula is Clayton; the mis-specified copulas are Gumbel and Frank.
% NOTE: This program uses Scenario 1 as an example. Users can modify parameters to complete Scenario 2 and 3.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure:
% Step 1. Generate random data of sizes mc_sample_size = [100,500,1000] from the true copula model
% Step 2. Compute the theoretical local Kendall’s tau for each region using formula (5) and the true copula
% Step 3. Fit Clayton, Gumbel, and Frank copula models separately to each generated dataset
% Step 4. Estimate local Kendall's tau for each region using the estimated copulas in Step 3 and Formula (5)
% Step 5. For each simulation, generate bootstrap samples of size mc_sample_size with B_times replications
%         Then, compute the bootstrap standard error and 95% confidence interval.
% Step 6. Repeat steps 1-5 for mc_rep_times iterations
% Step 7. Calculate the following five performance metrics for each estimator:
%         1) Average Bias(Bias);
%         2) Root Mean Square Error(RMSE);
%         3) Empirical Standard Error(SE_e);
%         4) Average Bootstrap Standard Error(ASE_b);
%         5) Empirical 95% Coverage Probability(CP%).
% Step 8. Compare performance metrics between the correctly specified and two mis-specified copula models.
%
% Written for paper "Generalized local Kendall’s τ: a novel framework for 
%                    uncovering nonlinear local dependence" (Huang & Zhang,2026)
%
% Output: All results are stored in 'experiment_results.xlsx'
%
% Author: Zaixin Huang
% Date: finished at 2025.09.22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Basic Settings
% Copula configuration
true_copula_type = 'frank';     % true copula: Clayton; you can change this parameter in Scenarios 2 and 3
true_copula_par = 7;            % copula parameter
measure_type = 'kendall';       % dependence measure type

% Two mis-specified copula models
mis_spec_copula_type_1 = 'gumbel';     % gumbel copula; you can change this parameter in Scenarios 2 and 3
mis_spec_copula_type_2 = 'clayton';    % frank copula; you can change this parameter in Scenarios 2 and 3

% Monte Carlo simulation parameters
mc_sample_size = 100;     % Monte Carlo sample size: 100, 500, 1000
mc_rep_times = 500;       % Monte Carlo simulation times
B_times = 500;            % bootstrap times for each simulation

%% Local Kendall's Tau Configuration
% local Kendall's tau type: 'general', 'll', 'uu'
ld_type = 'general';    

% NOTE: if ld_type = 'general', you need to set the following four parameters to construct a rectangular or square region
pl = 0.45;   % lower bound of U
pu = 0.55;   % upper bound of U
ql = 0.45;   % lower bound of V
qu = 0.55;   % upper bound of V

% If ld_type = 'll', 'uu', 'ul', 'lu', setting only the following two parameters is enough, which will improve the computational efficiency.
% Of course, you can also use the first configuration method. 
% for example, ld_type = 'general' and pl = 0, pu = 0.3, ql = 0, qu = 0.2 is equivalent to ld_type = 'll' and p = 0.3, q = 0.2 
p = 0.8;    % bound of U
q = 0.8;    % bound of V

% Clear unused variables based on the selected local Kendall's tau type
if strcmpi(ld_type, 'general')
    clear p q;
else
    clear pl pu ql qu;
end

%% Pre-allocate arrays used in the simulation
% Bootstrap standard errors for each simulation
b_se_true_copula = zeros(mc_rep_times,1);          % For estimator based on the true copula model
b_se_mis_spec_copula_1 = zeros(mc_rep_times,1);    % For estimator based on the mis-specified copula model 1
b_se_mis_spec_copula_2 = zeros(mc_rep_times,1);    % For estimator based on the mis-specified copula model 2

% 95% confidence intervals for each simulation
b_conf_interval_95_true_copula = zeros(mc_rep_times,2);         % For true copula model
b_conf_interval_95_mis_spec_copula_1 = zeros(mc_rep_times,2);   % For mis-specified copula model 1
b_conf_interval_95_mis_spec_copula_2 = zeros(mc_rep_times,2);   % For mis-specified copula model 2

% Estimated local Kendall's tau for each simulation
mc_true_copula_based_local_tau = zeros(mc_rep_times,1);         % Using true copula model
mc_mis_spec_copula_1_based_local_tau = zeros(mc_rep_times,1);   % Using mis-specified copula model 1
mc_mis_spec_copula_2_based_local_tau = zeros(mc_rep_times,1);   % Using mis-specified copula model 2

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Monte Carlo Simulation                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start Monte Carlo simulation: repeat mc_rep_times times
fprintf('Starting Monte Carlo simulation with %d repetitions...\n', mc_rep_times);

for mc_index = 1:1:mc_rep_times
    fprintf('Simulation %d is in progress...left %d simulations. This program may consume much time.\n', mc_index, mc_rep_times - mc_index);
    
    % generate random data from the true copula model for current simulation
    UV = copularnd(true_copula_type,true_copula_par,mc_sample_size);
    
    % estimate local Kendall's tau using the true copula model
    temp_true_copula_par = fun_copula_estimation(true_copula_type,UV);
    if strcmpi(ld_type, 'general') 
        % store the estimated local Kendall's tau in mc_true_copula_based_local_tau(mc_index,1)
        mc_true_copula_based_local_tau(mc_index,1) = fun_copulald_general(true_copula_type,0,0,temp_true_copula_par,0,0,measure_type,pl,pu,ql,qu);
    else 
        mc_true_copula_based_local_tau(mc_index,1) = fun_copulald_type_II(true_copula_type,0,0,temp_true_copula_par,0,0,measure_type,ld_type,p,q);
    end

    % estimate local Kendall's tau using the mis-specified copula model 1
    temp_mis_copula_1_par = fun_copula_estimation(mis_spec_copula_type_1,UV);
    % store local Kendall's tau estimated by the mis-specified copula model 1
    if strcmpi(ld_type, 'general')
        mc_mis_spec_copula_1_based_local_tau(mc_index,1) = fun_copulald_general(mis_spec_copula_type_1,0,0,temp_mis_copula_1_par,0,0,measure_type,pl,pu,ql,qu);
    else 
        mc_mis_spec_copula_1_based_local_tau(mc_index,1) = fun_copulald_type_II(mis_spec_copula_type_1,0,0,temp_mis_copula_1_par,0,0,measure_type,ld_type,p,q);    
    end
        
    % estimated local Kendall's tau using the mis-specified copula model 2
    temp_mis_copula_2_par = fun_copula_estimation(mis_spec_copula_type_2,UV);
    % store local Kendall's tau estimated by the mis-specified copula model 2 
    if strcmpi(ld_type, 'general')
       mc_mis_spec_copula_2_based_local_tau(mc_index,1) = fun_copulald_general(mis_spec_copula_type_2,0,0,temp_mis_copula_2_par,0,0,measure_type,pl,pu,ql,qu); 
    else 
       mc_mis_spec_copula_2_based_local_tau(mc_index,1) = fun_copulald_type_II(mis_spec_copula_type_2,0,0,temp_mis_copula_2_par,0,0,measure_type,ld_type,p,q);
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                     Start bootstrap procedure                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate bootstrap sample (order number)
    b_matrix_index = bootstrp(B_times,@(x)x,[1:1:mc_sample_size]);
    
    % pre-allocation for bootstrap results
    b_true_copula_par = zeros(B_times,1);                    % estimated true copula parameter for each bootstrap sample
    b_true_copula_based_local_tau = zeros(B_times,1);        % estimated local Kendall's tau using the true copula model for each bootstrap sample 
    
    b_mis_spec_copula_1_par = zeros(B_times,1);              % estimated mis-specified copula 1 parameter for each bootstrap sample
    b_mis_spec_copula_1_based_local_tau = zeros(B_times,1);  % estimated local Kendall's tau using the mis-specified copula model 1 for each bootstrap sample
    
    b_mis_spec_copula_2_par = zeros(B_times,1);              % estimated mis-specified copula 2 parameter for each bootstrap sample
    b_mis_spec_copula_2_based_local_tau = zeros(B_times,1);  % estimated local Kendall's tau using the mis-specified copula model 2 for each bootstrap sample
    
    % bootstrap iterations: repeat B_times times
    for b_index = 1:1:B_times
        % get the b_index-th bootstrap sample 
        temp_bootsample_UV = UV(b_matrix_index(b_index,:),:);
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        %  calculate and store estimated local Kendall's tau for each bootstrap sample  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        %% true copula model
        b_true_copula_par(b_index,1) = fun_copula_estimation(true_copula_type,temp_bootsample_UV);
        true_copula_par_temp = b_true_copula_par(b_index,1);
        if strcmpi(ld_type, 'general')
            b_true_copula_based_local_tau(b_index,1) = fun_copulald_general(true_copula_type,0,0,true_copula_par_temp,0,0,measure_type,pl,pu,ql,qu);
        else 
            b_true_copula_based_local_tau(b_index,1) = fun_copulald_type_II(true_copula_type,0,0,true_copula_par_temp,0,0,measure_type,ld_type,p,q);
        end
        
        %% mis-specified copula model 1
        b_mis_spec_copula_1_par(b_index,1) = fun_copula_estimation(mis_spec_copula_type_1,temp_bootsample_UV);
        mis_copula_1_par_temp = b_mis_spec_copula_1_par(b_index,1);
        if strcmpi(ld_type, 'general')
            b_mis_spec_copula_1_based_local_tau(b_index,1) = fun_copulald_general(mis_spec_copula_type_1,0,0,mis_copula_1_par_temp,0,0,measure_type,pl,pu,ql,qu);
        else 
            b_mis_spec_copula_1_based_local_tau(b_index,1) = fun_copulald_type_II(mis_spec_copula_type_1,0,0,mis_copula_1_par_temp,0,0,measure_type,ld_type,p,q);
        end
        
        %% mis-specified copula model 2
        b_mis_spec_copula_2_par(b_index,1) = fun_copula_estimation(mis_spec_copula_type_2,temp_bootsample_UV);
        mis_copula_2_par_temp = b_mis_spec_copula_2_par(b_index,1);
        if strcmpi(ld_type, 'general') 
            b_mis_spec_copula_2_based_local_tau(b_index,1) = fun_copulald_general(mis_spec_copula_type_2,0,0,mis_copula_2_par_temp,0,0,measure_type,pl,pu,ql,qu);
        else 
            b_mis_spec_copula_2_based_local_tau(b_index,1) = fun_copulald_type_II(mis_spec_copula_type_2,0,0,mis_copula_2_par_temp,0,0,measure_type,ld_type,p,q);
        end
        
    end
    % end bootstrap procedure
 
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       bootstrap standard error and 95% confidence interval       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% true copula model
    mean_b_true_copula = mean(b_true_copula_based_local_tau);   % mean of the B_times bootstrap samples
    % bootstrap standard error for each simulation
    per_b_se_true_copula = sqrt(sum((b_true_copula_based_local_tau(:,1) - mean_b_true_copula).^2)./(B_times - 1));
    b_se_true_copula(mc_index,1) = per_b_se_true_copula;
    % 95% confidence interval for each simulation
    b_conf_interval_95_true_copula(mc_index,1) = mc_true_copula_based_local_tau(mc_index,1) - 1.96 * per_b_se_true_copula;
    b_conf_interval_95_true_copula(mc_index,2) = mc_true_copula_based_local_tau(mc_index,1) + 1.96 * per_b_se_true_copula;  
    
    %% mis-specified copula model 1
    mean_b_mis_spec_copula_1 = mean(b_mis_spec_copula_1_based_local_tau);
    % bootstrap standard error for each simulation
    per_b_se_mis_spec_copula_1 = sqrt(sum((b_mis_spec_copula_1_based_local_tau(:,1) - mean_b_mis_spec_copula_1).^2)./(B_times - 1));
    b_se_mis_spec_copula_1(mc_index,1) = per_b_se_mis_spec_copula_1;
    % 95% confidence interval for each simulation
    b_conf_interval_95_mis_spec_copula_1(mc_index,1) = mc_mis_spec_copula_1_based_local_tau(mc_index,1) - 1.96 * per_b_se_mis_spec_copula_1;
    b_conf_interval_95_mis_spec_copula_1(mc_index,2) = mc_mis_spec_copula_1_based_local_tau(mc_index,1) + 1.96 * per_b_se_mis_spec_copula_1;  
    
    %% mis-specified copula model 2
    mean_b_mis_spec_copula_2 = mean(b_mis_spec_copula_2_based_local_tau);
    % bootstrap standard error for each simulation
    per_b_se_mis_spec_copula_2 = sqrt(sum((b_mis_spec_copula_2_based_local_tau(:,1) - mean_b_mis_spec_copula_2).^2)./(B_times - 1));
    b_se_mis_spec_copula_2(mc_index,1) = per_b_se_mis_spec_copula_2;
    % 95% confidence interval for each simulation
    b_conf_interval_95_mis_spec_copula_2(mc_index,1) = mc_mis_spec_copula_2_based_local_tau(mc_index,1) - 1.96 * per_b_se_mis_spec_copula_2;
    b_conf_interval_95_mis_spec_copula_2(mc_index,2) = mc_mis_spec_copula_2_based_local_tau(mc_index,1) + 1.96 * per_b_se_mis_spec_copula_2;  
end
% End Monte Carlo simulation

fprintf('Monte Carlo simulation completed.\n');

%% Compute the true value of local Kendall's tau using the true copula model
if strcmpi(ld_type, 'general')
    true_value = fun_copulald_general(true_copula_type,0,0,true_copula_par,0,0,measure_type,pl,pu,ql,qu);
else
    true_value = fun_copulald_type_II(true_copula_type,0,0,true_copula_par,0,0,measure_type,ld_type,p,q);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Comparison of performance metrics                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%
% Correctly specified copula model  
% Average Bias (Bias)
mean_mc_true_copula = mean(mc_true_copula_based_local_tau);
mc_Bias_true_copula = mean_mc_true_copula - true_value;

% Root Mean Square Error (RMSE)
mc_RMSE_true_copula = sqrt(sum((mc_true_copula_based_local_tau - true_value).^2)./mc_rep_times);

% Empirical Standard Error (SE_e) 
mc_ESE_true_copula = sqrt(sum((mc_true_copula_based_local_tau - mean_mc_true_copula).^2)./(mc_rep_times - 1));

% Average Bootstrap Standard Error (ASE_b)
ave_bootstrap_SE_true_copula = mean(b_se_true_copula);

% Empirical 95% Coverage Probability (CP%)
correct_number_true_copula = 0;
for i = 1:1:mc_rep_times
    if true_value >= b_conf_interval_95_true_copula(i,1) && true_value <= b_conf_interval_95_true_copula(i,2)
        correct_number_true_copula = correct_number_true_copula + 1;
    end
end
emp_95_CP_true_copula = correct_number_true_copula/mc_rep_times;


%% Mis-specified copula model 1 
% Average Bias (Bias)
mean_mc_mis_spec_1 = mean(mc_mis_spec_copula_1_based_local_tau);
mc_Bias_mis_spec_1 = mean_mc_mis_spec_1 - true_value;

% Root Mean Square Error (RMSE)
mc_RMSE_mis_spec_1 = sqrt(sum((mc_mis_spec_copula_1_based_local_tau - true_value).^2)./mc_rep_times);

% Empirical Standard Error (SE_e)
mc_ESE_mis_spec_1 = sqrt(sum((mc_mis_spec_copula_1_based_local_tau - mean_mc_mis_spec_1).^2)./(mc_rep_times - 1));

% Average Bootstrap Standard Error (ASE_b)
ave_bootstrap_SE_mis_spec_1 = mean(b_se_mis_spec_copula_1);

% Empirical 95% Coverage Probability (CP%)
correct_number_mis_spec_1 = 0;
for i = 1:1:mc_rep_times
    if true_value >= b_conf_interval_95_mis_spec_copula_1(i,1) && true_value <= b_conf_interval_95_mis_spec_copula_1(i,2)
        correct_number_mis_spec_1 = correct_number_mis_spec_1 + 1;
    end
end
emp_95_CP_mis_spec_1 = correct_number_mis_spec_1/mc_rep_times;


%% Mis-specified copula model 2  
% Average Bias (Bias)
mean_mc_mis_spec_2 = mean(mc_mis_spec_copula_2_based_local_tau);
mc_Bias_mis_spec_2 = mean_mc_mis_spec_2 - true_value;

% Root Mean Square Error (RMSE)
mc_RMSE_mis_spec_2 = sqrt(sum((mc_mis_spec_copula_2_based_local_tau - true_value).^2)./mc_rep_times);

% Empirical Standard Error (SE_e)
mc_ESE_mis_spec_2 = sqrt(sum((mc_mis_spec_copula_2_based_local_tau - mean_mc_mis_spec_2).^2)./(mc_rep_times - 1));

% Average Bootstrap Standard Error (ASE_b)
ave_bootstrap_SE_mis_spec_2 = mean(b_se_mis_spec_copula_2);

% Empirical 95% Coverage Probability (CP%)
correct_number_mis_spec_2 = 0;
for i = 1:1:mc_rep_times
    if true_value >= b_conf_interval_95_mis_spec_copula_2(i,1) && true_value <= b_conf_interval_95_mis_spec_copula_2(i,2)
        correct_number_mis_spec_2 = correct_number_mis_spec_2 + 1;
    end
end
emp_95_CP_mis_spec_2 = correct_number_mis_spec_2/mc_rep_times;


