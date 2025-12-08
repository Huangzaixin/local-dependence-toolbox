%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                      SIMULATION STUDIES                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: Evaluate and compare the performance of U-statistic-based and copula model-based 
%              estimators under different copula model specifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure:
% Step 1. Generate random datasets of sizes mc_sample_size = [100,200,500,1000,2000] from the specified copula model.
% Step 2. Estimate the specified copula model using each generated dataset.
% Step 3. Compute local Kendall’s tau in the selected region using both the U-statistic estimator and copula-based estimator, respectively.
% Step 4. For each simulation, generate bootstrap samples of size mc_sample_size with B_times replications.
%         Then, compute the bootstrap standard error and 95% confidence interval for each estimator.
% Step 5. Repeat steps 1-4 for mc_rep_times iterations.
% Step 6. Calculate the following performance metrics for each estimator:
%         1) Average Bias(Bias); 
%         2) Empirical Standard Error(SE_e);
%         3) Average Bootstrap Standard Error(ASE_b); 
%         4) Empirical 95% Coverage Probability(CP%)
% Step 7. Compare performance metrics between two estimators.
%
% Written for paper "Generalized local Kendall’s τ: a novel framework for 
%                     uncovering nonlinear local dependence" (Huang & Zhang,2026)
%
% Output: All results are stored in 'simulation_results.xlsx'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic Settings
% Copula configuration
copula_type = 'clayton';    % copula type: clayton, gumbel, frank 
copula_par = 1;             % copula parameter
measure_type = 'kendall';   % dependence measure type

% Monte Carlo simulation parameters
mc_sample_size = 100;       % sample size: 100, 200, 500, 1000, and 2000
mc_rep_times = 1000;        % number of Monte Carlo repetitions
B_times = 1000;             % bootstrap replications in each simulation

% Marginal distribution parameters
x_t_dof = 5;       % degree of freedom for X (t-distribution)
y_t_dof = 7;       % degree of freedom for Y (t-distribution)

%% Local Kendall's Tau Configuration
% local Kendall's tau type
ld_type = 'general';    % Options: 'general', 'll', 'uu'
                        % 'general': generalized local Kendall's tau
                        % 'll': lower-lower local Kendall's tau
                        % 'uu': upper-upper local Kendall's tau

% NOTE: if ld_type = 'general', define a rectangular region by setting the following four boundary parameters
Xl = -1;       % lower bound of X
Xu = 1;        % upper bound of X
Yl = -1;       % lower bound of Y
Yu = 1;        % upper bound of Y

% Convert bounds to quantiles
pl = tcdf(Xl,x_t_dof);     % quantile of lower bound
pu = tcdf(Xu,x_t_dof);     % quantile of upper bound
ql = tcdf(Yl,y_t_dof);     % quantile of lower bound
qu = tcdf(Yu,y_t_dof);     % quantile of upper bound

% If ld_type = 'll' or 'uu', set only the following two parameters is enough. 
% Of course, you can also use the first configuration method.
X = -2;     % bound of X
Y = -2;     % bound of Y

% Convert bounds to quantiles
p = tcdf(X,x_t_dof);      % quantile of X
q = tcdf(Y,y_t_dof);      % quantile of Y

% Clear unused variables based on the selected local Kendall's tau type
if strcmpi(ld_type, 'general')
    clear X Y p q;
else
    clear Xl Xu Yl Yu pl pu ql qu;
end

%% Pre-allocate arrays used in the simulation
% Bootstrap standard errors for each simulation
b_se_nonparametric = zeros(mc_rep_times,1);    % U-statistic-based estimator
b_se_copula_based = zeros(mc_rep_times,1);     % copula model-based estimator

% 95% confidence intervals for each simulation
b_conf_interval_95_nonparametric = zeros(mc_rep_times,2);   % U-statistic-based estimator
b_conf_interval_95_copula_based = zeros(mc_rep_times,2);    % copula model-based estimator

% Estimated local Kendall's tau for each simulation
mc_nonparametric_local_tau = zeros(mc_rep_times,1);   % estimated by the U-statistic-based estimator
mc_copula_based_local_tau = zeros(mc_rep_times,1);    % estimated by the copula model-based estimator

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Monte Carlo Simulation                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start Monte Carlo Simulation: repeat mc_rep_times times
fprintf('Starting Monte Carlo simulation with %d repetitions...\n', mc_rep_times);

for mc_index = 1:1:mc_rep_times
    fprintf('Simulation %d is in progress...left %d simulations. This program may consume much time.\n', mc_index, mc_rep_times - mc_index);
    
    % generate random data from the specified copula model
    UV = copularnd(copula_type,copula_par,mc_sample_size);
    
    % convert U, V into X，Y，X~t(5), Y~t(7)
    XY(:,1) = tinv(UV(:,1),x_t_dof);
    XY(:,2) = tinv(UV(:,2),y_t_dof);
    
    % select data points in the target region
    switch ld_type
        case 'general'
            selected_dataset = find(Xl <= XY(:,1) & XY(:,1) <= Xu & Yl <= XY(:,2) & XY(:,2) <= Yu);
        case 'll'
            selected_dataset = find(XY(:,1) <= X & XY(:,2) <= Y);
        case 'uu'
            selected_dataset = find(XY(:,1) >= X & XY(:,2) >= Y);
        otherwise
    end
    
    % if number of data points in the target region is less than 10, re-generate random data
    while (length(selected_dataset) < 10)
        % re-simulate
        UV = copularnd(copula_type,copula_par,mc_sample_size);
        % convert U, V into X，Y, X~t(5), Y~t(7)
        XY(:,1) = tinv(UV(:,1),x_t_dof);
        XY(:,2) = tinv(UV(:,2),y_t_dof);
        
        % dataset in the selected region, again
        switch ld_type
            case 'general'
                selected_dataset = find(Xl <= XY(:,1) & XY(:,1) <= Xu & Yl <= XY(:,2) & XY(:,2) <= Yu);
            case 'll'
                selected_dataset = find(XY(:,1) <= X & XY(:,2) <= Y);
            case 'uu'
                selected_dataset = find(XY(:,1) >= X & XY(:,2) >= Y);
            otherwise
        end
    end
    
    % Estimate local Kendall's tau using two estimators
    % U-statistic-based estimator
    if strcmpi(ld_type, 'general')  
        mc_nonparametric_local_tau(mc_index,1) = fun_u_statistic_based_ld_general(XY,Xl,Xu,Yl,Yu);
    else
        mc_nonparametric_local_tau(mc_index,1) = fun_u_statistic_based_ld_type_II(XY,ld_type,X,Y);
    end    
    
    % copula model-based estimator
    temp_copula_par = fun_copula_estimation(copula_type,UV);
    if strcmpi(ld_type, 'general')
        mc_copula_based_local_tau(mc_index,1) = fun_copulald_general(copula_type,0,0,temp_copula_par,0,0,measure_type,pl,pu,ql,qu);
    else
        mc_copula_based_local_tau(mc_index,1) = fun_copulald_type_II(copula_type,0,0,temp_copula_par,0,0,measure_type,ld_type,p,q);
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                     Start bootstrap procedure                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate bootstrap sample (order number)
    b_matrix_index = bootstrp(B_times,@(x)x,[1:1:mc_sample_size]);
    
    % pre-allocation for bootstrap results
    b_nonparametric_local_tau = zeros(B_times,1);
    b_copula_par = zeros(B_times,1);
    b_copula_based_local_tau = zeros(B_times,1);
    
    for b_index = 1:1:B_times
        % generate bootstrap sample
        temp_bootsample_UV = UV(b_matrix_index(b_index,:),:);
        XY_temp(:,1) = tinv(temp_bootsample_UV(:,1),x_t_dof);
        XY_temp(:,2) = tinv(temp_bootsample_UV(:,2),y_t_dof);
        
        switch ld_type
           case 'general'
               b_selected_dataset = find(Xl <= XY_temp(:,1) & XY_temp(:,1) <= Xu & Yl <= XY_temp(:,2) & XY_temp(:,2) <= Yu);
           case 'll'
               b_selected_dataset = find(XY_temp(:,1) <= X & XY_temp(:,2) <= Y);
           case 'uu'
               b_selected_dataset = find(XY_temp(:,1) >= X & XY_temp(:,2) >= Y);
           otherwise 
        end
        
        % if number of data points in the local region is less than 10, re-boostrap
        while (length(b_selected_dataset) < 10)
           % re-bootstrap
           temp_bootsample_UV = UV(bootstrp(1,@(x)x,[1:1:mc_sample_size]),:);
           XY_temp(:,1) = tinv(temp_bootsample_UV(:,1),x_t_dof);
           XY_temp(:,2) = tinv(temp_bootsample_UV(:,2),y_t_dof);
           
           switch ld_type
              case 'general'
                  b_selected_dataset = find(Xl <= XY_temp(:,1) & XY_temp(:,1) <= Xu & Yl <= XY_temp(:,2) & XY_temp(:,2) <= Yu);
              case 'll'
                  b_selected_dataset = find(XY_temp(:,1) <= X & XY_temp(:,2) <= Y);
              case 'uu'
                  b_selected_dataset = find(XY_temp(:,1) >= X & XY_temp(:,2) >= Y);
              otherwise 
           end
        end
        
        % U-statistic-based estimator
        if strcmpi(ld_type, 'general')
            b_nonparametric_local_tau(b_index,1) = u_statistic_based_ld_general(XY_temp,Xl,Xu,Yl,Yu);
        else
            b_nonparametric_local_tau(b_index,1) = fun_u_statistic_based_ld_type_II(XY_temp,ld_type,X,Y);
        end
        
        % copula model-based estimator
        b_copula_par(b_index,1) = fun_copula_estimation(copula_type,temp_bootsample_UV);
        copula_par_temp = b_copula_par(b_index,1);
        if strcmpi(ld_type, 'general')
            b_copula_based_local_tau(b_index,1) = fun_copulald_general(copula_type,0,0,copula_par_temp,0,0,measure_type,pl,pu,ql,qu);
        else
            b_copula_based_local_tau(b_index,1) = fun_copulald_type_II(copula_type,0,0,copula_par_temp,0,0,measure_type,ld_type,p,q);
        end
    end
    % end bootstrap procedure
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       bootstrap standard error and 95% confidence interval       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % U-statistic-based estimator
    mean_b_nonparametric = mean(b_nonparametric_local_tau);
    % standard error
    se_nonparametric = sqrt(sum((b_nonparametric_local_tau(:,1) - mean_b_nonparametric).^2)./(B_times - 1));
    b_se_nonparametric(mc_index,1) = se_nonparametric;
    % 95% confidence interval
    b_conf_interval_95_nonparametric(mc_index,1) = mc_nonparametric_local_tau(mc_index,1) - 1.96 * se_nonparametric;
    b_conf_interval_95_nonparametric(mc_index,2) = mc_nonparametric_local_tau(mc_index,1) + 1.96 * se_nonparametric;
    
    % copula model-based estimator
    mean_b_copula_based = mean(b_copula_based_local_tau);
    % standard error
    se_copula_based = sqrt(sum((b_copula_based_local_tau(:,1) - mean_b_copula_based).^2)./(B_times - 1));
    b_se_copula_based(mc_index,1) = se_copula_based;
    % 95% confidence interval
    b_conf_interval_95_copula_based(mc_index,1) = mc_copula_based_local_tau(mc_index,1) - 1.96 * se_copula_based;
    b_conf_interval_95_copula_based(mc_index,2) = mc_copula_based_local_tau(mc_index,1) + 1.96 * se_copula_based;  
end
% End Monte Carlo Simulation

fprintf('Monte Carlo simulation completed.\n');

%% Compute the true value of local Kendall's tau using the true copula model
if strcmpi(ld_type, 'general')
    true_value = fun_copulald_general(copula_type,0,0,copula_par,0,0,measure_type,pl,pu,ql,qu);
else
    true_value = fun_copulald_type_II(copula_type,0,0,copula_par,0,0,measure_type,ld_type,p,q);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Comparison of performance metrics                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%
% Average Bias (Bias)
% U-statistic estimator
mean_mc_nonparametric = mean(mc_nonparametric_local_tau);
mc_bias_nonparametric = mean_mc_nonparametric - true_value;
% copula model-based estimator
mean_mc_copula_based_estimator = mean(mc_copula_based_local_tau);
mc_bias_copula_based_estimator = mean_mc_copula_based_estimator - true_value;

% Empirical Standard Error (SE_e) 
% U-statistic estimator
mc_se_nonparametric_estimator = sqrt(sum((mc_nonparametric_local_tau(:,1) - mean_mc_nonparametric).^2)./(mc_rep_times - 1));
% copula model-based estimator
mc_se_copula_based_estimator = sqrt(sum((mc_copula_based_local_tau(:,1) - mean_mc_copula_based_estimator).^2)./(mc_rep_times - 1));

% Average Bootstrap Standard Error (ASE_b)
% U-statistic estimator
ave_bootstrap_se_nonparametric_estimator = mean(b_se_nonparametric);
% copula model-based estimator
ave_bootstrap_se_copula_based_estimator = mean(b_se_copula_based);

% Empirical 95% Coverage Probability (CP%)
% U-statistic estimator
correct_number = 0;
for i = 1:1:mc_rep_times
    % check whether true_value is in the 95% confidence interval
    if true_value >= b_conf_interval_95_nonparametric(i,1) && true_value <= b_conf_interval_95_nonparametric(i,2)
        correct_number = correct_number + 1;
    end
end
emp_95_cp_nonparametric = correct_number/mc_rep_times;

% copula model-based estimator
correct_number = 0;
for i = 1:1:mc_rep_times
    if true_value >= b_conf_interval_95_copula_based(i,1) && true_value <= b_conf_interval_95_copula_based(i,2)
        correct_number = correct_number + 1;
    end
end
emp_95_cp_copula_based = correct_number/mc_rep_times;


