%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       SIMULATION STUDIES                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  All results are stored in the file "simulation results.xlsx"  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialization
copula_type = 'clayton';    % copula type: clayton, gumbel, frank 
copula_par = 1;             % copula parameter
measure_type = 'kendall';
mc_sample_size = 100;       % Monte Carlo sample size: 100, 200, 500, 1000, 2000

%%
mc_rep_times = 1000;         % Monte Carlo simulation times
B_times = 1000;              % bootstrap times for each simulation sample
b_se_nonparametric = zeros(B_times,1);       % bootstrap standard error, U-statistic-based estimator
b_se_copula_based = zeros(B_times,1);        % bootstrap standard error, copula model-based estimator
b_conf_interval_95_nonparametric = zeros(mc_rep_times,2);        % 95% confidence interval, U-statistic-based estimator
b_conf_interval_95_copula_based = zeros(mc_rep_times,2);         % 95% confidence interval, copula model-based estimator
mc_nonparametric_local_tau = zeros(mc_rep_times,1);       % U-statistic-based estimator
mc_copula_based_local_tau = zeros(mc_rep_times,1);        % copula model-based estimator

%% Type II local Kendall's tau
ld_type = 'll';    % 'll': lower-lower local Kendall's tau; 'uu': upper-upper local Kendall's tau
X = -2;
Y = -2;
x_t_dof = 5;       % degree of freedom of X
y_t_dof = 7;       % degree of freedom of Y
x_quantile = tcdf(X,x_t_dof);      % quantile of X
y_quantile = tcdf(Y,y_t_dof);      % quantile of Y

%% For Type I local Kendall's tau, please use the following statements
% ld_type = 'general';
% Xl = -1;         % lower bound of X
% Xu = 1;          % upper bound of X
% Yl = -1;         % lower bound of Y
% Yu = 1;          % upper bound of Y
% x_t_dof = 5;     % degree of freedom of X
% y_t_dof = 7;     % degree of freedom of Y
% xl_quantile = tcdf(Xl,x_t_dof);     % quantile of Xl
% xu_quantile = tcdf(Xu,x_t_dof);     % quantile of Xu
% yl_quantile = tcdf(Yl,y_t_dof);     % quantile of Yl
% yu_quantile = tcdf(Yu,y_t_dof);     % quantile of Yu

%% Monte Carlo simulation: repeat 1000 times
for mc_index = 1:1:mc_rep_times
    fprintf('The %dth simulation is in progress...left %d simulations. This program may consume much time.\n', mc_index, mc_rep_times - mc_index);
    
    % generate data from the copula model
    UV = copularnd(copula_type,copula_par,mc_sample_size);
    
    % convert U, V into X，Y，X~t(5), Y~t(7)
    XY(:,1) = tinv(UV(:,1),x_t_dof);
    XY(:,2) = tinv(UV(:,2),y_t_dof);
    
    % dataset in the selected region
    switch ld_type
        case 'general'
            selected_dataset = find(Xl <= XY(:,1) & XY(:,1) <= Xu & Yl <= XY(:,2) & XY(:,2) <= Yu);
        case 'll'
            selected_dataset = find(XY(:,1) <= X & XY(:,2) <= Y);
        case 'uu'
            selected_dataset = find(XY(:,1) >= X & XY(:,2) >= Y);
        otherwise
    end
    
    % if number of data points in the selected region is less than 10, re-generate random data
    while (length(selected_dataset) < 10)
        % re-simulate
        UV = copularnd(copula_type,copula_par,mc_sample_size);
        % convert U, V into X，Y, X~t(5), Y~t(7)
        XY(:,1) = tinv(UV(:,1),x_t_dof);
        XY(:,2) = tinv(UV(:,2),y_t_dof);
        
        % dataset in the selected region again
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
    
    % U-statistic-based estimator
    mc_nonparametric_local_tau(mc_index,1) = fun_u_statistic_based_ld_type_II(XY,ld_type,X,Y);
    % mc_nonparametric_local_tau(mc_index,1) = fun_u_statistic_based_ld_general(XY,Xl,Xu,Yl,Yu);
    
    % copula model-based estimator
    temp_copula_par = fun_copula_estimation(copula_type,UV);
    mc_copula_based_local_tau(mc_index,1) = fun_copulald_type_II(copula_type,0,0,temp_copula_par,0,0,measure_type,ld_type,x_quantile,y_quantile);
    % mc_copula_based_local_tau(mc_index,1) = fun_copulald_general(copula_type,0,0,temp_copula_par,0,0,measure_type,xl_quantile,xu_quantile,yl_quantile,yu_quantile);
    
    %% start bootstrap: 1000 times
    b_nonparametric_local_tau = zeros(B_times,1);
    b_copula_par = zeros(B_times,1);
    b_copula_based_local_tau = zeros(B_times,1);
    b_matrix_index = bootstrp(B_times,@(x)x,[1:1:mc_sample_size]);
    
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
        b_nonparametric_local_tau(b_index,1) = fun_u_statistic_based_ld_type_II(XY_temp,ld_type,X,Y);
        % b_nonparametric_local_tau(b_index,1) = u_statistic_based_ld_general(XY_temp,Xl,Xu,Yl,Yu);
        
        % copula model-based estimator
        b_copula_par(b_index,1) = fun_copula_estimation(copula_type,temp_bootsample_UV);
        copula_par_temp = b_copula_par(b_index,1);
        b_copula_based_local_tau(b_index,1) = fun_copulald_type_II(copula_type,0,0,copula_par_temp,0,0,measure_type,ld_type,x_quantile,y_quantile);
        % b_copula_based_local_tau(b_index,1) = fun_copulald_general(copula_type,0,0,copula_par_temp,0,0,measure_type,xl_quantile,xu_quantile,yl_quantile,yu_quantile);
    end
    
    %% bootstrap standard error and 95% confidence interval
    %% U-statistic-based estimator
    mean_b_nonparametric = mean(b_nonparametric_local_tau);
    % standard error
    se_nonparametric = sqrt(sum((b_nonparametric_local_tau(:,1) - mean_b_nonparametric).^2)./(B_times - 1));
    b_se_nonparametric(mc_index,1) = se_nonparametric;
    % 95% confidence interval
    b_conf_interval_95_nonparametric(mc_index,1) = mc_nonparametric_local_tau(mc_index,1) - 1.96 * se_nonparametric;
    b_conf_interval_95_nonparametric(mc_index,2) = mc_nonparametric_local_tau(mc_index,1) + 1.96 * se_nonparametric;
    
    %% copula model-based estimator
    mean_b_copula_based = mean(b_copula_based_local_tau);
    % standard error
    se_copula_based = sqrt(sum((b_copula_based_local_tau(:,1) - mean_b_copula_based).^2)./(B_times - 1));
    b_se_copula_based(mc_index,1) = se_copula_based;
    % 95% confidence interval
    b_conf_interval_95_copula_based(mc_index,1) = mc_copula_based_local_tau(mc_index,1) - 1.96 * se_copula_based;
    b_conf_interval_95_copula_based(mc_index,2) = mc_copula_based_local_tau(mc_index,1) + 1.96 * se_copula_based;  
end

%% true value of local Kendall's tau
ture_value = fun_copulald_type_II(copula_type,0,0,copula_par,0,0,measure_type,ld_type,x_quantile,y_quantile);
% ture_value = fun_copulald_general(copula_type,0,0,copula_par,0,0,measure_type,xl_quantile,xu_quantile,yl_quantile,yu_quantile);

%% simulation-based bias
% U-statistic estimator
mean_mc_nonparametric = mean(mc_nonparametric_local_tau);
mc_bias_nonparametric = mean_mc_nonparametric - ture_value;
% copula model-based estimator
mean_mc_copula_based_estimator = mean(mc_copula_based_local_tau);
mc_bias_copula_based_estimator = mean_mc_copula_based_estimator - ture_value;

%% simulation-based standard error
% U-statistic estimator
mc_se_nonparametric_estimator = sqrt(sum((mc_nonparametric_local_tau(:,1) - mean_mc_nonparametric).^2)./(mc_rep_times - 1));
% copula model-based estimator
mc_se_copula_based_estimator = sqrt(sum((mc_copula_based_local_tau(:,1) - mean_mc_copula_based_estimator).^2)./(mc_rep_times - 1));

%% average bootstrap standard error
% U-statistic estimator
ave_bootstrap_se_nonparametric_estimator = mean(b_se_nonparametric);
% copula model-based estimator
ave_bootstrap_se_copula_based_estimator = mean(b_se_copula_based);

%% empirical 95% coverage probability
%% U-statistic estimator
correct_number = 0;
for i = 1:1:mc_rep_times
    % check whether ture_value is in the 95% confidence interval
    if ture_value >= b_conf_interval_95_nonparametric(i,1) && ture_value <= b_conf_interval_95_nonparametric(i,2)
        correct_number = correct_number + 1;
    end
end
emp_95_cp_nonparametric = correct_number/mc_rep_times;

%% copula model-based estimator
correct_number = 0;
for i = 1:1:mc_rep_times
    if ture_value >= b_conf_interval_95_copula_based(i,1) && ture_value <= b_conf_interval_95_copula_based(i,2)
        correct_number = correct_number + 1;
    end
end
emp_95_cp_copula_based = correct_number/mc_rep_times;


