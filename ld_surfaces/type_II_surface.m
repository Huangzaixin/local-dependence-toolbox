%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Type II Local Kendall's Tau Surface                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: This program computes and plots empirical and theoretical    %
%              Type II local Kendall's tau surfaces for copula models       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic Settings
measure_type = 'Kendall';    % dependence measure type: 'Kendall'
region_type = 'uu';          % region type: 
                             % "uu": upper-upper local Kendall's tau; 
                             % "ll": lower-lower local Kendall's tau;
                             % "ul": upper-lower local Kendall's tau; 
                             % 'lu': lower-upper local Kendall's tau;
copula_type = 'clayton';     % copula type: 'clayton', 'gumbel', 'frank', etc.  
copula_parameter = 1.3;      % copula parameter 
sample_size = 10000;         % number of samples to generate

%% 1. Empirical type II local Kendall's tau surface
% U = csvread('data_examples/example_2/data/u_cases.csv',1,0);
% V = csvread('data_examples/example_2/data/v_deaths.csv',1,0);

% Generate random data from the specified copula
UV = copularnd(copula_type,copula_parameter,sample_size);   
U = UV(:,1); 
V = UV(:,2);

ldMatrix1 = fun_sampleldsurf_type_II(U,V,measure_type,region_type); hold on;


%% 2. Copula Model-Based Type II Local Kendall's Tau Surface
copulatype = copula_type;     % copula type, see functions/fun_copulald_general.m
weight1 = 0;                  % the first copula function's weight in the mixture copula model
weight2 = 0;                  % the second copula function's weight in the mixture copula model
copulaparameter1 = 1.3;       % for mixture copula model, parameter of the first copula function; for SJC copula, the upper tail dependence coefficient   
copulaparameter2 = 0;         % for mixture copula model, parameter of the second copula function; for SJC copula, the lower tail dependence coefficient 
copulaparameter3 = 0;         % for mixture copula model, parameter of the third copula function

ldMatrix2 = fun_copulaldsurf_type_II(copulatype,weight1,weight2,copulaparameter1,copulaparameter2,copulaparameter3,measuretype,region_type);


