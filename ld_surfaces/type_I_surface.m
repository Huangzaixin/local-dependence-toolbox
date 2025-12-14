%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Type I Local Kendall's Tau Surface                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: This program computes and plots empirical and theoretical     %
%              Type I local Kendall's tau surfaces for copula models         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic Settings
measure_type = 'kendall';
quantile_interval = 0.05;    % width of the square region: 0.05, 0.1, 0.2
copula_type = 'clayton';     % copula type: 'clayton', 'gumbel', 'frank', etc.  
copula_parameter = 1.3;      % copula parameter 
sample_size = 1000000;       % number of samples

%% 1. Empirical Type I local Kendall's tau surface
% sample data
% X = csvread('data_examples/example_2/data/u_cases.csv',1,0);
% Y = csvread('data_examples/example_2/data/v_deaths.csv',1,0);

% Generate random data from the specified copula
data = copularnd(copula_type,copula_parameter,sample_size);
X = data(:,1);
Y = data(:,2);

% Draw empirical Type I local Kendall's tau surface
ldMatrix1 = fun_sampleldsurf_general(X,Y,measure_type,quantile_interval); hold on;


%% 2. Copula model-based Type I local Kendall's tau surface
% copula parameters
copulatype = copula_type;      % copula type, see functions/fun_copulald_general.m
weight1 = 0;                   % the first copula function's weight in the mixture copula model
weight2 = 0;                   % the second copula function's weight in the mixture copula model
copula_parameter1 = 1.3;       % for mixture copula model, parameter of the first copula function; for sjc copula, the upper tail dependence coefficient   
copula_parameter2 = 0;         % for mixture copula model, parameter of the second copula function; for sjc copula, the lower tail dependence coefficient 
copula_parameter3 = 0;         % for mixture copula model, parameter of the third copula function

ldMatrix2 = fun_copulaldsurf_general(copulatype,weight1,weight2,copulaparameter1,copulaparameter2,copulaparameter3,measure_type,quantile_interval); hold on;  


%% Add auxiliary lines
l = fun_copulald_general(copulatype,weight1,weight2,copula_parameter1,copula_parameter2,copula_parameter3,measure_type,0,quantile_interval,0,quantile_interval);
u = fun_copulald_general(copulatype,weight1,weight2,copula_parameter1,copula_parameter2,copula_parameter3,measure_type,1-quantile_interval,1,1-quantile_interval,1);

% quantile_interval = 0.05
plot3([21,0],[0,0],[0,0],'-k');
plot3([0,0],[21,0],[0,0],'-k');
plot3([1,1],[1,21],[l,l],'--k');
plot3([1,1],[1,1],[l,0],'--k');

plot3([20,20],[20,21],[u,u],'--k');
plot3([20,20],[20,20],[u,0],'--k');
plot3([20,21],[20,20],[u,u],'--k');
plot3([20,20],[20,20],[u,0],'--k');

% quantile_interval = 0.1
plot3([11,0],[0,0],[0,0],'-k');
plot3([0,0],[11,0],[0,0],'-k');
plot3([1,1],[1,11],[l,l],'--k');
plot3([1,1],[1,1],[l,0],'--k');

plot3([10,10],[10,11],[u,u],'--k');
plot3([10,10],[10,10],[u,0],'--k');
plot3([10,11],[10,10],[u,u],'--k');
plot3([10,10],[10,10],[u,0],'--k');

% quantile_interval = 0.2
plot3([6,0],[0,0],[0,0],'-k');
plot3([0,0],[6,0],[0,0],'-k');
plot3([1,1],[1,6],[l,l],'--k');
plot3([1,1],[1,1],[l,0],'--k');

plot3([5,5],[5,6],[u,u],'--k');
plot3([5,5],[5,5],[u,0],'--k');
plot3([5,6],[5,5],[u,u],'--k');
plot3([5,5],[5,5],[u,0],'--k');


