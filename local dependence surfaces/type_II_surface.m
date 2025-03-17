%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       Type II local Kendall's tau surface       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% basic configuration
measuretype = 'Kendall';  
% "uu": upper-upper local Kendall's tau 
% "ll": lower-lower local Kendall's tau 
% "ul": upper-lower local Kendall's tau
% "lu": lower-upper local Kendall's tau
region_type = 'uu';

%% 1. empirical Type II local Kendall's tau surface
% U = csvread('data examples/example_2/data/u_cases.csv',1,0);
% V = csvread('data examples/example_2/data/v_deaths.csv',1,0);

UV = copularnd("clayton",1.3,10000);   
U = UV(:,1); 
V = UV(:,2);

ldMatrix1 = fun_sampleldsurf_type_II(U,V,measuretype,region_type); hold on;

%% 2. copula model-based theoretical Type II local Kendall's tau surface
copulatype = 'clayton';       % copula type, see functions/fun_copulald_general.m
weight1 = 0;                  % the weight of the first copula function in the mixture copula model
weight2 = 0;                  % the weight of the second copula function in the mixture copula model
copulaparameter1 = 1.3;       % for mixture copula model, parameter of the first copula function; for SJC copula, the upper tail dependence coefficient   
copulaparameter2 = 0;         % for mixture copula model, parameter of the second copula function; for SJC copula, the lower tail dependence coefficient 
copulaparameter3 = 0;         % for mixture copula model, parameter of the third copula function

ldMatrix2 = fun_copulaldsurf_type_II(copulatype,weight1,weight2,copulaparameter1,copulaparameter2,copulaparameter3,measuretype,region_type);


