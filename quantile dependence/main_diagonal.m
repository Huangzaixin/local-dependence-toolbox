%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        empirical and theoretical quantile Kendall's tau along the main diagonal        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load uniform (0,1) dataset
u = csvread('quantile dependence/data/simulated_u.csv',1,0); 
v = csvread('quantile dependence/data/simulated_v.csv',1,0);  
UV = [u v];
Xtemp = UV(:,1);
Ytemp = UV(:,2);
 
%% 1. empirical quantile Kendall's tau along the main diagonal
quantile_parameter = 0.05:0.02:0.5;
for i=1:1:length(quantile_parameter)
       temp = find((UV(:,1)<=quantile(UV(:,1),quantile_parameter(i))&(UV(:,2)<=quantile(UV(:,2),quantile_parameter(i))))); 
       X = Xtemp(temp);
       Y = Ytemp(temp); 
       if ~isempty(X) && ~isempty(Y)
          mainlocalKendallTau_l(i) = corr(X,Y,'type','Kendall');
       end   
end
plot(quantile_parameter,mainlocalKendallTau_l,'-o','MarkerSize',8,'LineWidth',1); hold on;

quantile_parameter_u = 0.5:0.02:0.95;
for i=1:1:length(quantile_parameter_u)
       temp = find((UV(:,1)>=quantile(UV(:,1),quantile_parameter_u(i))&(UV(:,2)>=quantile(UV(:,2),quantile_parameter_u(i))))); 
       X = Xtemp(temp);
       Y = Ytemp(temp);
       if ~isempty(X) && ~isempty(Y)
          mainlocalKendallTau_u(i) = corr(X,Y,'type','Kendall');
       end    
end 
plot(quantile_parameter_u,mainlocalKendallTau_u,'-o','MarkerSize',8,'LineWidth',1); hold on;

y=[mainlocalKendallTau_l(23) mainlocalKendallTau_u(1)];
x=[0.5 0.5];
plot(x,y,'-','MarkerSize',10,'LineWidth',1,'Color',[0 0 1]); hold on;

% empirical global Kendall's tau
[emp_globalKendalltau,PVAL] = corr(Xtemp,Ytemp,'type','Kendall');

x=[0 1];
y=[emp_globalKendalltau emp_globalKendalltau];
plot(x,y,'-','MarkerSize',10,'LineWidth',2,'Color',[0 0 1]);


%% 2. copula model-based quantile Kendall's tau curve along the main diagonal
% copula parameters
copulatype = 'mixFG';                   % copula type, see "functions/fun_copulald_general.m"
weight1 = 0.51552102291413;             % the weight of the first copula function in the mixture copula model
weight2 = 0;                            % the weight of the second copula function in the mixture copula model
copulaparameter1 = 1.00109966476764;    % for mixture copula model, parameter of the first copula function; for sjc copula, the upper tail dependence coefficient    
copulaparameter2 = 2.03324976306063;    % for mixture copula model, parameter of the second copula function; for sjc copula, the lower tail dependence coefficient 
copulaparameter3 = 0;                   % for mixture copula model, parameter of the third copula function
measuretype = 'Kendall';                % Kendall's tau

% theoretical quantile Kendall's tau along the main diagonal
quantile_parameter_l = 0.05:0.01:0.5;
for i=1:1:length(quantile_parameter_l) 
       mainlocalKendallTau_l(i) = fun_copulald_type_II(copulatype,weight1,weight2,copulaparameter1,copulaparameter2,copulaparameter3,'kendall','ll',quantile_parameter_l(i),quantile_parameter_l(i));
end
plot(quantile_parameter_l,mainlocalKendallTau_l,'-','MarkerSize',8,'LineWidth',3,'Color',[1 0 0]);hold on;

quantile_parameter_u = 0.5:0.01:0.95;
for j=1:1:length(quantile_parameter_u)
       mainlocalKendallTau_u(j) = fun_copulald_type_II(copulatype,weight1,weight2,copulaparameter1,copulaparameter2,copulaparameter3,'kendall','uu',quantile_parameter_u(j),quantile_parameter_u(j));
end
plot(quantile_parameter_u,mainlocalKendallTau_u,'-','MarkerSize',8,'LineWidth',3,'Color',[1 0 0]);hold on;

y=[mainlocalKendallTau_l(46) mainlocalKendallTau_u(1)];
x=[0.5 0.5];
plot(x,y,'-','MarkerSize',10,'LineWidth',2,'Color',[1 0 0]);
hold on;

% theoretical global Kendall's tau
globalKendalltau = fun_copulald_type_II(copulatype,weight1,weight2,copulaparameter1,copulaparameter2,copulaparameter3,'kendall','ll',1,1);

x=[0 1];
y=[globalKendalltau globalKendalltau];
plot(x,y,'-','MarkerSize',10,'LineWidth',2,'Color',[1 0 0]);


