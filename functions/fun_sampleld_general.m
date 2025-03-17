function [ld] = fun_sampleld_general(X,Y,type,pl,pu,ql,qu)
% Description: Type I local Kendall's tau for sample data
% Inputs:  
%       1. X and Y: sample data
%       2. type:
%          'Kendall' : local Kendall's tau
%       3. pl: the lower bound quantile of variable X
%          pu: the upper bound quantile of variable X  
%          ql: the lower bound quantile of variable Y
%          qu: the upper bound quantile of variable Y
% Output: ld
%       1. ld: Type I local Kendall's tau of sample in the selected region
% Author: Zaixin Huang
% Date: finished at 2023.01.01; this version: 2025.03.16
% Bug reports and suggestions: 
%       if you find any bugs or have suggestions, please contact me at eric.huangzaixin@gmail.com. 
%       I will update them on GitHub and acknowledge your contribution. Thank you!
% The latest version can be downloaded from https://github.com/huangzaixin/local-dependence-toolbox
%%
Z(:,1) = X;
Z(:,2) = Y;

temp = find( quantile(Z(:,1),pu) >= Z(:,1) & Z(:,1) > quantile(Z(:,1),pl) ...,
               & quantile(Z(:,2),qu) >= Z(:,2) & Z(:,2) > quantile(Z(:,2),ql) ); 
 
Xtemp = X(temp);
Ytemp = Y(temp);

if isempty(Xtemp) || isempty(Ytemp)
    ld = 0;
    return;
end
        
switch lower(type)
    case 'kendall' 
        ld = corr(Xtemp,Ytemp,'type','Kendall'); 
    otherwise
        error('Error: unknown measure type.');
end


