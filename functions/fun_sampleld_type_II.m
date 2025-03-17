function [ld] = fun_sampleld_type_II(X,Y,type,region,Xquantile,Yquantile)
% Description: Type II local Kendall's tau for sample data
% Inputs:
%       1. X and Y: sample data
%       2. type: 
%          'Kendall' : local Kendall's tau
%       3. region:
%          'UU': upper-upper region; 'UL': upper-lower region 
%          'LU': lower-upper region; 'LL': lower-lower region 
%       4. Xquantile: quantile of variable X
%          Yquantile: quantile of variable Y
% Output: ld
%       1. ld: Type II local dependence of samples in the selected region
% Author: Zaixin Huang
% Date: finished at 2015.06.07; this version: 2025.03.16
% Bug reports and suggestions: 
%       if you find any bugs or have suggestions, please contact me at eric.huangzaixin@gmail.com. 
%       I will update them on GitHub and acknowledge your contribution. Thank you!
% The latest version can be downloaded from https://github.com/huangzaixin/local-dependence-toolbox
%%
Z(:,1) = X;
Z(:,2) = Y;
 
switch lower(region)
case 'uu'
    temp = find((Z(:,1) >= quantile(Z(:,1),Xquantile) & (Z(:,2) >= quantile(Z(:,2),Yquantile) ))); 
case 'ul'
    temp = find((Z(:,1) >= quantile(Z(:,1),Xquantile) & (Z(:,2) <= quantile(Z(:,2),Yquantile) ))); 
case 'lu'
    temp = find((Z(:,1) <= quantile(Z(:,1),Xquantile) & (Z(:,2) >= quantile(Z(:,2),Yquantile) ))); 
case 'll'
    temp = find((Z(:,1) <= quantile(Z(:,1),Xquantile) & (Z(:,2) <= quantile(Z(:,2),Yquantile) ))); 
otherwise
    error('Error: unknown region.');
end 
 
%% 
Xtemp = X(temp);
Ytemp = Y(temp); 

switch lower(type)
case 'kendall'
    if ~isempty(Xtemp)
        ld = corr(Xtemp,Ytemp,'type','Kendall');
    else 
        ld = 0;
    end
otherwise
    error('Error: unknown measure type.');
end



