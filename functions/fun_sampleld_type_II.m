function [ld] = fun_sampleld_type_II(X,Y,type,region,Xquantile,Yquantile)
% Description: Type II local Kendall's tau for sample data
% Inputs:
%       1. X and Y: sample data vectors
%       2. type: 
%          'Kendall' : local Kendall's tau
%       3. region:
%          'UU': upper-upper region; 'UL': upper-lower region 
%          'LU': lower-upper region; 'LL': lower-lower region 
%       4. Xquantile: quantile of variable X
%          Yquantile: quantile of variable Y
%
% Outputs: ld
%       1. ld: Type II local dependence of samples in the selected region
%
% Written for paper "Generalized local Kendall’s τ: a novel framework for uncovering nonlinear local dependence" (Huang & Zhang,2026)
%
% Author: Zaixin Huang
% Date: finished at 2015.06.07; current version: 2025.03.16
% Contact: For bug reports and suggestions, please contact me at eric.huangzaixin@gmail.com. 
%          I will update them on GitHub and acknowledge your contribution. Thank you!
% The latest version can be downloaded from https://github.com/huangzaixin/local-dependence-toolbox
%%
if Xquantile <= 0 || Xquantile >= 1 || Yquantile <= 0 || Yquantile >= 1
    error('Error: quantiles must be between 0 and 1.');
end

% Calculate quantile thresholds
X_threshold = quantile(X, Xquantile);
Y_threshold = quantile(Y, Yquantile);

% Find data points in the specified region
switch lower(region)
case 'uu'
    temp = find((X >= X_threshold) & (Y >= Y_threshold)); 
case 'ul'
    temp = find((X >= X_threshold) & (Y <= Y_threshold)); 
case 'lu'
    temp = find((X <= X_threshold) & (Y >= Y_threshold));
case 'll'
    temp = find((X <= X_threshold) & (Y <= Y_threshold)); 
otherwise
    error('Error: unknown region. Please use ''UU'', ''UL'', ''LU'', or ''LL''.');
end 
 
%% 
Xtemp = X(temp);
Ytemp = Y(temp); 

switch lower(type)
case 'kendall'
    try
        if ~isempty(Xtemp)
            ld = corr(Xtemp,Ytemp,'type','Kendall');
            if isnan(ld)
                ld = 0;
            end
        else 
            ld = 0;
        end
    catch
        warning('Error calculating Kendall''s tau. Returning 0.');
        ld = 0;
    end
otherwise
    error('Error: unknown measure type.');
end



