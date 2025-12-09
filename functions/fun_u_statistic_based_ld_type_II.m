function ld = fun_u_statistic_based_ld_type_II(data,ld_type,X,Y)
% Description: Type II local Kendall's tau based on U-statistic
%
% Inputs: 
%       1. data: sample data
%       2. ld_type: 
%          'll': lower-lower local Kendall's tau; 'uu': upper-upper local Kendall's tau
%       3. X: for 'll', upper bound of variable X
%             for 'uu', lower bound of variable X
%       4. Y: for 'll', upper bound of variable Y
%             for 'uu', lower bound of variable Y
%
% Written for paper "Generalized local Kendall’s τ: a novel framework for uncovering nonlinear local dependence" (Huang & Zhang,2026)
%
% Author: Zaixin Huang
% Date: finished at 2015.06.07; updated at 2018.01.04; current version: 2025.03.16
% Contact: For bug reports and suggestions, please contact me at eric.huangzaixin@gmail.com. 
%          I will update them on GitHub and acknowledge your contribution. Thank you!
% The latest version can be downloaded from https://github.com/huangzaixin/local-dependence-toolbox
%%
if isempty(data)
    ld = NaN;
    return;
end

d_X = data(:,1);
d_Y = data(:,2);
total_sum = 0;
M = 0;   % number of observation pairs in the selected region  

switch lower(ld_type)
    case 'uu'
        for i=1:1:length(d_X)-1
            for j=i+1:1:length(d_X)
                total_sum = total_sum + sign((d_X(i)-d_X(j)).*(d_Y(i)-d_Y(j))).*(d_X(i) >= X & d_X(j) >= X & d_Y(i) >= Y & d_Y(j) >= Y);
                M = M + (d_X(i) >= X & d_X(j) >= X & d_Y(i) >= Y & d_Y(j) >= Y);
            end
        end
    case 'll'
        for i=1:1:length(d_X)-1
            for j=i+1:1:length(d_X)
                total_sum = total_sum + sign((d_X(i)-d_X(j)).*(d_Y(i)-d_Y(j))).*(X >= d_X(i) & X >= d_X(j) & Y >= d_Y(i) & Y >= d_Y(j));
                M = M + (X >= d_X(i) & X >= d_X(j) & Y >= d_Y(i) & Y >= d_Y(j));
            end
        end
end

if M == 0
    ld = NaN;
else
    ld = total_sum/M;   
end

end

