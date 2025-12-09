function ld = fun_u_statistic_based_ld_general(data,Xl,Xu,Yl,Yu)
% Description: Type I local Kendall's tau based on U-statistic
% Inputs:  
%       1. data: sample data
%       2. Xl: lower bound of variable X
%          Xu: upper bound of variable X
%          Yl: lower bound of variable Y
%          Yu: upper bound of variable Y
%
% Written for paper "Generalized local Kendall’s τ: a novel framework for uncovering nonlinear local dependence" (Huang & Zhang,2026)
%
% Author: Zaixin Huang
% Date: finished at 2023.01.01; current version: 2025.03.16
% Contact: For bug reports and suggestions, please contact me at eric.huangzaixin@gmail.com. 
%          I will update them on GitHub and acknowledge your contribution. Thank you!
% The latest version can be downloaded from https://github.com/huangzaixin/local-dependence-toolbox
%%
if Xu <= Xl || Yu <= Yl
    error('Upper bounds must be greater than lower bounds.');
end

d_X = data(:,1);
d_Y = data(:,2);     
total_sum = 0;
M = 0;   % number of observation pairs in the selected region  

for i=1:length(d_X)-1
    for j=i+1:length(d_X)
        total_sum = total_sum + sign((d_X(i) - d_X(j)).*(d_Y(i) - d_Y(j))).*(Xu >= d_X(i) & d_X(i) >= Xl & Xu >= d_X(j) & d_X(j) >= Xl & Yu >= d_Y(i) & d_Y(i) >= Yl & Yu >= d_Y(j) & d_Y(j) >= Yl);
        M = M + (Xu >= d_X(i) & d_X(i) >= Xl & Xu >= d_X(j) & d_X(j) >= Xl & Yu >= d_Y(i) & d_Y(i) >= Yl & Yu >= d_Y(j) & d_Y(j) >= Yl);
    end
end

if M == 0
    ld = NaN;
else
    ld = total_sum / M;   
end   

end



