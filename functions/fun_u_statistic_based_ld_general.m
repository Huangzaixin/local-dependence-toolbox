function ld = fun_u_statistic_based_ld_general(data,Xl,Xu,Yl,Yu)
% Description: Type I local Kendall's tau based on U-statistic
% Inputs: 
%      1. data: sample data
%      2. Xl: lower bound of variable X
%         Xu: upper bound of variable X
%         Yl: lower bound of variable Y
%         Yu: upper bound of variable Y
% Author: Zaixin Huang
% Date: finished at 2023.01.01; this version: 2025.03.16
% Bug reports and suggestions: 
%      if you find any bugs or have suggestions, please contact me at eric.huangzaixin@gmail.com. 
%      I will update them on GitHub and acknowledge your contribution. Thank you!
% The latest version can be downloaded from https://github.com/huangzaixin/local-dependence-toolbox
%%
    d_X = data(:,1);
    d_Y = data(:,2);     
    sum = 0;
    M = 0;   % number of observation pairs in the selected region  
    
    for i=1:1:length(d_X)-1
        for j=i+1:1:length(d_X)
            sum = sum + sign((d_X(i) - d_X(j)).*(d_Y(i) - d_Y(j))).*(Xu >= d_X(i) & d_X(i) >= Xl & Xu >= d_X(j) & d_X(j) >= Xl & Yu >= d_Y(i) & d_Y(i) >= Yl & Yu >= d_Y(j) & d_Y(j) >= Yl);
            M = M + (Xu >= d_X(i) & d_X(i) >= Xl & Xu >= d_X(j) & d_X(j) >= Xl & Yu >= d_Y(i) & d_Y(i) >= Yl & Yu >= d_Y(j) & d_Y(j) >= Yl);
        end
    end
            
    ld = sum/M;   
end



