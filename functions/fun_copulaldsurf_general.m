function [ldM] = fun_copulaldsurf_general(copulatype,weight1,weight2,copulaparameter1,copulaparameter2,copulaparameter3,measuretype,quantile_interval)
% Description: Type I local Kendall's tau surface of a copula model
% Inputs: 
%      1. copulatype: see functions/fun_copulald_general.m
%      2. weight1: the weight of the first copula function in the mixture copula model
%         weight2: the weight of the second copula function in the mixture copula model
%      3. copulaparameter1: parameter of the first copula function
%         copulaparameter2: parameter of the second copula function
%         copulaparameter3: parameter of the third copula function
%      4. measuretype: 
%         'kendall' : local Kendall's tau surface
%      5. quantile_interval: 
%         the width (or height) of the square region
% Output: ldM
%      1. ldM: Type I local Kendall's tau matrix
% Author: Zaixin Huang
% Date: finished at 2023.01.01; this version: 2025.03.16
% Bug reports and suggestions: 
%       if you find any bugs or have suggestions, please contact me at eric.huangzaixin@gmail.com. 
%       I will update them on GitHub and acknowledge your contribution. Thank you!
% The latest version can be downloaded from https://github.com/huangzaixin/local-dependence-toolbox
%%
% the square region interval
switch quantile_interval
    case 0.05
        quantile_upper = 0.95;
    case 0.1
        quantile_upper = 0.9;
    case 0.2
        quantile_upper = 0.8;
    case 0.25
        quantile_upper = 0.75;
    case 0.5
        quantile_upper = 0.5;
end

quantile_X = 0:quantile_interval:quantile_upper;
quantile_Y = 0:quantile_interval:quantile_upper;
        
for i=1:1:length(quantile_X) 
    for j=1:1:length(quantile_Y)
        if quantile_X(i) < 1 || quantile_Y(j) < 1
             temp_ld = fun_copulald_general(copulatype,weight1,weight2,copulaparameter1,copulaparameter2,copulaparameter3,measuretype,...,
                                           quantile_X(i),quantile_X(i) + quantile_interval,...,
                                           quantile_Y(j),quantile_Y(j) + quantile_interval);
             localdependence(i,j) = temp_ld;
        end
    end
end

h = surf(localdependence); 
h.FaceColor = 'interp';
% colormap jet;  
colorbar; 
view(-30,30);

xlabel('u','FontSize',13); 
ylabel('v','FontSize',13);
set(gca,'FontSize',13);
set(gcf,'color','w') 
set(gcf,'unit','centimeters','position',[8 18 10 8]);

switch lower(measuretype)
    case 'kendall'
         zlabel('theoretical local Kendall''s ¦Ó','FontSize',12);        
    otherwise
       error(message('Unknown measure type.'));
end 

switch quantile_interval
    case 0.05
        set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21]);
        set(gca,'YTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21]);
        set(gca,'XTickLabel',{'0',' ',' ',' ','0.2',' ',' ',' ','0.4',' ',' ',' ','0.6',' ',' ',' ','0.8',' ',' ',' ','1'});
        set(gca,'YTickLabel',{'0',' ',' ',' ','0.2',' ',' ',' ','0.4',' ',' ',' ','0.6',' ',' ',' ','0.8',' ',' ',' ','1'});   
    case 0.1
        set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11]);
        set(gca,'YTick',[1 2 3 4 5 6 7 8 9 10 11]);
        set(gca,'XTickLabel',{'0',' ','0.2',' ','0.4',' ','0.6',' ','0.8',' ','1'});
        set(gca,'YTickLabel',{'0',' ','0.2',' ','0.4',' ','0.6',' ','0.8',' ','1' }); 
    case 0.2
        set(gca,'XTick',[1 2 3 4 5 6]);
        set(gca,'YTick',[1 2 3 4 5 6]);
        set(gca,'XTickLabel',{'0','0.2','0.4','0.6','0.8','1'});
        set(gca,'YTickLabel',{'0','0.2','0.4','0.6','0.8','1' });
    case 0.25
        set(gca,'XTick',[1 2 3 4 5]);
        set(gca,'YTick',[1 2 3 4 5]);
        set(gca,'XTickLabel',{'0','0.25','0.5','0.75','1'});
        set(gca,'YTickLabel',{'0','0.25','0.5','0.75','1'});   
    case 0.5
        set(gca,'XTick',[1 2 3]);
        set(gca,'YTick',[1 2 3]);
        set(gca,'XTickLabel',{'0','0.5','1'});
        set(gca,'YTickLabel',{'0','0.5','1'});    
end

ldM = localdependence;



