function [ldM] = fun_copulaldsurf_type_II(copulatype,weight1,weight2,copulaparameter1,copulaparameter2,copulaparameter3,measuretype,region)
% Description: Type II local Kendall's tau surface of a copula model
% Inputs: 
%      1. copulatype: see functions/fun_copulald_type_II.m
%      2. weight1: the weight of the first copula function in the mixture copula model
%         weight2: the weight of the second copula function in the mixture copula model
%      3. copulaparameter1: parameter of the first copula function
%         copulaparameter2: parameter of the second copula function
%         copulaparameter3: parameter of the third copula function
%      4. measuretype:
%         'kendall': local Kendall's tau surface
%      5. region:
%         'UU': upper-upper region; 'UL': upper-lower region;
%         'LU': lower-upper region; 'LL': lower-lower region;
% Output: ldM
%      1. ldM: local Kendall's tau matrix in the selected region
% Author: Zaixin Huang
% Date: finished at 2015.06.07; updated at 2018.01.04; this version: 2025.03.16
% Bug reports and suggestions: 
%       if you find any bugs or have suggestions, please contact me at eric.huangzaixin@gmail.com. 
%       I will update them on GitHub and acknowledge your contribution. Thank you!
% The latest version can be downloaded from https://github.com/huangzaixin/local-dependence-toolbox
%%
quantile_X = 0.05:0.05:0.95;
quantile_Y = 0.05:0.05:0.95;
      
for i=1:1:length(quantile_X)
    for j=1:1:length(quantile_Y)
        localdependence(i,j) = fun_copulald_type_II(copulatype,weight1,weight2,copulaparameter1,copulaparameter2,copulaparameter3,measuretype,region,quantile_X(i),quantile_Y(j));            
    end
end

surf(localdependence, 'FaceColor',[0.7100, 0.2549, 0.0588]); % [0.87,0.53,0.39]
xlabel('u','FontSize',13);
ylabel('v','FontSize',13);
set(gcf,'color','w');
set(gca,'FontSize',13);
zlim([0 1]);
set(gcf,'unit','centimeters','position',[8 18 10 8]);

%%
switch lower(region)
    case 'uu'        
        set(gca,'XTick',[1 4 7 10 13 16 19]);
        set(gca,'YTick',[1 4 7 10 13 16 19]);
        set(gca,'XTickLabel',{'0.05',' ','0.35',' ','0.65',' ','0.95'});
        set(gca,'YTickLabel',{'0.05',' ','0.35',' ','0.65',' ','0.95'});  
     
        switch lower(measuretype)
            case 'kendall'
                 zlabel('upper-upper local Kendall''s 而 surface','FontSize',13);        
        otherwise
           error(message('unknown measure type.'));
        end 
        
    case 'ul'     
        set(gca,'XTick',[1 4 7 10 13 16 19]);
        set(gca,'YTick',[1 4 7 10 13 16 19]);
        set(gca,'XTickLabel',{'0.95',' ','0.65',' ','0.35',' ','0.05'});
        set(gca,'YTickLabel',{'0.95',' ','0.65',' ','0.35',' ','0.05'});
        
        switch lower(measuretype) 
            case 'kendall'
                 zlabel('upper-lower local Kendall''s 而 surface','FontSize',13);        
        otherwise
           error(message('unknown measure type.'));
        end 
        
    case 'lu'      
        set(gca,'XTick',[1 4 7 10 13 16 19]);
        set(gca,'YTick',[1 4 7 10 13 16 19]);
        set(gca,'XTickLabel',{'0.95',' ','0.65',' ','0.35',' ','0.05'});
        set(gca,'YTickLabel',{'0.95',' ','0.65',' ','0.35',' ','0.05'});
    
        switch lower(measuretype)
            case 'kendall'
                 zlabel('lower-upper local Kendall''s 而 surface','FontSize',13);        
        otherwise
           error(message('unknown measure type.'));
        end 
    
    case 'll'    
        set(gca,'XTick',[1 4 7 10 13 16 19]);
        set(gca,'YTick',[1 4 7 10 13 16 19]);
        set(gca,'XTickLabel',{'0.05',' ','0.35',' ','0.65',' ','0.95'});
        set(gca,'YTickLabel',{'0.05',' ','0.35',' ','0.65',' ','0.95'});  
     
        switch lower(measuretype)
            case 'kendall'
                 zlabel('lower-lower local Kendall''s 而 surface','FontSize',13);        
            otherwise
               error(message('unknown measure type.'));
        end 
    
    otherwise
       error('unknown selected region.');
end 

disp('The program is running...');

ldM = localdependence;



