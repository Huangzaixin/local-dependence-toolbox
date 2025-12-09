function [ldM] = fun_sampleldsurf_type_II(X,Y,measuretype,region)
% Description: Type II local Kendall's tau surface for sample data
% Inputs:  
%       1. X and Y: sample data vectors
%       2. measuretype:  'Kendall' - local Kendall's tau
%       3. region:
%          'UU': upper-upper region; 'UL': upper-lower region
%          'LU': lower-upper region; 'LL': lower-lower region
%
% Outputs: ldM
%       1. ldM: Type II local Kendall's tau matrix
%
% Written for paper "Generalized local Kendall’s τ: a novel framework for uncovering nonlinear local dependence" (Huang & Zhang,2026)
%
% Author: Zaixin Huang
% Date: finished at 2015.06.07; updated at 2018.01.04; current version: 2025.03.16
% Contact: For bug reports and suggestions, please contact me at eric.huangzaixin@gmail.com. 
%          I will update them on GitHub and acknowledge your contribution. Thank you!
% The latest version can be downloaded from https://github.com/huangzaixin/local-dependence-toolbox
%% Sample data vectors
Xtemp = X;
Ytemp = Y;

%% Initialize quantile grid and output matrix
quantile_X = 0.05:0.05:0.95;
quantile_Y = 0.05:0.05:0.95;
ldM = zeros(length(quantile_X), length(quantile_Y));

for i=1:length(quantile_X)
    for j=1:length(quantile_Y)
        quantile_x = quantile_X(i);
        quantile_y = quantile_Y(j);
        ldM(i,j) = fun_sampleld_type_II(Xtemp,Ytemp,measuretype,region,quantile_x,quantile_y); 
        ldM(i,j);
    end
end

%% Local Kendall's tau surface 
surf(ldM, 'FaceColor',[0, 0.4471, 0.7412]); % [0.33,0.65,0.68]

xlabel('u','FontSize',13);
ylabel('v','FontSize',13);

switch lower(measuretype)
    case 'kendall'
         zlabel('empirical local Kendall''s tau','FontSize',13);     
    otherwise
       error('unknown measure type');
end 

%%
set(gcf,'color','w');
set(gca,'FontSize',13);
zlim([0 1]); 
set(gcf,'unit','centimeters','position',[8 18 10 8]);

switch lower(region)
    case 'uu'
        view(-30,15);
        set(gca,'XTick',[1 4 7 10 13 16 19]);
        set(gca,'YTick',[1 4 7 10 13 16 19]);
        set(gca,'XTickLabel',{'0.05',' ','0.35',' ','0.65',' ','0.95'});
        set(gca,'YTickLabel',{'0.05',' ','0.35',' ','0.65',' ','0.95'});  
  
    case 'ul'  
        view(150,15);
        set(gca,'XTick',[1 4 7 10 13 16 19]);
        set(gca,'YTick',[1 4 7 10 13 16 19]);
        set(gca,'XTickLabel',{'0.95',' ','0.65',' ','0.35',' ','0.05'});
        set(gca,'YTickLabel',{'0.95',' ','0.65',' ','0.35',' ','0.05'});
        
    case 'lu'
        view(150,15);
        set(gca,'XTick',[1 4 7 10 13 16 19]);
        set(gca,'YTick',[1 4 7 10 13 16 19]);
        set(gca,'XTickLabel',{'0.95',' ','0.65',' ','0.35',' ','0.05'});
        set(gca,'YTickLabel',{'0.95',' ','0.65',' ','0.35',' ','0.05'});
    
    case 'll'
        view(-30,15);
        set(gca,'XTick',[1 4 7 10 13 16 19]);
        set(gca,'YTick',[1 4 7 10 13 16 19]);
        set(gca,'XTickLabel',{'0.05',' ','0.35',' ','0.65',' ','0.95'});
        set(gca,'YTickLabel',{'0.05',' ','0.35',' ','0.65',' ','0.95'});  
        
    otherwise
       error('unknown region.');
end 

disp('The program is running, please wait for a few seconds...');


