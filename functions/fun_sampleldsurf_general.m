function [ldM] = fun_sampleldsurf_general(X,Y,measuretype,quantile_interval)
% Description: Type I local Kendall's tau surface for sample data
% Inputs:  
%       1. X and Y: sample
%       2. measuretype:
%          'Kendall' : local Kendall's tau
%       3. quantile_interval: the width (or height) of the local region
%          We set quantile_interval = 0.05, 0.1, 0.2
% Output: ldM
%       1. ldM: the Type I local Kendall's tau matrix
% Author: Zaixin Huang
% Date: finished at 2023.01.01; this version: 2025.03.16
% Bug reports and suggestions: if you find any bugs or have suggestions, please contact me at eric.huangzaixin@gmail.com. 
%                              I will update them on GitHub and acknowledge your contribution. Thank you!
% The latest version can be downloaded from https://github.com/huangzaixin/local-dependence-toolbox
%%
Xtemp = X;
Ytemp = Y;

%% select region
% Note, when I = 0.05, the upper bounds of quantile_X and quantile_Y are 0.95
%       when I = 0.1 (default), the upper bounds of quantile_X and quantile_Y are 0.9
%       when I = 0.2, the upper bounds of quantile_X and quantile_Y are 0.8
switch quantile_interval
    case 0.05
        quantile_upper = 0.95;
    case 0.1
        quantile_upper = 0.9;
    case 0.2
        quantile_upper = 0.8;     
    otherwise
       error('unsupported measure type.');
end

quantile_X = 0:quantile_interval:quantile_upper;
quantile_Y = 0:quantile_interval:quantile_upper;
        
for i=1:1:length(quantile_X)
    for j=1:1:length(quantile_Y)
        LocalKendallMarix(i,j) = fun_sampleld_general(Xtemp, Ytemp, measuretype, ...,
                                    quantile_X(i), quantile_X(i) + quantile_interval, ...,
                                    quantile_Y(j), quantile_Y(j) + quantile_interval);
        LocalKendallMarix(i,j);
        disp('The program is running, please wait for a few seconds...');
    end
end 

bar3(LocalKendallMarix);
% surf(LocalKendallMarix);
view(-115,35);

xlabel('v','FontSize',12);
ylabel('u','FontSize',12);

switch lower(measuretype)
    case 'kendall'
         zlabel('local Kendall''s ¦Ó','FontSize',12);     
    otherwise
       error('unsupported measure type.');
end

switch quantile_interval
    case 0.05
        set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21]);
        set(gca,'YTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21]);
        set(gca,'XTickLabel',{'0',' ',' ',' ','0.2',' ',' ',' ','0.4',' ',' ',' ','0.6',' ',' ',' ','0.8',' ',' ',' ','1'});
        set(gca,'YTickLabel',{'0',' ',' ',' ','0.2',' ',' ',' ','0.4',' ',' ',' ','0.6',' ',' ',' ','0.8',' ',' ',' ','1' });   
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
end

% set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19]);
% set(gca,'YTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19]);
% set(gca,'XTickLabel',{'0.05',' ',' ',' ',' ',' ','0.35',' ',' ',' ',' ',' ','0.65',' ',' ',' ',' ',' ','0.95'});
% set(gca,'YTickLabel',{'0.05',' ',' ',' ',' ',' ','0.35',' ',' ',' ',' ',' ','0.65',' ',' ',' ',' ',' ','0.95'});   
    
ldM = LocalKendallMarix;



