%% Comparison of global Kendall’s tau and generalized local Kendall’s tau for Gumbel, Clayton, Frank and FGM copula models
%% basic setup
clear;
measuretype = "kendall"; 
copulatype = 'fgm';          % "gumbel", "clayton", "frank", "fgm"
quantile_interval = 0.5;     % the square region interval
% copula parameter
% for Gumbel copula, theta = 1.1, 1.25, 1.5, 3
% for Clayton copula, theta = 0.1, 0.2, 0.5, 2
% for Frank copula, theta = 0.3, 3, 8, 15
% for FGM copula, theta = 0.1, 0.3, 0.5, 0.7
theta = 0.7;

%% theoretical global Kendall's tau
switch lower(copulatype)
case "gumbel"    % gumbel copula
    copula_globalKendalltau = 1 - 1/theta;    
case "clayton"   % clayton copula
    copula_globalKendalltau = theta/(theta+2);   
case "frank"     % frank copula
    copula_globalKendalltau = fun_copulald_general('frank',0,0,theta,0,0,measuretype,0,1,0,1);   
case "fgm"       % FGM copula
    copula_globalKendalltau = 2/9*theta;   
end

%% calculate generalized local Kendall’s tau in each region
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
             temp_ld = fun_copulald_general(copulatype,0,0,theta,0,0,measuretype,...,
                                           quantile_X(i),quantile_X(i) + quantile_interval,...,
                                           quantile_Y(j),quantile_Y(j) + quantile_interval);
             ldMaxtrix(i,j) = temp_ld;
        end
    end
end

%% 3d bar chart
b = bar3(ldMaxtrix,1);
view(-135,35);

for k = 1:length(b)
    zdata = b(k).ZData; 
    b(k).CData = zdata;  
    b(k).FaceColor = 'interp';
end

colorbar; % colormap jet;

zlim([0 0.2]);
zticks(0:0.05:0.2);

zlim([0 0.1]);
zticks(0:0.025:0.1);

zlim([0 0.16]);
zticks(0:0.04:0.16);

zlim([0 0.3]);
zticks(0:0.075:0.3);

zlim([0 0.04]);
zticks(0:0.01:0.04);

zlim([0 0.4]);
zticks(0:0.1:0.4);

zlim([0 0.44]);
zticks(0:0.11:0.44);

zlim([0 0.48]);
zticks(0:0.12:0.48);

zlim([0 0.6]);
zticks(0:0.15:0.6);

zlim([0 0.8]);
zticks(0:0.2:0.8);

zlim([0 1]);
zticks(0:0.25:1);

xlabel('v','FontSize',11); 
ylabel('u','FontSize',11);
set(gca,'FontSize',11);
set(gcf,'color','w') 
set(gcf,'unit','centimeters','position',[6 18 10 7]); 

%% XTickLabel, YTickLabel
switch quantile_interval
    case 0.05
        set(gca,'XTick',[0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5 11.5 12.5 13.5 14.5 15.5 16.5 17.5 18.5 19.5 20.5]);
        set(gca,'YTick',[0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5 11.5 12.5 13.5 14.5 15.5 16.5 17.5 18.5 19.5 20.5]);
        set(gca,'XTickLabel',{'0',' ',' ',' ','0.2',' ',' ',' ','0.4',' ',' ',' ','0.6',' ',' ',' ','0.8',' ',' ',' ','1'});
        set(gca,'YTickLabel',{'0',' ',' ',' ','0.2',' ',' ',' ','0.4',' ',' ',' ','0.6',' ',' ',' ','0.8',' ',' ',' ','1' });   
    case 0.1
        set(gca,'XTick',[0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5]);
        set(gca,'YTick',[0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5]);
        set(gca,'XTickLabel',{'0',' ','0.2',' ','0.4',' ','0.6',' ','0.8',' ','1'});
        set(gca,'YTickLabel',{'0',' ','0.2',' ','0.4',' ','0.6',' ','0.8',' ','1' }); 
    case 0.2
        set(gca,'XTick',[0.5 1.5 2.5 3.5 4.5 5.5]);
        set(gca,'YTick',[0.5 1.5 2.5 3.5 4.5 5.5]);
        set(gca,'XTickLabel',{'0','0.2','0.4','0.6','0.8','1'});
        set(gca,'YTickLabel',{'0','0.2','0.4','0.6','0.8','1' });
    case 0.25
        set(gca,'XTick',[0.5 1.5 2.5 3.5 4.5]);
        set(gca,'YTick',[0.5 1.5 2.5 3.5 4.5]);
        set(gca,'XTickLabel',{'0','0.25','0.5','0.75','1'});
        set(gca,'YTickLabel',{'0','0.25','0.5','0.75','1' }); 
    case 0.5
        set(gca,'XTick',[0.5 1.5 2.5]);
        set(gca,'YTick',[0.5 1.5 2.5]);
        set(gca,'XTickLabel',{'0','0.5','1'});
        set(gca,'YTickLabel',{'0','0.5','1' });     
end

%%
hold on;
switch quantile_interval
    case 0.05
        x = [0.5, 20.5, 20.5, 0.5];
        y = [0.5, 0.5, 20.5, 20.5];
        z = [copula_globalKendalltau, copula_globalKendalltau, copula_globalKendalltau, copula_globalKendalltau];
        fill3(x, y, z, [0.6, 0.6, 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

        plot3([0.5,20.5],[0.5,0.5],[copula_globalKendalltau,copula_globalKendalltau],'--k','LineWidth', 2, 'Color', [0.85,0.33,0.10]);
        plot3([0.5,0.5],[0.5,20.5],[copula_globalKendalltau,copula_globalKendalltau],'--k','LineWidth', 2, 'Color', [0.85,0.33,0.10]);
        plot3([0.5,20.5],[20.5,20.5],[copula_globalKendalltau,copula_globalKendalltau],'--k','LineWidth', 2, 'Color', [0.85,0.33,0.10]);
        plot3([20.5,20.5],[0.5,20.5],[copula_globalKendalltau,copula_globalKendalltau],'--k','LineWidth', 2, 'Color', [0.85,0.33,0.10]);
   
        %% lines: for gumbel copula model
        plot3([19.5,20.5],[19.5,19.5],[copula_globalKendalltau,copula_globalKendalltau],'-','LineWidth', 0.5, 'Color', [0,0,0]);
        plot3([19.5,19.5],[19.5,20.5],[copula_globalKendalltau,copula_globalKendalltau],'-','LineWidth', 0.5, 'Color', [0,0,0]);

        % lines: for Clayton copula model
        % plot3([0.5,1.5],[0.5,0.5],[copula_globalKendalltau,copula_globalKendalltau],'-','LineWidth', 0.5, 'Color', [0,0,0]);
        % plot3([0.5,0.5],[0.5,1.5],[copula_globalKendalltau,copula_globalKendalltau],'-','LineWidth', 0.5, 'Color', [0,0,0]);

    case 0.1
        x = [0.5, 10.5, 10.5, 0.5];
        y = [0.5, 0.5, 10.5, 10.5];
        z = [copula_globalKendalltau, copula_globalKendalltau, copula_globalKendalltau, copula_globalKendalltau];
        fill3(x, y, z, [0.6, 0.6, 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

        plot3([0.5,10.5],[0.5,0.5],[copula_globalKendalltau,copula_globalKendalltau],'--k','LineWidth', 2, 'Color', [0.85,0.33,0.10]);
        plot3([0.5,0.5],[0.5,10.5],[copula_globalKendalltau,copula_globalKendalltau],'--k','LineWidth', 2, 'Color', [0.85,0.33,0.10]);
        plot3([0.5,10.5],[10.5,10.5],[copula_globalKendalltau,copula_globalKendalltau],'--k','LineWidth', 2, 'Color', [0.85,0.33,0.10]);
        plot3([10.5,10.5],[0.5,10.5],[copula_globalKendalltau,copula_globalKendalltau],'--k','LineWidth', 2, 'Color', [0.85,0.33,0.10]);
    
        %% lines: for gumbel copula model
        plot3([9.5,10.5],[9.5,9.5],[copula_globalKendalltau,copula_globalKendalltau],'-','LineWidth', 0.5, 'Color', [0,0,0]);
        plot3([9.5,9.5],[9.5,10.5],[copula_globalKendalltau,copula_globalKendalltau],'-','LineWidth', 0.5, 'Color', [0,0,0]);
     
        % lines: for Clayton copula model
        % plot3([0.5,1.5],[0.5,0.5],[copula_globalKendalltau,copula_globalKendalltau],'-','LineWidth', 0.5, 'Color', [0,0,0]);
        % plot3([0.5,0.5],[0.5,1.5],[copula_globalKendalltau,copula_globalKendalltau],'-','LineWidth', 0.5, 'Color', [0,0,0]);

    case 0.2
        x = [0.5, 5.5, 5.5, 0.5];
        y = [0.5, 0.5, 5.5, 5.5];
        z = [copula_globalKendalltau, copula_globalKendalltau, copula_globalKendalltau, copula_globalKendalltau];
        fill3(x, y, z, [0.6, 0.6, 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

        plot3([0.5,5.5],[0.5,0.5],[copula_globalKendalltau,copula_globalKendalltau],'--k','LineWidth', 2, 'Color', [0.85,0.33,0.10]);
        plot3([0.5,0.5],[0.5,5.5],[copula_globalKendalltau,copula_globalKendalltau],'--k','LineWidth', 2, 'Color', [0.85,0.33,0.10]);
        plot3([0.5,5.5],[5.5,5.5],[copula_globalKendalltau,copula_globalKendalltau],'--k','LineWidth', 2, 'Color', [0.85,0.33,0.10]);
        plot3([5.5,5.5],[0.5,5.5],[copula_globalKendalltau,copula_globalKendalltau],'--k','LineWidth', 2, 'Color', [0.85,0.33,0.10]);

        %% lines: for gumbel copula model
        plot3([4.5,5.5],[4.5,4.5],[copula_globalKendalltau,copula_globalKendalltau],'-','LineWidth', 0.5, 'Color', [0,0,0]);
        plot3([4.5,4.5],[4.5,5.5],[copula_globalKendalltau,copula_globalKendalltau],'-','LineWidth', 0.5, 'Color', [0,0,0]);
     
        % lines: for Clayton copula model
        % plot3([0.5,1.5],[0.5,0.5],[copula_globalKendalltau,copula_globalKendalltau],'-','LineWidth', 0.5, 'Color', [0,0,0]);
        % plot3([0.5,0.5],[0.5,1.5],[copula_globalKendalltau,copula_globalKendalltau],'-','LineWidth', 0.5, 'Color', [0,0,0]);

    case 0.25
        x = [0.5, 4.5, 4.5, 0.5];
        y = [0.5, 0.5, 4.5, 4.5];
        z = [copula_globalKendalltau, copula_globalKendalltau, copula_globalKendalltau, copula_globalKendalltau];
        fill3(x, y, z, [0.6, 0.6, 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

        plot3([0.5,4.5],[0.5,0.5],[copula_globalKendalltau,copula_globalKendalltau],'--k','LineWidth', 2, 'Color', [0.85,0.33,0.10]);
        plot3([0.5,0.5],[0.5,4.5],[copula_globalKendalltau,copula_globalKendalltau],'--k','LineWidth', 2, 'Color', [0.85,0.33,0.10]);
        plot3([0.5,4.5],[4.5,4.5],[copula_globalKendalltau,copula_globalKendalltau],'--k','LineWidth', 2, 'Color', [0.85,0.33,0.10]);
        plot3([4.5,4.5],[0.5,4.5],[copula_globalKendalltau,copula_globalKendalltau],'--k','LineWidth', 2, 'Color', [0.85,0.33,0.10]);    
        
        %% lines: for gumbel copula model
        plot3([3.5,4.5],[3.5,3.5],[copula_globalKendalltau,copula_globalKendalltau],'-','LineWidth', 0.5, 'Color', [0,0,0]);
        plot3([3.5,3.5],[3.5,4.5],[copula_globalKendalltau,copula_globalKendalltau],'-','LineWidth', 0.5, 'Color', [0,0,0]);
        
        % lines: for Clayton copula model
        % plot3([0.5,1.5],[0.5,0.5],[copula_globalKendalltau,copula_globalKendalltau],'-','LineWidth', 0.5, 'Color', [0,0,0]);
        % plot3([0.5,0.5],[0.5,1.5],[copula_globalKendalltau,copula_globalKendalltau],'-','LineWidth', 0.5, 'Color', [0,0,0]);
   
    case 0.5
        x = [0.5, 2.5, 2.5, 0.5];
        y = [0.5, 0.5, 2.5, 2.5];
        z = [copula_globalKendalltau, copula_globalKendalltau, copula_globalKendalltau, copula_globalKendalltau];
        fill3(x, y, z, [0.6, 0.6, 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

        plot3([0.5,2.5],[0.5,0.5],[copula_globalKendalltau,copula_globalKendalltau],'--k','LineWidth', 2, 'Color', [0.85,0.33,0.10]);
        plot3([0.5,0.5],[0.5,2.5],[copula_globalKendalltau,copula_globalKendalltau],'--k','LineWidth', 2, 'Color', [0.85,0.33,0.10]);
        plot3([0.5,2.5],[2.5,2.5],[copula_globalKendalltau,copula_globalKendalltau],'--k','LineWidth', 2, 'Color', [0.85,0.33,0.10]);
        plot3([2.5,2.5],[0.5,2.5],[copula_globalKendalltau,copula_globalKendalltau],'--k','LineWidth', 2, 'Color', [0.85,0.33,0.10]);    
        
        %% lines: for gumbel copula model
        % plot3([1.5,2.5],[1.5,1.5],[copula_globalKendalltau,copula_globalKendalltau],'-','LineWidth', 0.5, 'Color', [0,0,0]);
        % plot3([1.5,1.5],[1.5,2.5],[copula_globalKendalltau,copula_globalKendalltau],'-','LineWidth', 0.5, 'Color', [0,0,0]);
        
        % lines: for Clayton copula model
        % plot3([0.5,1.5],[0.5,0.5],[copula_globalKendalltau,copula_globalKendalltau],'-','LineWidth', 0.5, 'Color', [0,0,0]);
        % plot3([0.5,0.5],[0.5,1.5],[copula_globalKendalltau,copula_globalKendalltau],'-','LineWidth', 0.5, 'Color', [0,0,0]);
end


