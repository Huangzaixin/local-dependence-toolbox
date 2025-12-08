%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%        Draw Figures for Simulation Results           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       Example 1: taking region 1 as an example       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Region 1: Empirical Standard Errors (SE_e)
% Read SE data from "simulation_results.xlsx"
region1_SE_u_statistic_par1 = xlsread('simulation_results.xlsx','Sheet1','E1:E7');
region1_SE_copula_model_par1 = xlsread('simulation_results.xlsx','Sheet1','G1:G7');

region1_SE_u_statistic_par2 = xlsread('simulation_results.xlsx','Sheet1','E8:E12');
region1_SE_copula_model_par2 = xlsread('simulation_results.xlsx','Sheet1','G8:G12');

region1_SE_u_statistic_par3 = xlsread('simulation_results.xlsx','Sheet1','E13:E17');
region1_SE_copula_model_par3 = xlsread('simulation_results.xlsx','Sheet1','G13:G17');

% Clayton copula: theta = 1; U-statistic-based estimator
plot(region1_SE_u_statistic_par1,'--o','Color',[0.00,0.45,0.74],'LineWidth',1,'MarkerSize',7.2); hold on; % [0.05,0.40,0.63]
% Clayton copula: theta = 1; Copula model-based estimator
plot(region1_SE_copula_model_par1,'-o','Color',[0.85,0.40,0.29],'LineWidth',2.5,'MarkerSize',7.2); hold on; % [0.78,0.31,0.19]

% Clayton copula: theta = 2; U-statistic-based estimator
plot(region1_SE_u_statistic_par2,'--x','Color',[0.00,0.45,0.74],'LineWidth',1,'MarkerSize',7.5); hold on;
% Clayton copula: theta = 2; Copula model-based estimator
plot(region1_SE_copula_model_par2,'-x','Color',[0.85,0.40,0.29],'LineWidth',2.5,'MarkerSize',7.5); hold on;

% Clayton copula: theta = 4; U-statistic-based estimator
plot(region1_SE_u_statistic_par3,'--^','Color',[0.00,0.45,0.74],'LineWidth',1,'MarkerSize',6.7); hold on;
% Clayton copula: theta = 4; Copula model-based estimator
plot(region1_SE_copula_model_par3,'-^','Color',[0.85,0.40,0.29],'LineWidth',2.5,'MarkerSize',6.7); hold on;

plot([1 5],[95 95],'-.','Color',[0,0,0],'LineWidth',1.5);

ylim([0 0.5])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5]);
set(gca,'Xticklabel',[100 200 500 1000 2000]);
set(gcf,'unit','centimeters','position',[15 15 9 10]);
set(gca,'Position', [.125 .17 .80 .74]);
set(gca,'lineWidth', 1);
legend('U-statistic, Clayton(u,v,1)','Copula model, Clayton(u,v,1)','U-statistic, Clayton(u,v,2)','Copula model, Clayton(u,v,2)','U-statistic, Clayton(u,v,4)','Copula model, Clayton(u,v,4)','fontsize',11.5,'Location','northwest','EdgeColor',[1 1 1])

xlabel('Region I   SE_{e} (Clayton copula)','fontsize',13.5);
set(gca,'FontSize',13.5);
set(gcf, 'color', [1 1 1]);


%% Region 1: 95% Coverage Probabilities (CP%)
% read CP data from "simulation_results.xlsx"
region1_CP_u_statistic_par1 = xlsread('simulation_results.xlsx','Sheet1','F1:F7');
region1_CP_copula_model_par1 = xlsread('simulation_results.xlsx','Sheet1','H1:H7');

region1_CP_u_statistic_par2 = xlsread('simulation_results.xlsx','Sheet1','F8:F12');
region1_CP_copula_model_par2 = xlsread('simulation_results.xlsx','Sheet1','H8:H12');

region1_CP_u_statistic_par3 = xlsread('simulation_results.xlsx','Sheet1','F13:F17');
region1_CP_copula_model_par3 = xlsread('simulation_results.xlsx','Sheet1','H13:H17');

plot(region1_CP_u_statistic_par1,'--o','Color',[0.00,0.45,0.74],'LineWidth',1,'MarkerSize',7.2); hold on;
plot(region1_CP_copula_model_par1,'-o','Color',[0.85,0.40,0.29],'LineWidth',2.5,'MarkerSize',7.2); hold on;

plot(region1_CP_u_statistic_par2,'--x','Color',[0.00,0.45,0.74],'LineWidth',1,'MarkerSize',7.5); hold on;
plot(region1_CP_copula_model_par2,'-x','Color',[0.85,0.40,0.29],'LineWidth',2.5,'MarkerSize',7.5); hold on;

plot(region1_CP_u_statistic_par3,'--^','Color',[0.00,0.45,0.74],'LineWidth',1,'MarkerSize',6.7); hold on;
plot(region1_CP_copula_model_par3,'-^','Color',[0.85,0.40,0.29],'LineWidth',2.5,'MarkerSize',6.7); hold on;

plot([1 5],[95 95],'-.','Color',[0,0,0],'LineWidth',2.5);

ylim([0 100]);
set(gca,'YTick',[0 20 40 60 80 100]);
set(gca,'Xticklabel',[100 200 500 1000 2000]);
set(gcf,'unit','centimeters','position',[15 15 9 10]);
set(gca,'Position', [.125 .17 .80 .74]);
set(gca,'lineWidth', 1);
legend('U-statistic, Clayton(u,v,1)','Copula model, Clayton(u,v,1)','U-statistic, Clayton(u,v,2)','Copula model, Clayton(u,v,2)','U-statistic, Clayton(u,v,4)','Copula model, Clayton(u,v,4)','CP=95%','fontsize',11.5,'Location','southwest','EdgeColor',[1 1 1]);

xlabel('Region I   CP% (Clayton copula)','fontsize',13.5);
set(gca,'FontSize',13.5);
set(gcf, 'color', [1 1 1]);



