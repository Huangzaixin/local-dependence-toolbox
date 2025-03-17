%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Example 3: Estimation of some single and Archimedean mixture copula models  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% basic setups
break;
format long g;
options = optimset('Display','iter','TolCon',10^-9,'TolFun',10^-8,'TolX',10^-8);

%% read data
u = csvread('data examples/example_3/data/u.csv',1,0);
v = csvread('data examples/example_3/data/v.csv',1,0);

%% Gaussion copula
par_gauss = corrcoef12(norminv(u),norminv(v));
LL_gauss = NormalCopula_CL(par_gauss,[u,v]);
AIC_gauss = -2*(-LL_gauss) + 2;
AIC_gauss   % -176.75

%% Student's t copula
lower = [-0.999, 2.01];
upper = [ 0.999, 100];
par_0 = [par_gauss;10];
[para_t LL_t] = fmincon('tcopulaCL',par_0,[],[],[],[],lower,upper,[],options,[u,v]);
AIC = -2*(-LL_t) + 2*2;
AIC   %  -212.50
para_t

%% Clayton copula
lower = 0.001;    % lower bound 
par_0 = 1.001;    % starting values
[par_clayton,LL] = fmincon('claytonCL',par_0,[],[],[],[],lower,[],[],options,[u,v]);
AIC = -2*(-LL) + 2;
AIC    %  -104.37
par_clayton

%% Rotated Clayton copula (180-degrees)
lower = 0.001;
par_0 = 1.001;
[par_rotclayton LL] = fmincon('claytonCL',par_0,[],[],[],[],lower,[],[],options,1-[u,v]);
AIC = -2*(-LL) + 2;
AIC    %  -192.60
par_rotclayton

%% Gumbel copula
lower = 1.001;
par_0 = 2;
[par_gumbel LL] = fmincon('gumbelCL',par_0,[],[],[],[],lower,[],[],options,[u,v]);
AIC = -2*(-LL) + 2;
AIC    % -217.49
par_gumbel

%% Rotated Gumbel copula (180-degrees)
lower = 1.001;
par_0 = 2;
[par_rotgumbel LL] = fmincon('gumbelCL',par_0,[],[],[],[],lower,[],[],options,1-[u,v]);
AIC = -2*(-LL) + 2;
AIC    % -152.19
par_rotgumbel

%% MixGC copula : Gumbel copula + Clayton copula
lower = [0.001 1.001 0.001];    % lower bound 
upper = [0.999 +Inf +Inf];      % upper bound 
par_0 = [0.5 2 1];              % starting values
[pars_gc LL] = fmincon('mixGCCL',par_0,[],[],[],[],lower,upper,[],options,[u,v]);
AIC = -2*(-LL) + 2*3;
AIC   % -229.55
pars_gc

%% MixrGrC copula : rotated Gumbel copula (180-degrees) + rotated Clayton copula (180-degrees)
lower = [0.001 1.001 0.001];
upper = [0.999 +Inf +Inf];
par_0 = [0.2 2 2];
[pars_rgrc LL] = fmincon('mixrG180rC180CL',par_0,[],[],[],[],lower,upper,[],options,[u,v]);
AIC = -2*(-LL) + 2*3;
AIC   % -221.12
pars_rgrc

%% MixGrG180 copula : Gumbel copula + rotated Gumbel copula (180-degrees)
lower = [0.001 1.001 1.001];
upper = [0.999 +Inf +Inf];
par_0 = [0.3 3 3];
[pars_grg LL] = fmincon('mixGrG180CL',par_0,[],[],[],[],lower,upper,[],options,[u,v]);
AIC = -2*(-LL) + 2*3;
AIC   % -223.17
pars_grg

%% MixCrC180 copula : Clayton copula + rotated Clayton copula (180-degrees)
lower = [0.001 1.001 1.001];
upper = [0.999 +Inf +Inf];
par_0 = [0.3 2 2];
[pars_crc LL] = fmincon('mixCrC180CL',par_0,[],[],[],[],lower,upper,[],options,[u,v]);
AIC = -2*(-LL) + 2*3;
AIC   % -201.66
pars_crc

%% MixFC copula : Frank copula + Clayton copula
lower = [0.001 -Inf 0.001];
upper = [0.999 +Inf +Inf];
par_0 = [0.5 3 3];
[pars_fc LL] = fmincon('mixFCCL',par_0,[],[],[],[],lower,upper,[],options,[u,v]);
AIC = -2*(-LL) + 2*3;
AIC   % -189.02
pars_fc

%% MixFrC copula : Frank copula + rotated Clayton copula (180-degrees)
lower = [0.001 -Inf 0.001];
upper = [0.999 +Inf +Inf];
par_0 = [0.5 3 2];
[pars_frc LL] = fmincon('mixFrCCL',par_0,[],[],[],[],lower,upper,[],options,[u,v]);
AIC = -2*(-LL) + 2*3;
AIC   % -220.63
pars_frc

%% MixFG copula : Frank copula + Gumbel copula
lower = [0.001 -Inf 1.001];
upper = [0.999 +Inf +Inf];
par_0 = [0.5 2 2];
[pars_fg LL] = fmincon('mixFGCL',par_0,[],[],[],[],lower,upper,[],options,[u,v]);
AIC = -2*(-LL) + 2*3;
AIC   % -230.26, the smallest AIC
pars_fg

%% MixFG copula : Frank copula + rotated Gumbel copula (180-degrees)
lower = [0.001 -Inf 1.01];
upper = [0.999 +Inf +Inf];
par_0 = [0.5 3 3];
[pars_frg LL] = fmincon('mixFrGCL',par_0,[],[],[],[],lower,upper,[],options,[u,v]);
AIC = -2*(-LL) + 2*3;
AIC   % -205.97
pars_frg



%% %%%%%%%%%%%%%%%%%%%%%  other copula models  %%%%%%%%%%%%%%%%%%%%%%%
% MixFGrG180 copula : Frank copula + Gumbel copula + rotated Gumbel copula (180-degrees)
lower = [0.01 0.01 -Inf 1.01 1.01];
upper = [0.99 0.99 +Inf +Inf +Inf];
par_0 = [0.3 0.3 3 3 3];
[pars_fgrg LL] = fmincon('mixFGrG180CL',par_0,[],[],[],[],lower,upper,[],options,[u,v]);
AIC = -2*(-LL) + 2*5;
AIC  % -227.51
pars_fgrg

% Patton's SJC copula
lower = 1e-44*ones(2,1);	
upper = 1-1e-4*ones(2,1);  
tauU = 0.3;   % upper
tauL = 0.6;   % lower
par_0 = [tauU;tauL];
[pars_sjc LL] = fmincon('sym_jc_CL',par_0,[],[],[],[],lower,upper,[],options,[u,v]);
AIC = -2 * -LL + 2*2;
AIC   % -205.21
pars_sjc;

% generating data
T = 10000;
data = sym_jc_rnd(tauU,tauL,T);	


