function [ld] = fun_copulald_general(copulatype,weight1,weight2,copulaparameter1,copulaparameter2,copulaparameter3,measuretype,pl,pu,ql,qu)
% Description: calculate the generalized local Kendall's tau of a copula model
% Inputs: 
%      1. copulatype  
%         'gumbel'         : Gumbel copula 
%         'rotgumbel180'   : rotated Gumbel copula£¨180-degrees£©
%         'clayton'        : Clayton copula
%         'rotclayton180'  : rotated Clayton copula£¨180-degrees£©
%         'frank'          : Frank copula
%         'fgm'            : FGM copula
%         'sjc'            : SJC copula
%         'mixgc'          : Gumbel copula + Clayton copula
%         'mixrg180rc180'  : rotated Gumbel copula (180-degrees) + rotated Clayton copula£¨180-degrees£©
%         'mixgrg180'      : Gumbel copula + rotated Gumbel copula£¨180-degrees£©
%         'mixcrc180'      : Clayton copula + rotated Clayton copula£¨180-degrees£©
%         'mixfc'          : Frank copula + Clayton copula
%         'mixfrc180'      : Frank copula + rotated Clayton copula£¨180-degrees£©
%         'mixfg'          : Frank copula + Gumbel copula
%         'mixfrg180'      : Frank copula + rotated Gumbel copula£¨180-degrees£©
%         'mixfgrg180'     : Frank copula + Gumbel copula + rotated Gumbel copula£¨180-degrees£©
%      2. weight1          : the weight of the first copula function in the mixture copula model
%         weight2          : the weight of the second copula function in the mixture copula model
%      3. copulaparameter1 : for mixture copula model, parameter of the first copula function; for sjc copula, the upper tail dependence coefficient 
%         copulaparameter2 : for mixture copula model, parameter of the second copula function; for sjc copula, the lower tail dependence coefficient
%         copulaparameter3 : for mixture copula model, parameter of the third copula function
%      4. measuretype: 
%         'kendall' : local Kendall's tau
%      5. pl: lower quantile of X
%         pu: upper quantile of X
%         ql: lower quantile of Y
%         qu: upper quantile of Y
% Output: ld
%      1. ld: Type I local Kendall's tau of the selected copula model
% Author: Zaixin Huang
% Date: finished at 2023.01.01; this version: 2025.03.16
% Bug reports and suggestions: 
%       if you find any bugs or have suggestions, please contact me at eric.huangzaixin@gmail.com. 
%       I will update them on GitHub and acknowledge your contribution. Thank you!
% The latest version can be downloaded from https://github.com/huangzaixin/local-dependence-toolbox
%%
switch lower(measuretype) 
    
case 'kendall'
%% 2. local Kendall's tau
    switch lower(copulatype)
    
    case 'gumbel'   % Gumbel copula
          alpha = copulaparameter1; 
          Copula = @(u,v) exp(-((-log(u)).^alpha + (-log(v)).^alpha).^(alpha.^(-1)));
          % Second partial derivative of Gumbel copula
          diffCopula_uv = @(u,v) (exp(-((-log(u)).^alpha + (-log(v)).^alpha).^(1./alpha)).*(-log(u)).^(alpha - 1).*(-log(v)).^(alpha - 1).*((-log(u)).^alpha + (-log(v)).^alpha).^(2./alpha - 2))./(u.*v) - (alpha.*exp(-((-log(u)).^alpha + (-log(v)).^alpha).^(1./alpha)).*(-log(u)).^(alpha - 1).*(-log(v)).^(alpha - 1).*(1./alpha - 1).*((-log(u)).^alpha + (-log(v)).^alpha).^(1./alpha - 2))./(u.*v);
    
    case 'rotgumbel180'   % rotated Gumbel copula£¨180-degrees£©
          alpha = copulaparameter1;
          Copula = @(u,v) u + v - 1 + exp(-((-log(1-u)).^alpha + (-log(1-v)).^alpha).^(alpha.^(-1)));
          % Second partial derivative of rotated Gumbel copula
          diffCopula_uv = @(u,v) (exp(-((-log(1 - u)).^alpha + (-log(1 - v)).^alpha).^(1./alpha)).*(-log(1 - u)).^(alpha - 1).*(-log(1 - v)).^(alpha - 1).*((-log(1 - u)).^alpha + (-log(1 - v)).^alpha).^(2./alpha - 2))./((u - 1).*(v - 1)) - (alpha.*exp(-((-log(1 - u)).^alpha + (-log(1 - v)).^alpha).^(1./alpha)).*(-log(1 - u)).^(alpha - 1).*(-log(1 - v)).^(alpha - 1).*(1./alpha - 1).*((-log(1 - u)).^alpha + (-log(1 - v)).^alpha).^(1./alpha - 2))./((u - 1).*(v - 1));
          
    case 'clayton'   % Clayton copula
          alpha = copulaparameter1;
          Copula = @(u,v) (u.^(-alpha) + v.^(-alpha) - 1).^(-1./alpha);
          % Second partial derivative of Clayton copula
          diffCopula_uv = @(u,v) (alpha.*(1./alpha + 1))./(u.^(alpha + 1).*v.^(alpha + 1).*(1./u.^alpha + 1./v.^alpha - 1).^(1./alpha + 2));
   
    case 'rotclayton180'   % rotated Clayton copula£¨180-degrees£©
          alpha = copulaparameter1;
          Copula = @(u,v) u + v - 1 + ((1-u).^(-alpha) + (1-v).^(-alpha) - 1).^(-1./alpha);
          % Second partial derivative of rotated Clayton copula
          diffCopula_uv = @(u,v) (alpha.*(1./alpha + 1))./((1 - u).^(alpha + 1).*(1 - v).^(alpha + 1).*(1./(1 - u).^alpha + 1./(1 - v).^alpha - 1).^(1./alpha + 2));
    
    case 'frank'   % Frank copula
          alpha = copulaparameter1;
          % Second partial derivative of Frank copula
          Copula = @(u,v) (-1./alpha).*log(1 + (exp(-alpha.*u)-1).*(exp(-alpha.*v)-1)./(exp(-alpha)-1));
          diffCopula_uv = @(u,v) (alpha.*exp(-alpha.*u).*exp(-alpha.*v).*(exp(-alpha.*u) - 1).*(exp(-alpha.*v) - 1))./((exp(-alpha) - 1).^2.*(((exp(-alpha.*u) - 1).*(exp(-alpha.*v) - 1))./(exp(-alpha) - 1) + 1).^2) - (alpha.*exp(-alpha.*u).*exp(-alpha.*v))./((exp(-alpha) - 1).*(((exp(-alpha.*u) - 1).*(exp(-alpha.*v) - 1))./(exp(-alpha) - 1) + 1));

    case 'fgm'   % FGM copula
          alpha = copulaparameter1;
          %Copula = @(u,v) (u.*v + alpha.*u.*v.*(1-u).*(1-v));
          % Second partial derivative of FGM copula
          %diffCopula_uv = @(u,v) alpha.*u.*(v - 1) + alpha.*v.*(u - 1) + alpha.*(u - 1).*(v - 1) + alpha.*u.*v + 1;
           
          ld = 2.*alpha.*(pu-pl).*(qu-ql)./(9.*(1+alpha.*(pl+pu-1).*(ql+qu-1)).^2);
          return;
   
    case 'sjc'   % SJC copula
          tauU = copulaparameter1;
          tauL = copulaparameter2;
          LOG2 = log(2);
          % Copula = @(u,v)  u./2 + v./2 - (1 - ((1 - u.^(log(2)./log(2 - tauL))).^(log(2)./log(tauU)) + (1 - v.^(log(2)./log(2 - tauL))).^(log(2)./log(tauU)) - 1).^(log(tauU)./log(2))).^(log(2 - tauL)./log(2))./2 - (1 - ((1 - (1 - u).^(log(2)./log(2 - tauU))).^(log(2)./log(tauL)) + (1 - (1 - v).^(log(2)./log(2 - tauU))).^(log(2)./log(tauL)) - 1).^(log(tauL)./log(2))).^(log(2 - tauU)./log(2))./2 + 1./2;
          Copula = @(u,v) u./2 + v./2 - (1 - ((1 - u.^(LOG2./log(2 - tauL))).^(LOG2./log(tauU)) + (1 - v.^(LOG2./log(2 - tauL))).^(LOG2./log(tauU)) - 1).^(log(tauU)./LOG2)).^(log(2 - tauL)./LOG2)./2 - (1 - ((1 - (1 - u).^(LOG2./log(2 - tauU))).^(LOG2./log(tauL)) + (1 - (1 - v).^(LOG2./log(2 - tauU))).^(LOG2./log(tauL)) - 1).^(log(tauL)./LOG2)).^(log(2 - tauU)./LOG2)./2 + 1./2;
          % the second derivative of SJC copula
          % diffCopula_uv = @(u,v) (log(2).^2.*(log(tauL)./log(2) - 1).*(1 - ((1 - (1 - u).^(log(2)./log(2 - tauU))).^(log(2)./log(tauL)) + (1 - (1 - v).^(log(2)./log(2 - tauU))).^(log(2)./log(tauL)) - 1).^(log(tauL)./log(2))).^(log(2 - tauU)./log(2) - 1).*(1 - (1 - u).^(log(2)./log(2 - tauU))).^(log(2)./log(tauL) - 1).*(1 - (1 - v).^(log(2)./log(2 - tauU))).^(log(2)./log(tauL) - 1).*(1 - u).^(log(2)./log(2 - tauU) - 1).*(1 - v).^(log(2)./log(2 - tauU) - 1).*((1 - (1 - u).^(log(2)./log(2 - tauU))).^(log(2)./log(tauL)) + (1 - (1 - v).^(log(2)./log(2 - tauU))).^(log(2)./log(tauL)) - 1).^(log(tauL)./log(2) - 2))./(2.*log(2 - tauU).*log(tauL)) - (log(2).*u.^(log(2)./log(2 - tauL) - 1).*v.^(log(2)./log(2 - tauL) - 1).*(log(2 - tauL)./log(2) - 1).*(1 - ((1 - u.^(log(2)./log(2 - tauL))).^(log(2)./log(tauU)) + (1 - v.^(log(2)./log(2 - tauL))).^(log(2)./log(tauU)) - 1).^(log(tauU)./log(2))).^(log(2 - tauL)./log(2) - 2).*(1 - u.^(log(2)./log(2 - tauL))).^(log(2)./log(tauU) - 1).*(1 - v.^(log(2)./log(2 - tauL))).^(log(2)./log(tauU) - 1).*((1 - u.^(log(2)./log(2 - tauL))).^(log(2)./log(tauU)) + (1 - v.^(log(2)./log(2 - tauL))).^(log(2)./log(tauU)) - 1).^((2.*log(tauU))./log(2) - 2))./(2.*log(2 - tauL)) - (log(2).*(log(2 - tauU)./log(2) - 1).*(1 - ((1 - (1 - u).^(log(2)./log(2 - tauU))).^(log(2)./log(tauL)) + (1 - (1 - v).^(log(2)./log(2 - tauU))).^(log(2)./log(tauL)) - 1).^(log(tauL)./log(2))).^(log(2 - tauU)./log(2) - 2).*(1 - (1 - u).^(log(2)./log(2 - tauU))).^(log(2)./log(tauL) - 1).*(1 - (1 - v).^(log(2)./log(2 - tauU))).^(log(2)./log(tauL) - 1).*(1 - u).^(log(2)./log(2 - tauU) - 1).*(1 - v).^(log(2)./log(2 - tauU) - 1).*((1 - (1 - u).^(log(2)./log(2 - tauU))).^(log(2)./log(tauL)) + (1 - (1 - v).^(log(2)./log(2 - tauU))).^(log(2)./log(tauL)) - 1).^((2*log(tauL))./log(2) - 2))./(2.*log(2 - tauU)) + (log(2).^2.*u.^(log(2)./log(2 - tauL) - 1).*v.^(log(2)./log(2 - tauL) - 1).*(1 - ((1 - u.^(log(2)./log(2 - tauL))).^(log(2)./log(tauU)) + (1 - v.^(log(2)./log(2 - tauL))).^(log(2)./log(tauU)) - 1).^(log(tauU)./log(2))).^(log(2 - tauL)./log(2) - 1).*(log(tauU)./log(2) - 1).*(1 - u.^(log(2)./log(2 - tauL))).^(log(2)./log(tauU) - 1).*(1 - v.^(log(2)./log(2 - tauL))).^(log(2)./log(tauU) - 1).*((1 - u.^(log(2)./log(2 - tauL))).^(log(2)./log(tauU)) + (1 - v.^(log(2)./log(2 - tauL))).^(log(2)./log(tauU)) - 1).^(log(tauU)./log(2) - 2))./(2.*log(2 - tauL).*log(tauU));
          diffCopula_uv = @(u,v) (LOG2.^2.*(log(tauL)./LOG2 - 1).*(1 - ((1 - (1 - u).^(LOG2./log(2 - tauU))).^(LOG2./log(tauL)) + (1 - (1 - v).^(LOG2./log(2 - tauU))).^(LOG2./log(tauL)) - 1).^(log(tauL)./LOG2)).^(log(2 - tauU)./LOG2 - 1).*(1 - (1 - u).^(LOG2./log(2 - tauU))).^(LOG2/log(tauL) - 1).*(1 - (1 - v).^(LOG2./log(2 - tauU))).^(LOG2./log(tauL) - 1).*(1 - u).^(LOG2./log(2 - tauU) - 1).*(1 - v).^(LOG2./log(2 - tauU) - 1).*((1 - (1 - u).^(LOG2./log(2 - tauU))).^(LOG2./log(tauL)) + (1 - (1 - v).^(LOG2./log(2 - tauU))).^(LOG2./log(tauL)) - 1).^(log(tauL)./LOG2 - 2))./(2.*log(2 - tauU).*log(tauL)) - (LOG2.*u.^(LOG2./log(2 - tauL) - 1).*v.^(LOG2./log(2 - tauL) - 1).*(log(2 - tauL)./LOG2 - 1).*(1 - ((1 - u.^(LOG2./log(2 - tauL))).^(LOG2./log(tauU)) + (1 - v.^(LOG2./log(2 - tauL))).^(LOG2./log(tauU)) - 1).^(log(tauU)./LOG2)).^(log(2 - tauL)./LOG2 - 2).*(1 - u.^(LOG2./log(2 - tauL))).^(LOG2./log(tauU) - 1).*(1 - v.^(LOG2./log(2 - tauL))).^(LOG2./log(tauU) - 1).*((1 - u.^(LOG2./log(2 - tauL))).^(LOG2./log(tauU)) + (1 - v.^(LOG2./log(2 - tauL))).^(LOG2./log(tauU)) - 1).^((2.*log(tauU))./LOG2 - 2))./(2.*log(2 - tauL)) - (LOG2.*(log(2 - tauU)./LOG2 - 1).*(1 - ((1 - (1 - u).^(LOG2./log(2 - tauU))).^(LOG2./log(tauL)) + (1 - (1 - v).^(LOG2./log(2 - tauU))).^(LOG2./log(tauL)) - 1).^(log(tauL)./LOG2)).^(log(2 - tauU)./LOG2 - 2).*(1 - (1 - u).^(LOG2./log(2 - tauU))).^(LOG2./log(tauL) - 1).*(1 - (1 - v).^(LOG2./log(2 - tauU))).^(LOG2./log(tauL) - 1).*(1 - u).^(LOG2./log(2 - tauU) - 1).*(1 - v).^(LOG2./log(2 - tauU) - 1).*((1 - (1 - u).^(LOG2./log(2 - tauU))).^(LOG2./log(tauL)) + (1 - (1 - v).^(LOG2./log(2 - tauU))).^(LOG2./log(tauL)) - 1).^((2.*log(tauL))./LOG2 - 2))./(2.*log(2 - tauU)) + (LOG2.^2.*u.^(LOG2./log(2 - tauL) - 1).*v.^(LOG2./log(2 - tauL) - 1).*(1 - ((1 - u.^(LOG2./log(2 - tauL))).^(LOG2./log(tauU)) + (1 - v.^(LOG2./log(2 - tauL))).^(LOG2./log(tauU)) - 1).^(log(tauU)./LOG2)).^(log(2 - tauL)./LOG2 - 1).*(log(tauU)./LOG2 - 1).*(1 - u.^(LOG2./log(2 - tauL))).^(LOG2./log(tauU) - 1).*(1 - v.^(LOG2./log(2 - tauL))).^(LOG2./log(tauU) - 1).*((1 - u.^(LOG2./log(2 - tauL))).^(LOG2./log(tauU)) + (1 - v.^(LOG2./log(2 - tauL))).^(LOG2./log(tauU)) - 1).^(log(tauU)./LOG2 - 2))./(2.*log(2 - tauL).*log(tauU));
   
    case 'mixgc'   % MixGC copula: Gumbel copula + Clayton copula
          w1 = weight1;                 % weight of Gumbel copula
          alpha1 = copulaparameter1;    % parameter of Gumbel copula
          alpha2 = copulaparameter2;    % parameter of Clayton copula
          Copula = @(u,v) w1.*(exp(-((-log(u)).^alpha1 + (-log(v)).^alpha1).^(alpha1.^(-1)))) + (1-w1).*((u.^(-alpha2) + v.^(-alpha2) - 1).^(-1./alpha2));
          diffCopula_uv = @(u,v) (w1*exp(-((-log(u)).^alpha1 + (-log(v)).^alpha1).^(1./alpha1)).*(-log(u)).^(alpha1 - 1).*(-log(v)).^(alpha1 - 1).*((-log(u)).^alpha1 + (-log(v)).^alpha1).^(2./alpha1 - 2))./(u.*v) - (alpha2.*(1./alpha2 + 1).*(w1 - 1))./(u.^(alpha2 + 1).*v.^(alpha2 + 1).*(1./u.^alpha2 + 1./v.^alpha2 - 1).^(1./alpha2 + 2)) - (alpha1.*w1.*exp(-((-log(u)).^alpha1 + (-log(v)).^alpha1).^(1./alpha1)).*(-log(u)).^(alpha1 - 1).*(-log(v)).^(alpha1 - 1).*(1./alpha1 - 1).*((-log(u)).^alpha1 + (-log(v)).^alpha1).^(1./alpha1 - 2))./(u.*v);
    
    case 'mixrg180rc180'   % MixrGrC copula: rotated Gumbel copula£¨180-degrees£©+ rotated Clayton copula£¨180-degrees£©
          w1 = weight1;                 % weight of rotated Gumbel copula£¨180-degrees£©
          alpha1 = copulaparameter1;    % parameter of rotated Gumbel copula£¨180-degree£©
          alpha2 = copulaparameter2;    % parameter of rotated Clayton copula (180-degree)
          Copula = @(u,v) w1.*(u + v - 1 + exp(-((-log(1-u)).^alpha1 + (-log(1-v)).^alpha1).^(alpha1.^(-1)))) + (1-w1).*(u + v - 1 + ((1-u).^(-alpha2) + (1-v).^(-alpha2) - 1).^(-1./alpha2));
          diffCopula_uv = @(u,v) w1.*((exp(-((-log(1 - u)).^alpha1 + (-log(1 - v)).^alpha1).^(1./alpha1)).*(-log(1 - u)).^(alpha1 - 1).*(-log(1 - v)).^(alpha1 - 1).*((-log(1 - u)).^alpha1 + (-log(1 - v)).^alpha1).^(2./alpha1 - 2))./((u - 1).*(v - 1)) - (alpha1.*exp(-((-log(1 - u)).^alpha1 + (-log(1 - v)).^alpha1).^(1./alpha1)).*(-log(1 - u)).^(alpha1 - 1).*(-log(1 - v)).^(alpha1 - 1).*(1./alpha1 - 1).*((-log(1 - u)).^alpha1 + (-log(1 - v)).^alpha1).^(1./alpha1 - 2))./((u - 1).*(v - 1))) - (alpha2.*(1./alpha2 + 1).*(w1 - 1))./((1 - u).^(alpha2 + 1).*(1 - v).^(alpha2 + 1).*(1./(1 - u).^alpha2 + 1./(1 - v).^alpha2 - 1).^(1./alpha2 + 2));
  
    case 'mixgrg180'   % MixGrG copula: Gumbel copula + rotated Gumbel copula£¨180-degrees£©
          w1 = weight1;                 % weight of Gumbel copula 
          alpha1 = copulaparameter1;    % parameter of Gumbel copula
          alpha2 = copulaparameter2;    % parameter of rotated Gumbel copula (180-degrees)
          Copula = @(u,v)  w1.*( exp(-((-log(u)).^alpha1 + (-log(v)).^alpha1).^(alpha1.^(-1))) ) + (1-w1).*( u + v - 1 + exp(-((-log(1-u)).^alpha2 + (-log(1-v)).^alpha2).^(alpha2.^(-1))));
          diffCopula_uv = @(u,v) (w1.*exp(-((-log(u)).^alpha1 + (-log(v)).^alpha1).^(1./alpha1)).*(-log(u)).^(alpha1 - 1).*(-log(v)).^(alpha1 - 1).*((-log(u)).^alpha1 + (-log(v)).^alpha1).^(2./alpha1 - 2))./(u.*v) - ((exp(-((-log(1 - u)).^alpha2 + (-log(1 - v)).^alpha2).^(1./alpha2)).*(-log(1 - u)).^(alpha2 - 1).*(-log(1 - v)).^(alpha2 - 1).*((-log(1 - u)).^alpha2 + (-log(1 - v)).^alpha2).^(2./alpha2 - 2))./((u - 1).*(v - 1)) - (alpha2.*exp(-((-log(1 - u)).^alpha2 + (-log(1 - v)).^alpha2).^(1./alpha2)).*(-log(1 - u)).^(alpha2 - 1).*(-log(1 - v)).^(alpha2 - 1).*(1./alpha2 - 1).*((-log(1 - u)).^alpha2 + (-log(1 - v)).^alpha2).^(1./alpha2 - 2))./((u - 1).*(v - 1))).*(w1 - 1) - (alpha1.*w1.*exp(-((-log(u)).^alpha1 + (-log(v)).^alpha1).^(1./alpha1)).*(-log(u)).^(alpha1 - 1).*(-log(v)).^(alpha1 - 1).*(1./alpha1 - 1).*((-log(u)).^alpha1 + (-log(v)).^alpha1).^(1./alpha1 - 2))./(u.*v);

    case 'mixcrc180'   % MixCrC copula: Clayton copula + rotated Clayton copula£¨180-degrees£©
          w1 = weight1;                 % weight of Clayton copula
          alpha1 = copulaparameter1;    % parameter of Clayton copula
          alpha2 = copulaparameter2;    % parameter of rotated Clayton copula£¨180-degrees£©
          Copula = @(u,v) w1.*((u.^(-alpha1) + v.^(-alpha1) - 1).^(-1./alpha1)) + (1-w1).*(u + v - 1 + ((1-u).^(-alpha2) + (1-v).^(-alpha2) - 1).^(-1./alpha2));
          diffCopula_uv = @(u,v) (alpha1.*w1.*(1./alpha1 + 1))./(u.^(alpha1 + 1).*v.^(alpha1 + 1).*(1./u.^alpha1 + 1./v.^alpha1 - 1).^(1./alpha1 + 2)) - (alpha2.*(1./alpha2 + 1).*(w1 - 1))./((1 - u).^(alpha2 + 1).*(1 - v).^(alpha2 + 1).*(1./(1 - u).^alpha2 + 1./(1 - v).^alpha2 - 1).^(1./alpha2 + 2));

    case 'mixfc'   % MixFC copula: Frank copula + Clayton copula
          w1 = weight1;                 % weight of Frank copula
          alpha1 = copulaparameter1;    % parameter of Frank copula
          alpha2 = copulaparameter2;    % parameter of Clayton copula
          Copula = @(u,v) w1.*((-1./alpha1).*log(1 + (exp(-alpha1.*u)-1).*(exp(-alpha1.*v)-1)./(exp(-alpha1)-1))) + (1-w1).* ( (u.^(-alpha2) + v.^(-alpha2) - 1).^(-1./alpha2) );
          diffCopula_uv = @(u,v) (alpha1.*w1.*exp(-alpha1.*u).*exp(-alpha1.*v).*(exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./((exp(-alpha1) - 1).^2.*(((exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./(exp(-alpha1) - 1) + 1).^2) - (alpha1.*w1.*exp(-alpha1.*u).*exp(-alpha1.*v))./((exp(-alpha1) - 1).*(((exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./(exp(-alpha1) - 1) + 1)) - (alpha2.*(1./alpha2 + 1).*(w1 - 1))./(u.^(alpha2 + 1).*v.^(alpha2 + 1).*(1./u.^alpha2 + 1./v.^alpha2 - 1).^(1./alpha2 + 2));
   
    case 'mixfrc180'   % MixFrC180 copula: Frank copula + rotated Clayton copula£¨180-degrees£©
          w1 = weight1;                 % weight of Frank copula
          alpha1 = copulaparameter1;    % parameter of Frank copula
          alpha2 = copulaparameter2;    % parameter of rotated Clayton copula£¨180-degrees£© 
          Copula = @(u,v) w1.*((-1./alpha1).*log(1 + (exp(-alpha1.*u)-1).*(exp(-alpha1.*v)-1)./(exp(-alpha1)-1))) + (1-w1).* ( u + v - 1 + ((1-u).^(-alpha2) + (1-v).^(-alpha2) - 1).^(-1./alpha2) );
          diffCopula_uv = @(u,v) (alpha1.*w1.*exp(-alpha1.*u).*exp(-alpha1.*v).*(exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./((exp(-alpha1) - 1).^2.*(((exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./(exp(-alpha1) - 1) + 1).^2) - (alpha1.*w1.*exp(-alpha1.*u).*exp(-alpha1.*v))./((exp(-alpha1) - 1).*(((exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./(exp(-alpha1) - 1) + 1)) - (alpha2.*(1./alpha2 + 1).*(w1 - 1))./((1 - u).^(alpha2 + 1).*(1 - v).^(alpha2 + 1).*(1./(1 - u).^alpha2 + 1./(1 - v).^alpha2 - 1).^(1./alpha2 + 2));
          
    case 'mixfg'   % MixFG copula: Frank copula + Gumbel copula
          w1 = weight1;                 % weight of Frank copula
          alpha1 = copulaparameter1;    % parameter of Frank copula
          alpha2 = copulaparameter2;    % parameter of Gumbel copula 
          Copula = @(u,v) w1.*((-1./alpha1).*log(1 + (exp(-alpha1.*u)-1).*(exp(-alpha1.*v)-1)./(exp(-alpha1)-1))) + (1-w1).* (exp(-((-log(u)).^alpha2 + (-log(v)).^alpha2).^(alpha2.^(-1))));
          diffCopula_uv = @(u,v) (alpha1.*w1.*exp(-alpha1.*u).*exp(-alpha1.*v).*(exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./((exp(-alpha1) - 1).^2.*(((exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./(exp(-alpha1) - 1) + 1).^2) - (exp(-((-log(u)).^alpha2 + (-log(v)).^alpha2).^(1./alpha2)).*(-log(u)).^(alpha2 - 1).*(-log(v)).^(alpha2 - 1).*((-log(u)).^alpha2 + (-log(v)).^alpha2).^(2./alpha2 - 2).*(w1 - 1))./(u.*v) - (alpha1.*w1.*exp(-alpha1.*u).*exp(-alpha1.*v))./((exp(-alpha1) - 1).*(((exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./(exp(-alpha1) - 1) + 1)) + (alpha2.*exp(-((-log(u)).^alpha2 + (-log(v)).^alpha2).^(1./alpha2)).*(-log(u)).^(alpha2 - 1).*(-log(v)).^(alpha2 - 1).*(1./alpha2 - 1).*((-log(u)).^alpha2 + (-log(v)).^alpha2).^(1./alpha2 - 2).*(w1 - 1))./(u.*v);
  
    case 'mixfrg180'   % MixFrG180 copula: Frank copula + rotated Gumbel copula£¨180-degrees£©
          w1 = weight1;                 % weight of Frank copula
          alpha1 = copulaparameter1;    % parameter of Frank copula
          alpha2 = copulaparameter2;    % parameter of rotated Gumbel copula (180-degree)
          Copula = @(u,v) w1.*((-1./alpha1).*log(1 + (exp(-alpha1.*u)-1).*(exp(-alpha1.*v)-1)./(exp(-alpha1)-1))) + (1-w1).* ( u + v - 1 + exp(-((-log(1-u)).^alpha2 + (-log(1-v)).^alpha2).^(alpha2.^(-1))) );
          diffCopula_uv = @(u,v) (alpha1.*w1.*exp(-alpha1.*u).*exp(-alpha1.*v).*(exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./((exp(-alpha1) - 1).^2.*(((exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./(exp(-alpha1) - 1) + 1).^2) - (alpha1.*w1.*exp(-alpha1.*u).*exp(-alpha1.*v))./((exp(-alpha1) - 1).*(((exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./(exp(-alpha1) - 1) + 1)) - ((exp(-((-log(1 - u)).^alpha2 + (-log(1 - v)).^alpha2).^(1./alpha2)).*(-log(1 - u)).^(alpha2 - 1).*(-log(1 - v)).^(alpha2 - 1).*((-log(1 - u)).^alpha2 + (-log(1 - v)).^alpha2).^(2./alpha2 - 2))./((u - 1).*(v - 1)) - (alpha2.*exp(-((-log(1 - u)).^alpha2 + (-log(1 - v)).^alpha2).^(1./alpha2)).*(-log(1 - u)).^(alpha2 - 1).*(-log(1 - v)).^(alpha2 - 1).*(1./alpha2 - 1).*((-log(1 - u)).^alpha2 + (-log(1 - v)).^alpha2).^(1./alpha2 - 2))./((u - 1).*(v - 1))).*(w1 - 1);
     
    case 'mixfgrg180'   % MixFGrG180 copula: Frank copula + Gumbel copula + rotated Gumbel copula£¨180-degrees£©
          w1 = weight1;                 % weight of Frank copula
          w2 = weight2;                 % weight of Gumbel copula
          alpha1 = copulaparameter1;    % parameter of Frank copula
          alpha2 = copulaparameter2;    % parameter of Gumbel copula
          alpha3 = copulaparameter3;    % parameter of rotated Gumbel copula (180-degrees)
          Copula = @(u,v) w1.*((-1./alpha1).*log(1 + (exp(-alpha1.*u)-1).*(exp(-alpha1.*v)-1)./(exp(-alpha1)-1)))  +  w2.*( exp(-((-log(u)).^alpha2 + (-log(v)).^alpha2).^(alpha2.^(-1))) ) + (1-w1-w2).*( u + v - 1 + exp(-((-log(1-u)).^alpha3 + (-log(1-v)).^alpha3).^(alpha3.^(-1))));
          diffCopula_uv = @(u,v) (w2.*exp(-((-log(u)).^alpha2 + (-log(v)).^alpha2).^(1./alpha2)).*(-log(u)).^(alpha2 - 1).*(-log(v)).^(alpha2 - 1).*((-log(u)).^alpha2 + (-log(v)).^alpha2).^(2./alpha2 - 2))./(u.*v) - (alpha1.*w1.*exp(-alpha1.*u).*exp(-alpha1.*v))./((exp(-alpha1) - 1).*(((exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./(exp(-alpha1) - 1) + 1)) - ((exp(-((-log(1 - u)).^alpha3 + (-log(1 - v)).^alpha3).^(1./alpha3)).*(-log(1 - u)).^(alpha3 - 1).*(-log(1 - v)).^(alpha3 - 1).*((-log(1 - u)).^alpha3 + (-log(1 - v)).^alpha3).^(2./alpha3 - 2))./((u - 1).*(v - 1)) - (alpha3.*exp(-((-log(1 - u)).^alpha3 + (-log(1 - v)).^alpha3).^(1./alpha3)).*(-log(1 - u)).^(alpha3 - 1).*(-log(1 - v)).^(alpha3 - 1).*(1./alpha3 - 1).*((-log(1 - u)).^alpha3 + (-log(1 - v)).^alpha3).^(1./alpha3 - 2))./((u - 1).*(v - 1))).*(w1 + w2 - 1) + (alpha1.*w1.*exp(-alpha1.*u).*exp(-alpha1.*v).*(exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./((exp(-alpha1) - 1).^2.*(((exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./(exp(-alpha1) - 1) + 1).^2) - (alpha2.*w2.*exp(-((-log(u)).^alpha2 + (-log(v)).^alpha2).^(1./alpha2)).*(-log(u)).^(alpha2 - 1).*(-log(v)).^(alpha2 - 1).*(1./alpha2 - 1).*((-log(u)).^alpha2 + (-log(v)).^alpha2).^(1./alpha2 - 2))./(u.*v);
    
    % Here you can add other copula models
    
    otherwise
          error(message('ERROR: unsupported copula model.'));
    end
    
    %% copula model-based expression of the generalized local Kendall's tau
    numerator_internal_term = @(u,v) (2 .* Copula(u,v) - Copula(u,ql) - Copula(u,qu) - Copula(pl,v) - Copula(pu,v)) .* diffCopula_uv(u,v);
    % numerator_term = 2 .* (integral2(numerator_internal_term,pl,pu,ql,qu,'RelTol',1e-16));  
    numerator_term = 2 .* (quad2d(numerator_internal_term,pl,pu,ql,qu,'AbsTol', 1e-16, 'RelTol', 2e-10));  
    denominator_term = (Copula(pu,qu) - Copula(pl,qu) - Copula(pu,ql) + Copula(pl,ql)).^2;
    constant_term = (Copula(pu,qu) + Copula(pl,qu) + Copula(pu,ql) + Copula(pl,ql))./(Copula(pu,qu) - Copula(pl,qu) - Copula(pu,ql) + Copula(pl,ql));
    ld = numerator_term ./ denominator_term + constant_term;
    
    disp("Calculating local Kendall's tau using copula model... output: " + ld);
    
otherwise
    error(message('ERROR: unknown measure type. The current toolbox is only support Kendall''s tau.'));
end


