function [ld] = fun_copulald_type_II(copulatype,weight1,weight2,copulaparameter1,copulaparameter2,copulaparameter3,measuretype,region,Xquantile,Yquantile)
% Description: Type II local Kendall's tau of a copula model
% Note: Formulas (11)-(14) demonstrate higher accuracy compared to formula (6), primarily due to their reduced number of double integrals.
% Inputs:
%      1. copulatype
%         'gumbel'         : Gumbel copula 
%         'clayton'        : Clayton copula
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
%      2. weight1          : the first copula weight of the mixture copula
%         weight2          : the second copula weight of the mixture copula 
%      3. copulaparameter1 : parameter of the first copula
%         copulaparameter2 : parameter of the second copula
%         copulaparameter3 : parameter of the third copula
%      4. measuretype
%         'Kendall' : local Kendall's tau
%      5. region:
%         'UU': upper-upper region; 
%         'UL': upper-lower region;
%         'LU': lower-upper region; 
%         'LL': lower-lower region;
%      6. Xquantile: quantile of X
%         Yquantile: quantile of Y
% Output: ld
%      1. ld: Type II local Kendall's tau of the selected copula model
% Author: Zaixin Huang
% Date: finished at 2015.06.07; updated at 2018.01.04, 2023.01.01; this version: 2025.03.16
% Bug reports and suggestions: 
%       if you find any bugs or have suggestions, please contact me at eric.huangzaixin@gmail.com. 
%       I will update them on GitHub and acknowledge your contribution. Thank you!
% The latest version can be downloaded from https://github.com/huangzaixin/local-dependence-toolbox
%%
format long;
current_p = Xquantile;
current_q = Yquantile;

%%
switch lower(measuretype)
case 'kendall'
    switch lower(copulatype)
    case 'gumbel'   % Gumbel copula
          alpha = copulaparameter1; 
          Copula = @(u,v) exp(-((-log(u)).^alpha + (-log(v)).^alpha).^(alpha.^(-1)));
          % Second partial derivative of Gumbel copula
          diffCopula_uv = @(u,v) (exp(-((-log(u)).^alpha + (-log(v)).^alpha).^(1./alpha)).*(-log(u)).^(alpha - 1).*(-log(v)).^(alpha - 1).*((-log(u)).^alpha + (-log(v)).^alpha).^(2./alpha - 2))./(u.*v) - (alpha.*exp(-((-log(u)).^alpha + (-log(v)).^alpha).^(1./alpha)).*(-log(u)).^(alpha - 1).*(-log(v)).^(alpha - 1).*(1./alpha - 1).*((-log(u)).^alpha + (-log(v)).^alpha).^(1./alpha - 2))./(u.*v);
         
    case 'clayton'  % Clayton copula
          alpha = copulaparameter1;
          Copula = @(u,v) (u.^(-alpha) + v.^(-alpha) - 1).^(-1./alpha);
          % Second partial derivative of Clayton copula
          diffCopula_uv = @(u,v) (alpha.*(1./alpha + 1))./(u.^(alpha + 1).*v.^(alpha + 1).*(1./u.^alpha + 1./v.^alpha - 1).^(1./alpha + 2));
      
    case 'frank'    % Frank copula
          alpha = copulaparameter1;
          Copula = @(u,v) (-1./alpha).*log(1 + (exp(-alpha.*u)-1).*(exp(-alpha.*v)-1)./(exp(-alpha)-1));
          % Second partial derivative of Frank copula
          diffCopula_uv = @(u,v) (alpha.*exp(-alpha.*u).*exp(-alpha.*v).*(exp(-alpha.*u) - 1).*(exp(-alpha.*v) - 1))./((exp(-alpha) - 1).^2.*(((exp(-alpha.*u) - 1).*(exp(-alpha.*v) - 1))./(exp(-alpha) - 1) + 1).^2) - (alpha.*exp(-alpha.*u).*exp(-alpha.*v))./((exp(-alpha) - 1).*(((exp(-alpha.*u) - 1).*(exp(-alpha.*v) - 1))./(exp(-alpha) - 1) + 1));

    case 'fgm'      % FGM copula
          alpha = copulaparameter1;
          Copula = @(u,v) (u.*v + alpha.*u.*v.*(1-u).*(1-v));
          % Second partial derivative of FGM copula
          diffCopula_uv = @(u,v) alpha.*u.*(v - 1) + alpha.*v.*(u - 1) + alpha.*(u - 1).*(v - 1) + alpha.*u.*v + 1;
   
    case 'mixgc'    % MixGC copula : Gumbel copula + Clayton copula
          w1 = weight1;                 % weight of Gumbel copula
          alpha1 = copulaparameter1;    % parameter of Gumbel copula
          alpha2 = copulaparameter2;    % parameter of Clayton copula
          Copula = @(u,v) w1.*(exp(-((-log(u)).^alpha1 + (-log(v)).^alpha1).^(alpha1.^(-1)))) + (1-w1).*((u.^(-alpha2) + v.^(-alpha2) - 1).^(-1./alpha2));
          diffCopula_uv = @(u,v) (w1*exp(-((-log(u)).^alpha1 + (-log(v)).^alpha1).^(1./alpha1)).*(-log(u)).^(alpha1 - 1).*(-log(v)).^(alpha1 - 1).*((-log(u)).^alpha1 + (-log(v)).^alpha1).^(2./alpha1 - 2))./(u.*v) - (alpha2.*(1./alpha2 + 1).*(w1 - 1))./(u.^(alpha2 + 1).*v.^(alpha2 + 1).*(1./u.^alpha2 + 1./v.^alpha2 - 1).^(1./alpha2 + 2)) - (alpha1.*w1.*exp(-((-log(u)).^alpha1 + (-log(v)).^alpha1).^(1./alpha1)).*(-log(u)).^(alpha1 - 1).*(-log(v)).^(alpha1 - 1).*(1./alpha1 - 1).*((-log(u)).^alpha1 + (-log(v)).^alpha1).^(1./alpha1 - 2))./(u.*v);
   
    case 'sjc'     % SJC copula
          tauU = copulaparameter1;
          tauL = copulaparameter2;
          LOG2 = log(2);
          % Copula = @(u,v)  u./2 + v./2 - (1 - ((1 - u.^(log(2)./log(2 - tauL))).^(log(2)./log(tauU)) + (1 - v.^(log(2)./log(2 - tauL))).^(log(2)./log(tauU)) - 1).^(log(tauU)./log(2))).^(log(2 - tauL)./log(2))./2 - (1 - ((1 - (1 - u).^(log(2)./log(2 - tauU))).^(log(2)./log(tauL)) + (1 - (1 - v).^(log(2)./log(2 - tauU))).^(log(2)./log(tauL)) - 1).^(log(tauL)./log(2))).^(log(2 - tauU)./log(2))./2 + 1./2;
          Copula = @(u,v) u./2 + v./2 - (1 - ((1 - u.^(LOG2./log(2 - tauL))).^(LOG2./log(tauU)) + (1 - v.^(LOG2./log(2 - tauL))).^(LOG2./log(tauU)) - 1).^(log(tauU)./LOG2)).^(log(2 - tauL)./LOG2)./2 - (1 - ((1 - (1 - u).^(LOG2./log(2 - tauU))).^(LOG2./log(tauL)) + (1 - (1 - v).^(LOG2./log(2 - tauU))).^(LOG2./log(tauL)) - 1).^(log(tauL)./LOG2)).^(log(2 - tauU)./LOG2)./2 + 1./2;
          % the second derivative of SJC copula
          % diffCopula_uv = @(u,v) (log(2).^2.*(log(tauL)./log(2) - 1).*(1 - ((1 - (1 - u).^(log(2)./log(2 - tauU))).^(log(2)./log(tauL)) + (1 - (1 - v).^(log(2)./log(2 - tauU))).^(log(2)./log(tauL)) - 1).^(log(tauL)./log(2))).^(log(2 - tauU)./log(2) - 1).*(1 - (1 - u).^(log(2)./log(2 - tauU))).^(log(2)./log(tauL) - 1).*(1 - (1 - v).^(log(2)./log(2 - tauU))).^(log(2)./log(tauL) - 1).*(1 - u).^(log(2)./log(2 - tauU) - 1).*(1 - v).^(log(2)./log(2 - tauU) - 1).*((1 - (1 - u).^(log(2)./log(2 - tauU))).^(log(2)./log(tauL)) + (1 - (1 - v).^(log(2)./log(2 - tauU))).^(log(2)./log(tauL)) - 1).^(log(tauL)./log(2) - 2))./(2.*log(2 - tauU).*log(tauL)) - (log(2).*u.^(log(2)./log(2 - tauL) - 1).*v.^(log(2)./log(2 - tauL) - 1).*(log(2 - tauL)./log(2) - 1).*(1 - ((1 - u.^(log(2)./log(2 - tauL))).^(log(2)./log(tauU)) + (1 - v.^(log(2)./log(2 - tauL))).^(log(2)./log(tauU)) - 1).^(log(tauU)./log(2))).^(log(2 - tauL)./log(2) - 2).*(1 - u.^(log(2)./log(2 - tauL))).^(log(2)./log(tauU) - 1).*(1 - v.^(log(2)./log(2 - tauL))).^(log(2)./log(tauU) - 1).*((1 - u.^(log(2)./log(2 - tauL))).^(log(2)./log(tauU)) + (1 - v.^(log(2)./log(2 - tauL))).^(log(2)./log(tauU)) - 1).^((2.*log(tauU))./log(2) - 2))./(2.*log(2 - tauL)) - (log(2).*(log(2 - tauU)./log(2) - 1).*(1 - ((1 - (1 - u).^(log(2)./log(2 - tauU))).^(log(2)./log(tauL)) + (1 - (1 - v).^(log(2)./log(2 - tauU))).^(log(2)./log(tauL)) - 1).^(log(tauL)./log(2))).^(log(2 - tauU)./log(2) - 2).*(1 - (1 - u).^(log(2)./log(2 - tauU))).^(log(2)./log(tauL) - 1).*(1 - (1 - v).^(log(2)./log(2 - tauU))).^(log(2)./log(tauL) - 1).*(1 - u).^(log(2)./log(2 - tauU) - 1).*(1 - v).^(log(2)./log(2 - tauU) - 1).*((1 - (1 - u).^(log(2)./log(2 - tauU))).^(log(2)./log(tauL)) + (1 - (1 - v).^(log(2)./log(2 - tauU))).^(log(2)./log(tauL)) - 1).^((2*log(tauL))./log(2) - 2))./(2.*log(2 - tauU)) + (log(2).^2.*u.^(log(2)./log(2 - tauL) - 1).*v.^(log(2)./log(2 - tauL) - 1).*(1 - ((1 - u.^(log(2)./log(2 - tauL))).^(log(2)./log(tauU)) + (1 - v.^(log(2)./log(2 - tauL))).^(log(2)./log(tauU)) - 1).^(log(tauU)./log(2))).^(log(2 - tauL)./log(2) - 1).*(log(tauU)./log(2) - 1).*(1 - u.^(log(2)./log(2 - tauL))).^(log(2)./log(tauU) - 1).*(1 - v.^(log(2)./log(2 - tauL))).^(log(2)./log(tauU) - 1).*((1 - u.^(log(2)./log(2 - tauL))).^(log(2)./log(tauU)) + (1 - v.^(log(2)./log(2 - tauL))).^(log(2)./log(tauU)) - 1).^(log(tauU)./log(2) - 2))./(2.*log(2 - tauL).*log(tauU));
          diffCopula_uv = @(u,v) (LOG2.^2.*(log(tauL)./LOG2 - 1).*(1 - ((1 - (1 - u).^(LOG2./log(2 - tauU))).^(LOG2./log(tauL)) + (1 - (1 - v).^(LOG2./log(2 - tauU))).^(LOG2./log(tauL)) - 1).^(log(tauL)./LOG2)).^(log(2 - tauU)./LOG2 - 1).*(1 - (1 - u).^(LOG2./log(2 - tauU))).^(LOG2/log(tauL) - 1).*(1 - (1 - v).^(LOG2./log(2 - tauU))).^(LOG2./log(tauL) - 1).*(1 - u).^(LOG2./log(2 - tauU) - 1).*(1 - v).^(LOG2./log(2 - tauU) - 1).*((1 - (1 - u).^(LOG2./log(2 - tauU))).^(LOG2./log(tauL)) + (1 - (1 - v).^(LOG2./log(2 - tauU))).^(LOG2./log(tauL)) - 1).^(log(tauL)./LOG2 - 2))./(2.*log(2 - tauU).*log(tauL)) - (LOG2.*u.^(LOG2./log(2 - tauL) - 1).*v.^(LOG2./log(2 - tauL) - 1).*(log(2 - tauL)./LOG2 - 1).*(1 - ((1 - u.^(LOG2./log(2 - tauL))).^(LOG2./log(tauU)) + (1 - v.^(LOG2./log(2 - tauL))).^(LOG2./log(tauU)) - 1).^(log(tauU)./LOG2)).^(log(2 - tauL)./LOG2 - 2).*(1 - u.^(LOG2./log(2 - tauL))).^(LOG2./log(tauU) - 1).*(1 - v.^(LOG2./log(2 - tauL))).^(LOG2./log(tauU) - 1).*((1 - u.^(LOG2./log(2 - tauL))).^(LOG2./log(tauU)) + (1 - v.^(LOG2./log(2 - tauL))).^(LOG2./log(tauU)) - 1).^((2.*log(tauU))./LOG2 - 2))./(2.*log(2 - tauL)) - (LOG2.*(log(2 - tauU)./LOG2 - 1).*(1 - ((1 - (1 - u).^(LOG2./log(2 - tauU))).^(LOG2./log(tauL)) + (1 - (1 - v).^(LOG2./log(2 - tauU))).^(LOG2./log(tauL)) - 1).^(log(tauL)./LOG2)).^(log(2 - tauU)./LOG2 - 2).*(1 - (1 - u).^(LOG2./log(2 - tauU))).^(LOG2./log(tauL) - 1).*(1 - (1 - v).^(LOG2./log(2 - tauU))).^(LOG2./log(tauL) - 1).*(1 - u).^(LOG2./log(2 - tauU) - 1).*(1 - v).^(LOG2./log(2 - tauU) - 1).*((1 - (1 - u).^(LOG2./log(2 - tauU))).^(LOG2./log(tauL)) + (1 - (1 - v).^(LOG2./log(2 - tauU))).^(LOG2./log(tauL)) - 1).^((2.*log(tauL))./LOG2 - 2))./(2.*log(2 - tauU)) + (LOG2.^2.*u.^(LOG2./log(2 - tauL) - 1).*v.^(LOG2./log(2 - tauL) - 1).*(1 - ((1 - u.^(LOG2./log(2 - tauL))).^(LOG2./log(tauU)) + (1 - v.^(LOG2./log(2 - tauL))).^(LOG2./log(tauU)) - 1).^(log(tauU)./LOG2)).^(log(2 - tauL)./LOG2 - 1).*(log(tauU)./LOG2 - 1).*(1 - u.^(LOG2./log(2 - tauL))).^(LOG2./log(tauU) - 1).*(1 - v.^(LOG2./log(2 - tauL))).^(LOG2./log(tauU) - 1).*((1 - u.^(LOG2./log(2 - tauL))).^(LOG2./log(tauU)) + (1 - v.^(LOG2./log(2 - tauL))).^(LOG2./log(tauU)) - 1).^(log(tauU)./LOG2 - 2))./(2.*log(2 - tauL).*log(tauU));
  
    case 'mixrg180rc180'   % MixrGrC copula : rotated Gumbel copula£¨180-degrees£©+ rotated Clayton copula£¨180-degrees£©
          w1 = weight1;                 % weight of rotated Gumbel copula£¨180-degrees£©
          alpha1 = copulaparameter1;    % parameter of rotated Gumbel copula£¨180-degree£©
          alpha2 = copulaparameter2;    % parameter of rotated Clayton copula (180-degree)
          Copula = @(u,v) w1.*(u + v - 1 + exp(-((-log(1-u)).^alpha1 + (-log(1-v)).^alpha1).^(alpha1.^(-1)))) + (1-w1).*(u + v - 1 + ((1-u).^(-alpha2) + (1-v).^(-alpha2) - 1).^(-1./alpha2));
          diffCopula_uv = @(u,v) w1.*((exp(-((-log(1 - u)).^alpha1 + (-log(1 - v)).^alpha1).^(1./alpha1)).*(-log(1 - u)).^(alpha1 - 1).*(-log(1 - v)).^(alpha1 - 1).*((-log(1 - u)).^alpha1 + (-log(1 - v)).^alpha1).^(2./alpha1 - 2))./((u - 1).*(v - 1)) - (alpha1.*exp(-((-log(1 - u)).^alpha1 + (-log(1 - v)).^alpha1).^(1./alpha1)).*(-log(1 - u)).^(alpha1 - 1).*(-log(1 - v)).^(alpha1 - 1).*(1./alpha1 - 1).*((-log(1 - u)).^alpha1 + (-log(1 - v)).^alpha1).^(1./alpha1 - 2))./((u - 1).*(v - 1))) - (alpha2.*(1./alpha2 + 1).*(w1 - 1))./((1 - u).^(alpha2 + 1).*(1 - v).^(alpha2 + 1).*(1./(1 - u).^alpha2 + 1./(1 - v).^alpha2 - 1).^(1./alpha2 + 2));
  
    case 'mixgrg180'   % MixGrG copula : Gumbel copula + rotated Gumbel copula£¨180-degrees£©
          w1 = weight1;                 % weight of Gumbel copula 
          alpha1 = copulaparameter1;    % parameter of Gumbel copula
          alpha2 = copulaparameter2;    % parameter of rotated Gumbel copula (180-degrees)
          Copula = @(u,v)  w1.*( exp(-((-log(u)).^alpha1 + (-log(v)).^alpha1).^(alpha1.^(-1))) ) + (1-w1).*( u + v - 1 + exp(-((-log(1-u)).^alpha2 + (-log(1-v)).^alpha2).^(alpha2.^(-1))));
          diffCopula_uv = @(u,v) (w1.*exp(-((-log(u)).^alpha1 + (-log(v)).^alpha1).^(1./alpha1)).*(-log(u)).^(alpha1 - 1).*(-log(v)).^(alpha1 - 1).*((-log(u)).^alpha1 + (-log(v)).^alpha1).^(2./alpha1 - 2))./(u.*v) - ((exp(-((-log(1 - u)).^alpha2 + (-log(1 - v)).^alpha2).^(1./alpha2)).*(-log(1 - u)).^(alpha2 - 1).*(-log(1 - v)).^(alpha2 - 1).*((-log(1 - u)).^alpha2 + (-log(1 - v)).^alpha2).^(2./alpha2 - 2))./((u - 1).*(v - 1)) - (alpha2.*exp(-((-log(1 - u)).^alpha2 + (-log(1 - v)).^alpha2).^(1./alpha2)).*(-log(1 - u)).^(alpha2 - 1).*(-log(1 - v)).^(alpha2 - 1).*(1./alpha2 - 1).*((-log(1 - u)).^alpha2 + (-log(1 - v)).^alpha2).^(1./alpha2 - 2))./((u - 1).*(v - 1))).*(w1 - 1) - (alpha1.*w1.*exp(-((-log(u)).^alpha1 + (-log(v)).^alpha1).^(1./alpha1)).*(-log(u)).^(alpha1 - 1).*(-log(v)).^(alpha1 - 1).*(1./alpha1 - 1).*((-log(u)).^alpha1 + (-log(v)).^alpha1).^(1./alpha1 - 2))./(u.*v);

    case 'mixcrc180'   % MixCrC copula : Clayton copula + rotated Clayton copula£¨180-degrees£©
          w1 = weight1;                 % weight of Clayton copula
          alpha1 = copulaparameter1;    % parameter of Clayton copula
          alpha2 = copulaparameter2;    % parameter of rotated Clayton copula£¨180-degrees£©
          Copula = @(u,v) w1.*((u.^(-alpha1) + v.^(-alpha1) - 1).^(-1./alpha1)) + (1-w1).*(u + v - 1 + ((1-u).^(-alpha2) + (1-v).^(-alpha2) - 1).^(-1./alpha2));
          diffCopula_uv = @(u,v) (alpha1.*w1.*(1./alpha1 + 1))./(u.^(alpha1 + 1).*v.^(alpha1 + 1).*(1./u.^alpha1 + 1./v.^alpha1 - 1).^(1./alpha1 + 2)) - (alpha2.*(1./alpha2 + 1).*(w1 - 1))./((1 - u).^(alpha2 + 1).*(1 - v).^(alpha2 + 1).*(1./(1 - u).^alpha2 + 1./(1 - v).^alpha2 - 1).^(1./alpha2 + 2));

    case 'mixfc'   % MixFC copula : Frank copula + Clayton copula
          w1 = weight1;                 % weight of Frank copula
          alpha1 = copulaparameter1;    % parameter of Frank copula
          alpha2 = copulaparameter2;    % parameter of Clayton copula
          Copula = @(u,v) w1.*((-1./alpha1).*log(1 + (exp(-alpha1.*u)-1).*(exp(-alpha1.*v)-1)./(exp(-alpha1)-1))) + (1-w1).* ( (u.^(-alpha2) + v.^(-alpha2) - 1).^(-1./alpha2) );
          diffCopula_uv = @(u,v) (alpha1.*w1.*exp(-alpha1.*u).*exp(-alpha1.*v).*(exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./((exp(-alpha1) - 1).^2.*(((exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./(exp(-alpha1) - 1) + 1).^2) - (alpha1.*w1.*exp(-alpha1.*u).*exp(-alpha1.*v))./((exp(-alpha1) - 1).*(((exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./(exp(-alpha1) - 1) + 1)) - (alpha2.*(1./alpha2 + 1).*(w1 - 1))./(u.^(alpha2 + 1).*v.^(alpha2 + 1).*(1./u.^alpha2 + 1./v.^alpha2 - 1).^(1./alpha2 + 2));
   
    case 'mixfrc180'   % MixFrC180 copula : Frank copula + rotated Clayton copula£¨180-degrees£©
          w1 = weight1;                 % weight of Frank copula
          alpha1 = copulaparameter1;    % parameter of Frank copula
          alpha2 = copulaparameter2;    % parameter of rotated Clayton copula£¨180-degrees£© 
          Copula = @(u,v) w1.*((-1./alpha1).*log(1 + (exp(-alpha1.*u)-1).*(exp(-alpha1.*v)-1)./(exp(-alpha1)-1))) + (1-w1).* ( u + v - 1 + ((1-u).^(-alpha2) + (1-v).^(-alpha2) - 1).^(-1./alpha2) );
          diffCopula_uv = @(u,v) (alpha1.*w1.*exp(-alpha1.*u).*exp(-alpha1.*v).*(exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./((exp(-alpha1) - 1).^2.*(((exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./(exp(-alpha1) - 1) + 1).^2) - (alpha1.*w1.*exp(-alpha1.*u).*exp(-alpha1.*v))./((exp(-alpha1) - 1).*(((exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./(exp(-alpha1) - 1) + 1)) - (alpha2.*(1./alpha2 + 1).*(w1 - 1))./((1 - u).^(alpha2 + 1).*(1 - v).^(alpha2 + 1).*(1./(1 - u).^alpha2 + 1./(1 - v).^alpha2 - 1).^(1./alpha2 + 2));
          
    case 'mixfg'   % MixFG copula : Frank copula + Gumbel copula
          w1 = weight1;                 % weight of Frank copula
          alpha1 = copulaparameter1;    % parameter of Frank copula
          alpha2 = copulaparameter2;    % parameter of Gumbel copula 
          Copula = @(u,v) w1.*((-1./alpha1).*log(1 + (exp(-alpha1.*u)-1).*(exp(-alpha1.*v)-1)./(exp(-alpha1)-1))) + (1-w1).* (exp(-((-log(u)).^alpha2 + (-log(v)).^alpha2).^(alpha2.^(-1))));
          diffCopula_uv = @(u,v) (alpha1.*w1.*exp(-alpha1.*u).*exp(-alpha1.*v).*(exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./((exp(-alpha1) - 1).^2.*(((exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./(exp(-alpha1) - 1) + 1).^2) - (exp(-((-log(u)).^alpha2 + (-log(v)).^alpha2).^(1./alpha2)).*(-log(u)).^(alpha2 - 1).*(-log(v)).^(alpha2 - 1).*((-log(u)).^alpha2 + (-log(v)).^alpha2).^(2./alpha2 - 2).*(w1 - 1))./(u.*v) - (alpha1.*w1.*exp(-alpha1.*u).*exp(-alpha1.*v))./((exp(-alpha1) - 1).*(((exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./(exp(-alpha1) - 1) + 1)) + (alpha2.*exp(-((-log(u)).^alpha2 + (-log(v)).^alpha2).^(1./alpha2)).*(-log(u)).^(alpha2 - 1).*(-log(v)).^(alpha2 - 1).*(1./alpha2 - 1).*((-log(u)).^alpha2 + (-log(v)).^alpha2).^(1./alpha2 - 2).*(w1 - 1))./(u.*v);
  
    case 'mixfrg180'   % MixFrG180 copula : Frank copula + rotated Gumbel copula£¨180-degrees£©
          w1 = weight1;                 % weight of Frank copula
          alpha1 = copulaparameter1;    % parameter of Frank copula
          alpha2 = copulaparameter2;    % parameter of rotated Gumbel copula (180-degree)
          Copula = @(u,v) w1.*((-1./alpha1).*log(1 + (exp(-alpha1.*u)-1).*(exp(-alpha1.*v)-1)./(exp(-alpha1)-1))) + (1-w1).* ( u + v - 1 + exp(-((-log(1-u)).^alpha2 + (-log(1-v)).^alpha2).^(alpha2.^(-1))) );
          diffCopula_uv = @(u,v) (alpha1.*w1.*exp(-alpha1.*u).*exp(-alpha1.*v).*(exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./((exp(-alpha1) - 1).^2.*(((exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./(exp(-alpha1) - 1) + 1).^2) - (alpha1.*w1.*exp(-alpha1.*u).*exp(-alpha1.*v))./((exp(-alpha1) - 1).*(((exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./(exp(-alpha1) - 1) + 1)) - ((exp(-((-log(1 - u)).^alpha2 + (-log(1 - v)).^alpha2).^(1./alpha2)).*(-log(1 - u)).^(alpha2 - 1).*(-log(1 - v)).^(alpha2 - 1).*((-log(1 - u)).^alpha2 + (-log(1 - v)).^alpha2).^(2./alpha2 - 2))./((u - 1).*(v - 1)) - (alpha2.*exp(-((-log(1 - u)).^alpha2 + (-log(1 - v)).^alpha2).^(1./alpha2)).*(-log(1 - u)).^(alpha2 - 1).*(-log(1 - v)).^(alpha2 - 1).*(1./alpha2 - 1).*((-log(1 - u)).^alpha2 + (-log(1 - v)).^alpha2).^(1./alpha2 - 2))./((u - 1).*(v - 1))).*(w1 - 1);
     
    case 'mixfgrg180'   % MixFGrG180 copula : Frank copula + Gumbel copula + rotated Gumbel copula£¨180-degrees£©
          w1 = weight1;                 % weight of Frank copula
          w2 = weight2;                 % weight of Gumbel copula
          alpha1 = copulaparameter1;    % parameter of Frank copula
          alpha2 = copulaparameter2;    % parameter of Gumbel copula
          alpha3 = copulaparameter3;    % parameter of rotated Gumbel copula (180-degrees)
          Copula = @(u,v) w1.*((-1./alpha1).*log(1 + (exp(-alpha1.*u)-1).*(exp(-alpha1.*v)-1)./(exp(-alpha1)-1)))  +  w2.*( exp(-((-log(u)).^alpha2 + (-log(v)).^alpha2).^(alpha2.^(-1))) ) + (1-w1-w2).*(u + v - 1 + exp(-((-log(1-u)).^alpha3 + (-log(1-v)).^alpha3).^(alpha3.^(-1))));
          diffCopula_uv = @(u,v) (w2.*exp(-((-log(u)).^alpha2 + (-log(v)).^alpha2).^(1./alpha2)).*(-log(u)).^(alpha2 - 1).*(-log(v)).^(alpha2 - 1).*((-log(u)).^alpha2 + (-log(v)).^alpha2).^(2./alpha2 - 2))./(u.*v) - (alpha1.*w1.*exp(-alpha1.*u).*exp(-alpha1.*v))./((exp(-alpha1) - 1).*(((exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./(exp(-alpha1) - 1) + 1)) - ((exp(-((-log(1 - u)).^alpha3 + (-log(1 - v)).^alpha3).^(1./alpha3)).*(-log(1 - u)).^(alpha3 - 1).*(-log(1 - v)).^(alpha3 - 1).*((-log(1 - u)).^alpha3 + (-log(1 - v)).^alpha3).^(2./alpha3 - 2))./((u - 1).*(v - 1)) - (alpha3.*exp(-((-log(1 - u)).^alpha3 + (-log(1 - v)).^alpha3).^(1./alpha3)).*(-log(1 - u)).^(alpha3 - 1).*(-log(1 - v)).^(alpha3 - 1).*(1./alpha3 - 1).*((-log(1 - u)).^alpha3 + (-log(1 - v)).^alpha3).^(1./alpha3 - 2))./((u - 1).*(v - 1))).*(w1 + w2 - 1) + (alpha1.*w1.*exp(-alpha1.*u).*exp(-alpha1.*v).*(exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./((exp(-alpha1) - 1).^2.*(((exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./(exp(-alpha1) - 1) + 1).^2) - (alpha2.*w2.*exp(-((-log(u)).^alpha2 + (-log(v)).^alpha2).^(1./alpha2)).*(-log(u)).^(alpha2 - 1).*(-log(v)).^(alpha2 - 1).*(1./alpha2 - 1).*((-log(u)).^alpha2 + (-log(v)).^alpha2).^(1./alpha2 - 2))./(u.*v); 

    % Here you can add other copula models

    otherwise
          error('unsupported copula model.');
    end   
    
    %% type II local Kendall's tau
    switch lower(region)
          case 'uu'  % upper-upper local Kendall's tau
               numerator_internal_term = @(u,v) (2 .* Copula(u,v) - Copula(u,current_q) - Copula(current_p,v) - u - v) .* diffCopula_uv(u,v); 
               numerator_term = 2 .* (integral2(numerator_internal_term,current_p,1,current_q,1,'RelTol',1e-6)); 
               denominator_term = (1 - current_p - current_q + Copula(current_p,current_q)).^2;
               constant_term = (1 + current_p + current_q + Copula(current_p,current_q))./(1 - current_p - current_q + Copula(current_p,current_q));
               ld = numerator_term./denominator_term + constant_term;
               
          case 'ul'  % upper-lower local Kendall's tau
               internal_term = @(u,v) (2 .* Copula(u,v) - Copula(current_p,v) - v) .* diffCopula_uv(u,v);
               numerator_term = 2 .* integral2(internal_term,current_p,1,0,current_q,'RelTol',1e-6);
               denominator_term = (current_q - Copula(current_p,current_q)).^2;
               ld = numerator_term ./ denominator_term;
               
          case 'lu'  % lower-upper local Kendall's tau
               internal_term = @(u,v) (2 .* Copula(u,v) - Copula(u,current_q) - u) .* diffCopula_uv(u,v);
               numerator_term = 2 .* integral2(internal_term, 0,current_p,current_q,1,'RelTol',1e-6);
               denominator_term = (current_p - Copula(current_p,current_q)).^2;
               ld = numerator_term ./ denominator_term;
               
          case 'll'  % lower-lower local Kendall's tau
               internal_term = @(u,v) Copula(u,v) .* diffCopula_uv(u,v); 
               ld = (4 .* integral2(internal_term,0,current_p,0,current_q,'RelTol',1e-6))./((Copula(current_p,current_q)).^2) - 1; 
              
          otherwise
               error('Error region type.');
    end
    
    disp("Calculating local Kendall's tau using copula model... output: " + ld);

otherwise
    error(message('Error measure type.'));
end


