function CL = mixGrG180CL(theta,data)
% Description: The negative log-likelihood of MixGrG180 copula model
%
% Written for paper "Generalized local Kendall’s τ: a novel framework for uncovering nonlinear local dependence" (Huang & Zhang,2026)
%
% Author: Zaixin Huang
% Date: finished at 2015.06.07; updated at 2018.01.04; current version: 2025.03.16
% Contact: For bug reports and suggestions, please contact me at eric.huangzaixin@gmail.com. 
%          I will update them on GitHub and acknowledge your contribution. Thank you!
% The latest version can be downloaded from https://github.com/huangzaixin/local-dependence-toolbox
%%
u = data(:,1);
v = data(:,2);

w1 = theta(1);      % weight of Gumbel copula
alpha1 = theta(2);  % Gumbel copula parameter
alpha2 = theta(3);  % rotated Gumbel copula parameter

CL = -sum(log((w1.*exp(-((-log(u)).^alpha1 + (-log(v)).^alpha1).^(1./alpha1)).*(-log(u)).^(alpha1 ...
     - 1).*(-log(v)).^(alpha1 - 1).*((-log(u)).^alpha1 + (-log(v)).^alpha1).^(2./alpha1 - 2))./(u.*v)...
     - ((exp(-((-log(1 - u)).^alpha2 + (-log(1 - v)).^alpha2).^(1./alpha2)).*(-log(1 - u)).^(alpha2 ...
     - 1).*(-log(1 - v)).^(alpha2 - 1).*((-log(1 - u)).^alpha2 + (-log(1 - v)).^alpha2).^(2./alpha2 ...
     - 2))./((u - 1).*(v - 1)) - (alpha2.*exp(-((-log(1 - u)).^alpha2 + (-log(1 ...
     - v)).^alpha2).^(1./alpha2)).*(-log(1 - u)).^(alpha2 - 1).*(-log(1 - v)).^(alpha2 ...
     - 1).*(1./alpha2 - 1).*((-log(1 - u)).^alpha2 + (-log(1 - v)).^alpha2).^(1./alpha2 ...
     - 2))./((u - 1).*(v - 1))).*(w1 - 1) - (alpha1.*w1.*exp(-((-log(u)).^alpha1 ...
     + (-log(v)).^alpha1).^(1./alpha1)).*(-log(u)).^(alpha1 - 1).*(-log(v)).^(alpha1 ...
     - 1).*(1./alpha1 - 1).*((-log(u)).^alpha1 + (-log(v)).^alpha1).^(1./alpha1 - 2))./(u.*v)));

 
