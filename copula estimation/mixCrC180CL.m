function CL = mixCrC180CL(theta,data)
% Description: The negative log-likelihood of MixCrC180 copula model
%
% Written for paper "Generalized local Kendall’s τ: a novel framework 
%                    for uncovering nonlinear local dependence" (Huang & Zhang,2026)
%
% Author: Zaixin Huang
% Date: finished at 2015.06.07; current version: 2025.03.16
% Contact: For bug reports and suggestions, please contact me at eric.huangzaixin@gmail.com. 
%          I will update them on GitHub and acknowledge your contribution. Thank you!
% The latest version can be downloaded from https://github.com/huangzaixin/local-dependence-toolbox
%%
u = data(:,1);
v = data(:,2);
w1 = theta(1);       % weight of Clayton copula
alpha1 = theta(2);   % Clayton copula parameter
alpha2 = theta(3);   % rotated Clayton copula parameter (180-degrees)

CL = -sum(log((alpha1.*w1.*(1./alpha1 + 1))./(u.^(alpha1 + 1).*v.^(alpha1 ...
     + 1).*(1./u.^alpha1 + 1./v.^alpha1 - 1).^(1./alpha1 + 2)) ...
     - (alpha2.*(1./alpha2 + 1).*(w1 - 1))./((1 - u).^(alpha2 ...
     + 1).*(1 - v).^(alpha2 + 1).*(1./(1 - u).^alpha2 + 1./(1 ...
     - v).^alpha2 - 1).^(1./alpha2 + 2))));
 
 