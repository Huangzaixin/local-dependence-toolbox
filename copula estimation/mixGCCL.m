function CL = mixGCCL(theta,data)
% Description: The negative log-likelihood of MixGC copula model
% Author: Zaixin Huang
% Date: finished at 2015.06.07; updated at 2018.01.04; this version: 2025.03.16
% Bug reports and suggestions: if you find any bugs or have suggestions, please contact me at eric.huangzaixin@gmail.com. 
%                              I will update them on GitHub and acknowledge your contribution. Thank you!
% The latest version can be downloaded from https://github.com/huangzaixin/local-dependence-toolbox
%%
u = data(:,1);
v = data(:,2);

w1 = theta(1);       % weight of gumbel copula
alpha1 = theta(2);   % gumbel copula parameter
alpha2 = theta(3);   % clayton copula parameter

CL = -sum(log((w1.*exp(-((-log(u)).^alpha1 + (-log(v)).^alpha1).^(1./alpha1)).*(-log(u)).^(alpha1 ...
     - 1).*(-log(v)).^(alpha1 - 1).*((-log(u)).^alpha1 + (-log(v)).^alpha1).^(2./alpha1 - 2))./(u.*v) ...
     - (alpha2.*(1./alpha2 + 1).*(w1 - 1))./(u.^(alpha2 + 1).*v.^(alpha2 + 1).*(1./u.^alpha2 ...
     + 1./v.^alpha2 - 1).^(1./alpha2 + 2)) - (alpha1.*w1.*exp(-((-log(u)).^alpha1 ...
     + (-log(v)).^alpha1).^(1./alpha1)).*(-log(u)).^(alpha1 - 1).*(-log(v)).^(alpha1 ...
     - 1).*(1./alpha1 - 1).*((-log(u)).^alpha1 + (-log(v)).^alpha1).^(1./alpha1 - 2))./(u.*v)));
    


