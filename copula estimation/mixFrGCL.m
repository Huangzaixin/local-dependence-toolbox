function CL = mixFrGCL(theta,data)
% Description: The negative log-likelihood of MixFrG copula model
% Author: Zaixin Huang
% Date: finished at 2015.06.07; updated at 2018.01.04; this version: 2025.03.16
% Bug reports and suggestions: if you find any bugs or have suggestions, please contact me at eric.huangzaixin@gmail.com. 
%                              I will update them on GitHub and acknowledge your contribution. Thank you!
% The latest version can be downloaded from https://github.com/huangzaixin/local-dependence-toolbox
%%
u = data(:,1);
v = data(:,2);

w1 = theta(1);       % weight of frank copula 
alpha1 = theta(2);   % frank copula parameter
alpha2 = theta(3);   % rotated clayton copula parameter

CL = -sum(log( (alpha1.*w1.*exp(-alpha1.*u).*exp(-alpha1.*v).*(exp(-alpha1.*u) - 1).*(exp(-alpha1.*v)...
     - 1))./((exp(-alpha1) - 1).^2.*(((exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./(exp(-alpha1) - 1) + 1).^2)...
     - (alpha1.*w1.*exp(-alpha1.*u).*exp(-alpha1.*v))./((exp(-alpha1) - 1).*(((exp(-alpha1.*u) - 1).*(exp(-alpha1.*v)...
     - 1))./(exp(-alpha1) - 1) + 1)) - ((exp(-((-log(1 - u)).^alpha2 + (-log(1 - v)).^alpha2).^(1./alpha2)).*(-log(1 ...
     - u)).^(alpha2 - 1).*(-log(1 - v)).^(alpha2 - 1).*((-log(1 - u)).^alpha2 + (-log(1 - v)).^alpha2).^(2./alpha2 ...
     - 2))./((u - 1).*(v - 1)) - (alpha2.*exp(-((-log(1 - u)).^alpha2 + (-log(1 - v)).^alpha2).^(1./alpha2)).*(-log(1 ...
     - u)).^(alpha2 - 1).*(-log(1 - v)).^(alpha2 - 1).*(1./alpha2 - 1).*((-log(1 - u)).^alpha2 + (-log(1 ...
     - v)).^alpha2).^(1./alpha2 - 2))./((u - 1).*(v - 1))).*(w1 - 1) ));

% Frank copula model
% f = -1./alpha1.*log(1 + (exp(-alpha1.*u)-1).*(exp(-alpha1.*v)-1)./(exp(-alpha1)-1))

% Rotated Gumbel copula model (180-degrees)
% f = u + v - 1 + ((1-u).^(-alpha2) + (1-v).^(alpha2)-1).^(-1./alpha2)

% MixFrG copula model
% f = w1.*( -1./alpha1.*log(1 + (exp(-alpha1.*u)-1).*(exp(-alpha1.*v)-1)./(exp(-alpha1)-1)) ) + (1-w1).*( u + v - 1 + ((1-u).^(-alpha2) + (1-v).^(alpha2)-1).^(-1./alpha2) )
% dudv = diff(diff(f, u),v) = (alpha1.*w1.*exp(-alpha1.*u).*exp(-alpha1.*v).*(exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./((exp(-alpha1) - 1).^2.*(((exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./(exp(-alpha1) - 1) + 1).^2) - (alpha1.*w1.*exp(-alpha1.*u).*exp(-alpha1.*v))./((exp(-alpha1) - 1).*(((exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) - 1))./(exp(-alpha1) - 1) + 1)) - ((exp(-((-log(1 - u)).^alpha2 + (-log(1 - v)).^alpha2).^(1./alpha2)).*(-log(1 - u)).^(alpha2 - 1).*(-log(1 - v)).^(alpha2 - 1).*((-log(1 - u)).^alpha2 + (-log(1 - v)).^alpha2).^(2./alpha2 - 2))./((u - 1).*(v - 1)) - (alpha2.*exp(-((-log(1 - u)).^alpha2 + (-log(1 - v)).^alpha2).^(1./alpha2)).*(-log(1 - u)).^(alpha2 - 1).*(-log(1 - v)).^(alpha2 - 1).*(1./alpha2 - 1).*((-log(1 - u)).^alpha2 + (-log(1 - v)).^alpha2).^(1./alpha2 - 2))./((u - 1).*(v - 1))).*(w1 - 1)


