function CL = mixFCrC180CL(theta,data)
% Description: The negative log-likelihood of MixFCrC180 copula model
%
% Written for paper "Generalized local Kendall’s τ: a novel framework for uncovering nonlinear local dependence" (Huang & Zhang,2026)
%
% Author: Zaixin Huang
% Date: 2025.03.16
% Contact: For bug reports and suggestions, please contact me at eric.huangzaixin@gmail.com. 
%          I will update them on GitHub and acknowledge your contribution. Thank you!
% The latest version can be downloaded from https://github.com/huangzaixin/local-dependence-toolbox
%%
u = data(:,1);
v = data(:,2);

w1 = theta(1);      % weight of Frank copula 
w2 = theta(2);      % weight of Clayton copula 

alpha1 = theta(3);  % Frank copula parameter
alpha2 = theta(4);  % Clayton copula parameter
alpha3 = theta(5);  % rotated Clayton copula parameter (180-degrees)

CL = -sum(log( (alpha2.*w2.*(1./alpha2 + 1))./(u.^(alpha2 + 1).*v.^(alpha2 + 1).*(1./u.^alpha2...
     + 1./v.^alpha2 - 1).^(1./alpha2 + 2)) - (alpha3.*(1./alpha3 + 1).*(w1 + w2 - 1))./((1 - u).^(alpha3...
     + 1).*(1 - v).^(alpha3 + 1).*(1./(1 - u).^alpha3 + 1./(1 - v).^alpha3 - 1).^(1./alpha3 + 2)) ...
     - (alpha1.*w1.* exp(-alpha1.*u).* exp(-alpha1.*v))./((exp(-alpha1) - 1).*(((exp(-alpha1.*u) ...
     - 1).*(exp(-alpha1.*v) - 1))./(exp(-alpha1) - 1) + 1)) + (alpha1.*w1.*exp(-alpha1.*u).*exp(-alpha1.*v).*(exp(-alpha1.*u)...
     - 1).*(exp(-alpha1.*v) - 1))./((exp(-alpha1) - 1).^2.*(((exp(-alpha1.*u) - 1).*(exp(-alpha1.*v) ...
     - 1))./(exp(-alpha1) - 1) + 1).^2) ));

% f = w1.*( -1./alpha1.*log(1 + (exp(-alpha1.*u)-1).*(exp(-alpha1.*v)-1)./(exp(-alpha1)-1)) ) + w2.*(((u).^(-alpha2) + (v).^(-alpha2) - 1).^(-1./alpha2)) + (1-w1-w2).*( u + v - 1 + ((1-u).^(-alpha3) + (1-v).^(-alpha3) - 1).^(-1./alpha3) )

