function CL = rotClayton180CL(theta,data)
% Description: The negative log-likelihood of rotated Clayton copula model.
%
% Written for paper "Generalized local Kendall’s τ: a novel framework for
% uncovering nonlinear local dependence", published in Biometrics (Huang & Zhang, 2026).
%
% Author: Zaixin Huang
% Date: completed on 2015-06-07; updated on 2023-01-01; current version: 2025-03-16
% Contact: For bug reports and suggestions, please contact me at eric.huangzaixin@gmail.com. 
%          I will update them on GitHub and acknowledge your contribution.
% Repository: https://github.com/huangzaixin/local-dependence-toolbox
%%
u = data(:,1);
v = data(:,2);

alpha = theta;       % rotated Clayton copula parameter (180-degrees)

CL = -sum(log( (alpha.*(1./alpha + 1))./((1 - u).^(alpha + 1).*(1 ...
        - v).^(alpha + 1).*(1./(1 - u).^alpha + 1./(1 - v).^alpha - 1).^(1./alpha + 2))  ));

% rotated clayton copula model (180-degrees)
% f = u + v - 1 + ((1-u).^(-alpha) + (1-v).^(-alpha) - 1).^(-1./alpha)
% dudv = diff(diff(f, u),v) = (alpha.*(1./alpha + 1))./((1 - u).^(alpha + 1).*(1 - v).^(alpha + 1).*(1./(1 - u).^alpha + 1./(1 - v).^alpha - 1).^(1./alpha + 2))


