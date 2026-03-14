function copula_par = fun_copula_estimation(copulatype,data)
% Description: Estimate parameters for Clayton, Gumbel, and Frank copula models.
% Inputs:
%     copulatype: the copula model type. Must be one of: 'clayton', 'gumbel', 'frank'.
%     data: An N-by-2 matrix of values in the open interval (0,1) representing the transformed uniform margins.
%
% Outputs:
%     copula_par: The estimated copula parameter.
%
% Written for paper "Generalized local Kendall’s τ: a novel framework for
% uncovering nonlinear local dependence", published in Biometrics (Huang & Zhang, 2026).
%
% Author: Zaixin Huang
% Date: Completed on 2023-01-01
% Contact: For bug reports and suggestions, please contact me at eric.huangzaixin@gmail.com. 
%          I will update them on GitHub and acknowledge your contribution. 
% Repository: https://github.com/huangzaixin/local-dependence-toolbox
%%
% Set Optimization Options
options = optimset('Display','off','TolCon',10^-12,'TolFun',10^-9,'TolX',10^-9);

u = data(:,1);
v = data(:,2);

if any(u < 0 | u > 1) || any(v < 0 | v > 1)
    error('Data values must be in the range [0, 1]');
end
    
% Estimate different copula models
switch lower(copulatype)
    case 'clayton'
        lower_par = 0.001;
        theta0 = 1.001;
        [copula_par] = fmincon('claytonCL',theta0,[],[],[],[],lower_par,[],[],options,[u,v]);
    case 'gumbel'
        lower_par = 1.001;
        theta0 = 2;
        [copula_par] = fmincon('gumbelCL',theta0,[],[],[],[],lower_par,[],[],options,[u,v]);
    case 'frank'
        lower_par = 0.001;
        theta0 = 1;
        [copula_par] = fmincon('frankCL',theta0,[],[],[],[],lower_par,[],[],options,[u,v]);  
    otherwise
        error('Unsupported copula type: %s. Supported types are: ''Clayton'', ''Gumbel'', ''Frank''', copulatype);    
end

end
