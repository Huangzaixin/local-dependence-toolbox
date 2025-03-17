function copula_par = fun_copula_estimation(copulatype,data)
% Description: estimate clayton, gumbel and frank copula models
% Author: Zaixin Huang
% Date: finished at 2023.01.01
% Bug reports and suggestions: 
%       if you find any bugs or have suggestions, please contact me at eric.huangzaixin@gmail.com. 
%       I will update them on GitHub and acknowledge your contribution. Thank you!
% The latest version can be downloaded from https://github.com/huangzaixin/local-dependence-toolbox
%%
options = optimset('Display','off','TolCon',10^-12,'TolFun',10^-9,'TolX',10^-9);
u = data(:,1);
v = data(:,2);
    switch copulatype
        case 'clayton'
            lower = 0.001;
            theta0 = 1.001;
            [copula_par] = fmincon('claytonCL',theta0,[],[],[],[],lower,[],[],options,[u,v]);
        case 'gumbel'
            lower = 1.001;
            theta0 = 2;
            [copula_par] = fmincon('gumbelCL',theta0,[],[],[],[],lower,[],[],options,[u,v]);
        case 'frank'
            lower = 0.001;
            theta0 = 1;
            [copula_par] = fmincon('frankCL',theta0,[],[],[],[],lower,[],[],options,[u,v]);   
    end
end