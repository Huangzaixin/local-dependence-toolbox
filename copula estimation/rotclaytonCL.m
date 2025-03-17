function CL = rotclaytonCL(theta,data)
u = data(:,1);
v = data(:,2);

alpha = theta;       % rotated clayton copula parameter (180-degrees)

CL = -sum(log( (alpha.*(1./alpha + 1))./((1 - u).^(alpha + 1).*(1 ...
     - v).^(alpha + 1).*(1./(1 - u).^alpha + 1./(1 - v).^alpha - 1).^(1./alpha + 2))  ));

% rotated clayton copula model (180-degrees)
% f = u + v - 1 + ((1-u).^(-alpha) + (1-v).^(-alpha) - 1).^(-1./alpha)
% dudv = diff(diff(f, u),v) = (alpha.*(1./alpha + 1))./((1 - u).^(alpha + 1).*(1 - v).^(alpha + 1).*(1./(1 - u).^alpha + 1./(1 - v).^alpha - 1).^(1./alpha + 2))

