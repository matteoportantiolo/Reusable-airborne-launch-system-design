function[m0, m1, m2, m3, ms1, ms2, ms3] = Staging_optimization_3st(Is, eps, m_pay, DV)

% INPUT PARAMETERS
% Is = Vector of specific impulses of propulsion systems per stage starting from 1 [s]
% eps = Vector of structural mass index of stages starting from 1 [-]
% m_pay = Payload mass [kg]
% Delta V [m/s]

% Gravitational costant
g0 = 9.80665;       %[m^2/s]

% Effective exhaust velocity
c = Is*g0;

% Lagrangina multiplier computation
f = @(lam) c(1)*log((c(1)*lam-1)/(c(1)*eps(1)*lam)) + c(2)*log((c(2)*lam-1)/(c(2)*eps(2)*lam)) +...
           + c(3)*log((c(3)*lam-1)/(c(3)*eps(3)*lam)) - DV;

options = optimoptions('fsolve','TolFun',1e-6,'Display','None');
lambda = fsolve(f, 1, options);

% Optimal mass fractions
n = (c*lambda-1)./(c.*eps*lambda);

% Stage gross lift off masses
m3 = (n(3)-1)/(1-n(3)*eps(3))*m_pay;
m2 = (n(2)-1)/(1-n(2)*eps(2))*(m_pay + m3);
m1 = (n(1)-1)/(1-n(1)*eps(1))*(m_pay + m3 +m2);

ms1 = eps(1)*m1;
ms2 = eps(2)*m2;
ms3 = eps(3)*m3;

% Total gross lift off mass
m0 = m1 + m2 + m3 + m_pay;

end