function [MASS, m_struct,m_prop] = OptimalStagingFun(eps1, eps2, isp1, isp2, DV, mpay)

% DV in m/s  !!!!
% mpay in Kg !!!

% DATA USED SO FAR

% eps1 = 0.15;
% eps2 = 0.13;
% From CEA analysis
% isp1 = 323.16;
% isp2 = 336.58;
% DV = 9.17e3;
% mpay = 250;

% computation of the effective exaust velocities
g0 = 9.80665;  
c1 = isp1*g0;
c2 = isp2*g0;

% Solving for lambda:=l
f = @(l) DV-c1*log((c1*l-1)/(eps1*c1*l))-c2*log((c2*l-1)/(eps2*c2*l));
options = optimoptions('fsolve','FunctionTolerance',1e-12,'Display','off');
l =fsolve(f,1,options);

% recovering n1,n2,n3
% defining eps and c vectors
eps = [eps1; eps2];
c   = [c1; c2];
n = zeros(2,1);
for i=1:2
    n(i) = (c(i)*l-1)/(eps(i)*c(i)*l);
end

% once n1, n2 have been computed, it's possible to recover the stage
% masses according to the following equations:
% STAGE MASS: Propellant mass of the stage+structural mass of the stage
MASS.m2 = ((n(2)-1)/(1-n(2)*eps2))*mpay; % second stage mass
MASS.m1 = ((n(1)-1)/(1-n(1)*eps1))*(mpay+MASS.m2); % first stage mass
% GROSS LIFT OFF MASS 
MASS.glom = MASS.m1+MASS.m2+mpay; % [Kg]

m_struct = zeros(2,1);
m_prop   = zeros(2,1);
m = [MASS.m1;MASS.m2];
for i=1:2
    m_struct(i) = m(i)*eps(i);
    m_prop(i) = m(i)*(1-eps(i));
end

fprintf('\n2STO LIQUID')
fprintf('\n gross liftoff mass = %f kg \n',MASS.glom );
fprintf('\n first stage mass = %f Kg ',MASS.m1);
fprintf('\n second stage mass = %f kg \n',MASS.m2);
fprintf('\n first stage propellant mass = %f kg',m_prop(1));
fprintf('\n second stage propellant mass = %f kg \n',m_prop(2));


end


