clear
close all
clc

%% STAGE 1
clear
close all
clc

% Data stage 1
m0 = 30e3;
n = 7;
Is = 330;
T = 300e3;
d = 5;
A = pi*d^2/4;
T2W = 1.4;
Cd = 0.5;
Cl = 0;
delta = 0;

g0 = 9.81;
m_dot = T/(Is*g0);
mf = m0/n;
m_prop = m0-mf;

t_burn = m_prop/m_dot;

t0 = 0;
tf = t_burn;
tspan = [t0, tf];

v0 = 2000;
gamma0 = deg2rad(49);
h0 = 20e3;
omega0 = deg2rad(0);
y0 = [v0;gamma0;h0;omega0];

param.T = T;
param.A = A;
param.Cd = Cd;
param.Cl = Cl;
param.m0 = m0;
param.t_burn = t_burn;
param.m_dot = m_dot;
param.delta = delta;

% Data stage 2
m02 = 5500;
T2 = 280e3;
Is2 = 350;
n2 = 3;
mf2 = m02/n2;
m_prop2 = m02-mf2;
m_dot2 = T2/(Is2*g0);
t_burn2 = m_prop2/m_dot2;

t02 = 0;
tf2 = t_burn2;
tspan2 = [t02, tf2];

param2.T = T2;
param2.A = A;
param2.Cd = Cd;
param2.Cl = Cl;
param2.m0 = m02;
param2.t_burn = t_burn2;
param2.m_dot = m_dot2;
param2.delta = delta;

% Perform integration
%options = odeset('reltol', 1e-12, 'abstol', 1e-12);
[tt, xx] = ode45(@(t,x) dynamics(t,x,param), tspan, y0);
[tt2, xx2] = ode45(@(t,x) dynamics(t,x,param2), tspan2, xx(end,:));

tt_update = [tt;tt(end)+tt2];
xx_update = [xx;xx2];

figure
plot(tt_update,xx_update(:,1),'LineWidth',2)
xlabel('t [-]',FontSize=14)
ylabel('y [-]',FontSize=14)
grid on
title('Velocity')

figure
plot(tt_update,rad2deg(xx_update(:,2)),'LineWidth',2)
xlabel('t [-]',FontSize=14)
ylabel('gamma [-]',FontSize=14)
grid on
title('Flight path angle')

figure
plot(tt_update,xx_update(:,3),'LineWidth',2)
xlabel('t [-]',FontSize=14)
ylabel('h [-]',FontSize=14)
grid on
title('Altitude')

figure
plot(tt_update,xx_update(:,4),'LineWidth',2)
xlabel('t [-]',FontSize=14)
ylabel('omega [-]',FontSize=14)
grid on
title('Angular position')




%% FUNCTIONS

function dydt = dynamics(t,y,param)

    dydt = NaN(length(y),1);

    v =     y(1); %[m/s]
    gamma = y(2); %[rad]
    h =     y(3); %[m]
    %omega = y(4);%[rad]

    T = param.T;
    A = param.A;
    Cd = param.Cd;
    Cl = param.Cl;
    m0 = param.m0;
    delta = param.delta;
    t_burn = param.t_burn;
    m_dot = param.m_dot;

    R_e = 6378.1366e3;
    g0 = 9.81;
    g = g0/(1+h/R_e)^2;

    rho0 = 1.225;
    h_scale = 7.5e3;
    rho = rho0*exp(-h/h_scale);
    D = 0.5*rho*v^2*Cd*A;
    L = 0.5*rho*v^2*Cl*A;

    if t <= t_burn

        m = m0-m_dot*t;

    elseif t > t_burn

        T = 0;
        m = m0-m_dot*t_burn;

    end

    dydt(1) = T/m*cos(delta) - D/m - g*sin(gamma);
    dydt(2) = (v*cos(gamma))/(R_e+h) + (T*sin(delta))/(m*v) + L/m - g*cos(gamma)/v;
    dydt(3) = v*sin(gamma);
    dydt(4) = v*cos(gamma)/(R_e+h);

end




