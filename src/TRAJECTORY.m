clear
clc
close all
format long


% OPTIMIZATION OF 2 + ITERATION OF 1, WITH ENGINES ON/OFF


% 2STO LIQUID
%  gross liftoff mass = 16435.352189 kg 
% 
%  first stage mass = 14265.908201 Kg 
%  second stage mass = 1769.443988 kg 
% 
%  first stage propellant mass = 12553.999217 kg
%  second stage propellant mass = 1486.332950 kg 
% 
% 3ST0 LIQUID ANALYSIS
%  gross liftoff mass = 13239.950384 kg 
% 
%  first stage mass = 9188.370409 Kg 
%  second stage mass = 2916.643427 kg
%  third stage mass = 734.936547 kg 
% 
%  first stage propellant mass = 8085.765960 kg
%  second stage propellant mass = 2508.313347 kg
%  third stage propellant mass = 609.997334 kg 


% TRAJECTORY PHASES
% 
% 1) Release: h = 39kfeet, M = 0.82, t = 5s
% 2) Stage 1 ignition: Tmax = 163kpounds, up to 53.9km
% 3) Stage 1 separation: after 16sec
% 3) Stage 2 ignition: after 1sec, 72.1km 
% 4) Fairing separation: 116km, T = 44kpounds
% 5) Coasting phase: 4min
% 6) Stage 2 separation
% 7) Stage 3 ignition: 69sec
% TOT) t_tot = 11min20sec


% Consider different delta for different leg (max 3/4)


%% INITIALIZATION

t_i = 0;

D = 1.5;                    % [m]
        
param_1st.A = pi*(D/2)^2;   % [m^2]
param_1st.Cd = 0.5;         % [-]
param_1st.Cl = 0;           % [-]
param_1st.m0 = 7188.370409+400;    % [kg]
param_1st.I_sp = 320;       % [s] 
param_1st.delta = 0;        % [rad]
param_1st.m_prop = 6085.765960;
param_1st.n = 1/(-param_1st.m_prop/param_1st.m0 + 1);
param_1st.T= param_1st.m0*9.81;  % [N] %T/W=1

param_2st.A = pi*(D/2)^2;   % [m^2]
param_2st.Cd = 0.5;         % [-]
param_2st.Cl = 0;           % [-]
param_2st.m0 = 2916.643427+400;    % [kg]
param_2st.I_sp = 320;       % [s] 
param_2st.delta = 0;        % [rad]
param_2st.m_prop = 2508.313347;
param_2st.n = 1/(-param_2st.m_prop/param_2st.m0 + 1);
param_2st.T= param_2st.m0*9.81;  % [N] %T/W=1

param_3st.A = pi*(D/2)^2;   % [m^2]
param_3st.Cd = 0.5;         % [-]
param_3st.Cl = 0;           % [-]
param_3st.m0 = 734.936547+400;     % [kg]
param_3st.I_sp = 320;       % [s] 
param_3st.delta = 0;        % [rad]
param_3st.m_prop = 609.997334;
param_3st.n = 1/(-param_3st.m_prop/param_3st.m0 + 1);
param_3st.T= param_3st.m0*9.81;  % [N] %T/W=1

t_burn1 = param_1st.m_prop*param_1st.I_sp*9.81/param_1st.T;
param_1st.t_burn = t_burn1;
t_burn2 = param_2st.m_prop*param_2st.I_sp*9.81/param_2st.T;
param_2st.t_burn = t_burn2;
t_burn3 = param_3st.m_prop*param_3st.I_sp*9.81/param_3st.T;
param_3st.t_burn = t_burn3;

%% Perform integration (initial concept)

% v0 = 2000;  % [m/s]
% gamma0 = deg2rad(33); % [rad]
% omega0 =     0; % [m]
% h0 =     80e3; % [m]
% 
% y0 = [v0;gamma0;h0;omega0];
%
% [tt, xx] = ode78(@(t,x) fun(t,x,param_2st), [t_i t_burn2], y0);
% [tt2, xx2] = ode78(@(t,x) fun(t,x,param_3st), [0 t_burn3], xx(end,:));
% 
% tt_update = [tt;tt2+tt(end)];
% xx_update = [xx;xx2];
% 
% figure(1);
% subplot(2,2,1)
% plot(tt_update,xx_update(:,1))
% xlabel('t [s]','Interpreter','latex',FontSize=14)
% ylabel('v [m/s]','Interpreter','latex',FontSize=14)
% title("Variation of velocity in time ",'Interpreter','latex');
% grid on;
% 
% figure(1);
% subplot(2,2,2)
% plot(tt_update,rad2deg(xx_update(:,2)) )
% xlabel('t [s]','Interpreter','latex',FontSize=14)
% ylabel('$\gamma$ [deg]','Interpreter','latex',FontSize=14)
% title("Variation of $\gamma$ in time ",'Interpreter','latex');
% grid on;
% 
% subplot(2,2,3)
% plot(tt_update,xx_update(:,3)./1000)
% xlabel('t [s]','Interpreter','latex',FontSize=14)
% ylabel('h [km]','Interpreter','latex',FontSize=14)
% title("Altitude as function of time ",'Interpreter','latex');
% grid on;
% 
% subplot(2,2,4)
% plot(tt_update,xx_update(:,4))
% xlabel('t [s]','Interpreter','latex',FontSize=14)
% ylabel('$\omega$ [km]','Interpreter','latex',FontSize=14)
% title("Angular position as function of time ",'Interpreter','latex');
% grid on;

%% Optimization process of second leg

v0 = 2000;  % [m/s]
gamma0 = deg2rad(33); % [rad]
omega0 =     0; % [m]
h0 =     100e3; % [m]

y0 = [v0;gamma0;h0;omega0];

delta_t0 = 1;
y0_new = [y0(1:3);delta_t0];

options = optimoptions("fmincon","Algorithm","sqp",'Display','iter-detailed',...
                      'SpecifyObjectiveGradient',false,'SpecifyConstraintGradient',false);

%options.StepTolerance = 1.000000e-10;
options.ConstraintTolerance = 1e-10;
options.MaxFunctionEvaluations = 5e3;
options.MaxIterations = 800;
options.StepTolerance = 1e-08;
lb = [ 1200, 0.001, 80e3, 0];
ub = [ 2500, deg2rad(70), 120e3, 600]; 

[x_min2,f] = fmincon(@(x)ObjFun(x,param_2st,param_3st), y0_new,...
                                  [],[],[],[],lb,ub, @(x)myconstraints(x,param_2st,param_3st),...
                                  options);

delta_v_opt = f;

xx_0_NEW  = [x_min2(1:3);0];
delta_t_min = x_min2(4);

fprintf([' \n The optimized initial state is :\n initial velocity : %+.5e [m/s] \n initial pitch angle : %+.5e [deg]'  ...
    '\n initial height : %+.5e [km] \n \n delta_t off : %+.5e [s] \n'],xx_0_NEW(1),rad2deg(xx_0_NEW(2)),xx_0_NEW(3)./1e3,delta_t_min)

[tt2, xx2] = ode78(@(t,x) fun(t,x,param_2st), [t_i t_burn2], xx_0_NEW);
[tt_off, xx_off] = ode78(@(t,x) fun_off(t,x,param_2st), [0 delta_t_min], xx2(end,:));
[tt3, xx3] = ode78(@(t,x) fun(t,x,param_3st), [0 t_burn3], xx_off(end,:));

tt_update = [tt2;tt_off + tt2(end);tt3 + tt_off(end)+tt2(end)];
xx_update = [xx2;xx_off;xx3];

figure(2);
subplot(2,2,1)
plot(tt_update,xx_update(:,1))
xlabel('t [s]','Interpreter','latex',FontSize=14)
ylabel('v [m/s]','Interpreter','latex',FontSize=14)
title("Variation of velocity in time ",'Interpreter','latex');
grid on;

subplot(2,2,2)
plot(tt_update,rad2deg(xx_update(:,2)) )
xlabel('t [s]','Interpreter','latex',FontSize=14)
ylabel('$\gamma$ [deg]','Interpreter','latex',FontSize=14)
title("Variation of $\gamma$ in time ",'Interpreter','latex');
grid on;

subplot(2,2,3)
plot(tt_update,xx_update(:,3)./1000)
xlabel('t [s]','Interpreter','latex',FontSize=14)
ylabel('h [km]','Interpreter','latex',FontSize=14)
title("Altitude as function of time ",'Interpreter','latex');
grid on;

subplot(2,2,4)
plot(tt_update,xx_update(:,4))
xlabel('t [s]','Interpreter','latex',FontSize=14)
ylabel('$\omega$ [deg]','Interpreter','latex',FontSize=14)
title("Angular position as function of time ",'Interpreter','latex');
grid on;

%% Iteration process of first leg

% matching condition
param_1st.v_f = x_min2(1);
param_1st.gamma_f = x_min2(2);
param_1st.h_f = x_min2(3);

% guess
v0 = 150;  % [m/s]
gamma0 = deg2rad(1); % [rad]
omega0 =     0; % [m]
h0 =     11e3; % [m]
y0 = [v0;gamma0;h0;omega0];

t_bo = 180;
y0_new = [y0(1:3);t_bo];

options = optimoptions("fmincon","Algorithm","sqp",'Display','iter-detailed',...
                      'SpecifyObjectiveGradient',false,'SpecifyConstraintGradient',false);

%options.StepTolerance = 1.000000e-10;
options.ConstraintTolerance = 1e-10;
options.MaxFunctionEvaluations = 5e3;
options.MaxIterations = 800;
options.StepTolerance = 1e-08;
lb = [ 100, 0.001, 7e3, 0];
ub = [ 250, deg2rad(4), 13e3, 250]; 

[x_min,f,exflag,output] = fmincon(@(x)ObjFun_pitchup(x,param_1st), y0_new,...
                                  [],[],[],[],lb,ub, @(x)myconstraints_pitchup(x,param_1st),...
                                  options);

xx_0_NEW  = [x_min(1:3);0];
t_bo = x_min(4);

fprintf([' \n The optimized initial state is :\n initial velocity : %+.5e [m/s] \n initial pitch angle : %+.5e [deg]'  ...
    '\n initial height : %+.5e [km] \n \n born out time: %+.5e [s] \n'],x_min(1),rad2deg(x_min(2)),x_min(3)./1e3,t_bo)

[tt1, xx1] = ode78(@(t,x) fun_pitchup(t,x,param_1st,t_bo), [t_i t_bo], xx_0_NEW);

v_vect1 =       xx1(:,1);
gamma_vect1 =   xx1(:,2);
h_vect1 =       xx1(:,3);

% delta_v losses
Cb1 = param_1st.m0/(param_1st.A*param_1st.Cd);
H_scale = 7500; %[m]
R_E = 6.37813660e+06; %[m]

rho =   @(h) 1.225*exp(-h./H_scale);
m1 =    @(t) param_1st.m0 - param_1st.T./(param_1st.I_sp*9.81).*t;
g =     @(h) 9.81/(1+h./R_E).^2;

term1_d = -0.5./Cb1.*rho(h_vect1).*v_vect1.^2.*param_1st.m0./m1(tt1);
term1_g = -g(h_vect1)'.*sin(gamma_vect1);

delta_v_drag1 = abs(trapz(tt1,term1_d));
delta_v_grav1 = abs(trapz(tt1,term1_g));
delta_v_tot1 = delta_v_drag1 + delta_v_grav1;

delta_v_TOT = delta_v_tot1 + delta_v_opt;

fprintf(' \n The total delta_v losses are : %+.5e [m/s] \n',delta_v_TOT)

% plot
figure(3);
subplot(2,2,1)
plot(tt1,xx1(:,1))
xlabel('t [s]','Interpreter','latex',FontSize=14)
ylabel('v [m/s]','Interpreter','latex',FontSize=14)
title("Variation of velocity in time ",'Interpreter','latex');
grid on;

subplot(2,2,2)
plot(tt1,rad2deg(xx1(:,2)) )
xlabel('t [s]','Interpreter','latex',FontSize=14)
ylabel('$\gamma$ [deg]','Interpreter','latex',FontSize=14)
title("Variation of $\gamma$ in time ",'Interpreter','latex');
grid on;

subplot(2,2,3)
plot(tt1,xx1(:,3)./1000)
xlabel('t [s]','Interpreter','latex',FontSize=14)
ylabel('h [km]','Interpreter','latex',FontSize=14)
title("Altitude as function of time ",'Interpreter','latex');
grid on;

subplot(2,2,4)
plot(tt1,xx1(:,4))
xlabel('t [s]','Interpreter','latex',FontSize=14)
ylabel('$\omega$ [deg]','Interpreter','latex',FontSize=14)
title("Angular position as function of time ",'Interpreter','latex');
grid on;

%% FINAL PLOT

tt_tot = [tt1;tt_update+tt1(end)];
xx_tot = [xx1;xx_update];

figure(4);
subplot(2,2,1)
plot(tt_tot,xx_tot(:,1))
xlabel('t [s]','Interpreter','latex',FontSize=14)
ylabel('v [m/s]','Interpreter','latex',FontSize=14)
title("Variation of velocity in time ",'Interpreter','latex');
grid on;

subplot(2,2,2)
plot(tt_tot,rad2deg(xx_tot(:,2)) )
xlabel('t [s]','Interpreter','latex',FontSize=14)
ylabel('$\gamma$ [deg]','Interpreter','latex',FontSize=14)
title("Variation of $\gamma$ in time ",'Interpreter','latex');
grid on;

subplot(2,2,3)
plot(tt_tot,xx_tot(:,3)./1000)
xlabel('t [s]','Interpreter','latex',FontSize=14)
ylabel('h [km]','Interpreter','latex',FontSize=14)
title("Altitude as function of time ",'Interpreter','latex');
grid on;

subplot(2,2,4)
plot(tt_tot,xx_tot(:,4))
xlabel('t [s]','Interpreter','latex',FontSize=14)
ylabel('$\omega$ [deg]','Interpreter','latex',FontSize=14)
title("Angular position as function of time ",'Interpreter','latex');
grid on;

%% ---------------------------FUNCTIONS--------------------------------

function [dydt] = fun(t,y,param)

    v =     y(1); % [m/s]
    gamma = y(2); % [rad]
    h =     y(3); % [m]
    
    T =  param.T;
    A = param.A;
    Cd = param.Cd;
    Cl = param.Cl;
    m0 = param.m0;
    n = param.n;
    I_sp = param.I_sp; 
    % gamma_targ = param.gamma_targ; % [rad]
    % gamma_rate = param.gamma_rate; % [rad]
    delta = param.delta; % [rad]
    
    H_scale = 7500; %[m]
    R_E = 6.37813660e+06; %[m]
    g = 9.81/(1+h/R_E)^2;
    rho = 1.225*exp(-h/H_scale);
    D = 0.5*rho*(v^2)*A*Cd;
    L = 0.5*rho*(v^2)*A*Cl;
    m_prop0 = m0*(1-1/n);
    t_burn = m_prop0*I_sp*9.81/T;
    
    % Initializing the memory
    dydt = NaN(length(y),1);

    if t < t_burn
        m = m0 - T/(I_sp*9.81)*t;

    else
        T = 0;
        m = m0 -m_prop0;
    end

    dydt(1) = T/m*cos(delta) -D/m -g*sin(gamma);
    dydt(2) = v*cos(gamma)/(R_E+h) + T*sin(delta)/(m*v) + L/m -g*cos(gamma)/v ;
    dydt(3) = v*sin(gamma);
    dydt(4) = v*cos(gamma)/(R_E+h);
        
end

function [dydt] = fun_pitchup(t,y,param,t_bo)

    v =     y(1); % [m/s]
    gamma = y(2); % [rad]
    h =     y(3); % [m]
    
    T =  param.T;
    A = param.A;
    Cd = param.Cd;
    %Cl = param.Cl;
    m0 = param.m0;
    n = param.n;
    I_sp = param.I_sp; 
    % gamma_targ = param.gamma_targ; % [rad]
    %gamma_rate = param.gamma_rate; % [rad]
    delta = param.delta; % [rad]
    
    H_scale = 7500; %[m]
    R_E = 6.37813660e+06; %[m]
    g = 9.81/(1+h/R_E)^2;
    rho = 1.225*exp(-h/H_scale);
    D = 0.5*rho*(v^2)*A*Cd;
    %L = 0.5*rho*(v^2)*A*Cl;
    m_prop0 = m0*(1-1/n);
    t_burn = m_prop0*I_sp*9.81/T;

    gamma_rate = deg2rad(4.46154e+01)/t_bo;
    
    % Initializing the memory
    dydt = NaN(length(y),1);

    if t < t_burn
        m = m0 - T/(I_sp*9.81)*t;

    else
        T = 0;
        m = m0 -m_prop0;
    end

    dydt(1) = T/m*cos(delta) -D/m -g*sin(gamma);
    dydt(2) = gamma_rate ;
    dydt(3) = v*sin(gamma);
    dydt(4) = v*cos(gamma)/(R_E+h);
        
end

function [dydt] = fun_off(~,y,param)

    v =     y(1); % [m/s]
    gamma = y(2); % [rad]
    h =     y(3); % [m]
    
    A = param.A;
    Cd = param.Cd;
    Cl = param.Cl;
    m0 = param.m0;
    m_prop = param.m_prop;
    % gamma_targ = param.gamma_targ; % [rad]
    % gamma_rate = param.gamma_rate; % [rad]
    delta = param.delta; % [rad]
    
    H_scale = 7500; %[m]
    R_E = 6.37813660e+06; %[m]
    g = 9.81/(1+h/R_E)^2;
    rho = 1.225*exp(-h/H_scale);
    D = 0.5*rho*(v^2)*A*Cd;
    L = 0.5*rho*(v^2)*A*Cl;
    
    % Initializing the memory
    dydt = NaN(length(y),1);

    m = m0-m_prop;
    T = 0;

    dydt(1) = T/m*cos(delta) -D/m -g*sin(gamma);
    dydt(2) = v*cos(gamma)/(R_E+h) + T*sin(delta)/(m*v) + L/m -g*cos(gamma)/v ;
    dydt(3) = v*sin(gamma);
    dydt(4) = v*cos(gamma)/(R_E+h);
        
end

function [F,G] = ObjFun(x_vect,param_2st,param_3st)
% Objective function to which is subject of a constrained minimization
    % INPUTS : 
    %   The states on which the fmincon will act to find the local minimum
    %   v = x_vect(1)
    %   gamma0 = x_vect(2)
    %   h0 = x_vect(3)
    
    % OUTPUT :
    %   F     scalar function of the deltaV losses
    %   G     analytic formula of the gradient of the function

    % Extract variables
    v0 = x_vect(1);     % [m/s]
    gamma0 = x_vect(2); % [rad]
    h0 = x_vect(3);     % [m]
    delta_t = x_vect(4); %[t]
    omega0 = 0;         % [rad]
    y0  = [v0;gamma0;h0;omega0];  

    t_burn2 = param_2st.t_burn ;
    t_burn3 =param_3st.t_burn;

    % Propagating the trajectory
    options = odeset('reltol', 1e-5, 'abstol', 1e-6);
    [tt2, xx2] =      ode78(@(t,x) fun(t,x,param_2st), [0 t_burn2], y0,options);
    [tt_off, xx_off] =    ode78(@(t,x) fun_off(t,x,param_2st), [0 delta_t], xx2(end,:),options);
    [tt3, xx3] =    ode78(@(t,x) fun(t,x,param_3st), [0 t_burn3], xx_off(end,:),options);
    
    v_vect2 =       xx2(:,1);
    gamma_vect2 =   xx2(:,2);
    h_vect2 =       xx2(:,3);

    v_vect_off = xx_off(:,1);
    gamma_vect_off = xx_off(:,2);
    h_vect_off = xx_off(:,3);

    v_vect3 =       xx3(:,1);
    gamma_vect3 =   xx3(:,2);
    h_vect3 =       xx3(:,3);

    % Manouvre intensity with sign
    Cb2 = param_2st.m0/(param_2st.A*param_2st.Cd);
    Cb_off = (param_2st.m0-param_2st.m_prop)/(param_2st.A*param_2st.Cd);
    Cb3 = param_3st.m0/(param_3st.A*param_3st.Cd);
    H_scale = 7500; %[m]
    R_E = 6.37813660e+06; %[m]
    
    rho =   @(h) 1.225*exp(-h./H_scale);
    m2 =    @(t) param_2st.m0 - param_2st.T./(param_2st.I_sp*9.81).*t;
    m3 =    @(t) param_3st.m0 - param_3st.T./(param_3st.I_sp*9.81).*t;
    g =     @(h) 9.81/(1+h./R_E).^2;

    term2_d = -0.5./Cb2.*rho(h_vect2).*v_vect2.^2.*param_2st.m0./m2(tt2);
    term_off_d = -0.5./Cb_off.*rho(h_vect_off).*v_vect_off.^2.;
    term3_d = -0.5./Cb3.*rho(h_vect3).*v_vect3.^2.*param_3st.m0./m3(tt3);

    term2_g = -g(h_vect2)'.*sin(gamma_vect2);
    term_off_g = -g(h_vect_off)'.*sin(gamma_vect_off);
    term3_g = -g(h_vect3)'.*sin(gamma_vect3);

    delta_v_drag = abs(trapz(tt2,term2_d)) + abs(trapz(tt_off,term_off_d) + abs(trapz(tt3,term3_d)));
    delta_v_grav = abs(trapz(tt2,term2_g)) + abs(trapz(tt_off,term_off_g) + abs(trapz(tt3,term3_g)));
    
    F =  delta_v_drag + delta_v_grav;

        if nargout>1
            G = [];
        end 

end

function [c,ceq,G,Geq] = myconstraints(x_vect,param_2st,param_3st)
% Non linear contraint function
    
    % Extract variables
    v0 = x_vect(1);     % [m/s]
    gamma0 = x_vect(2); % [rad]
    h0 = x_vect(3);     % [m]
    delta_t = x_vect(4); %[s]
    omega0 = 0;         % [rad]
    y0  = [v0;gamma0;h0;omega0];

    t_burn2 = param_2st.t_burn ;
    t_burn3 =param_3st.t_burn;

    mu = 398600.435; % gravitational constant of the Earth [km^3/s^2]
    R_E = 6.37813660e+06;   % [m]
    h_f = 400e3;            % [m] target final altitude
    r_f = R_E+h_f;          % [m] final target radius

    % Propagating the trajectory
    options = odeset('reltol', 1e-5, 'abstol', 1e-6);
    [~, xx2] =      ode78(@(t,x) fun(t,x,param_2st), [0 t_burn2], y0,options);
    [~, xx_off] =    ode78(@(t,x) fun_off(t,x,param_2st), [0 delta_t], xx2(end,:),options);
    [~, xx3] =    ode78(@(t,x) fun(t,x,param_3st), [0 t_burn3], xx_off(end,:),options);

    xx_f = xx3(end,:);
    v_f = xx_f(1);
    gamma_f =xx_f(2);
    h_f = xx_f(3);

    c = [];

    ceq(1) = sqrt(mu/(r_f/1e3))*1e3 - v_f;
    ceq(2) = gamma_f-0;
    ceq(3) = h_f - 400e3;

    if nargout>2
        G = [];
        Geq =[];
    end

end

function [F,G] = ObjFun_pitchup(x_vect,param_1st)
% Objective function to which is subject of a constrained minimization
    % INPUTS : 
    %   The states on which the fmincon will act to find the local minimum
    %   v = x_vect(1)
    %   gamma0 = x_vect(2)
    %   h0 = x_vect(3)
    
    % OUTPUT :
    %   F     scalar function of the deltaV losses
    %   G     analytic formula of the gradient of the function

    % Extract variables
    v0 = x_vect(1);     % [m/s]
    gamma0 = x_vect(2); % [rad]
    h0 = x_vect(3);     % [m]
    t_bo = x_vect(4);

    omega0 = 0;         % [rad]
    y0  = [v0;gamma0;h0;omega0];  

    %t_burn1 = param_1st.t_burn ;

    % Matching condition
    v_target = param_1st.v_f;
    gamma_target = param_1st.gamma_f;
    h_target = param_1st.h_f;

    % Propagating the trajectory
    options = odeset('reltol', 1e-5, 'abstol', 1e-6);
    [~, xx1] =      ode78(@(t,x) fun_pitchup(t,x,param_1st,t_bo), [0 t_bo], y0,options);
    
    v_f =       xx1(end,1);
    gamma_f =   xx1(end,2);
    h_f =       xx1(end,3);

    err_v = abs(v_f - v_target);
    err_gamma = abs(gamma_f - gamma_target);
    err_h = abs(h_f - h_target);
    
    F =  err_v + err_gamma + err_h;

        if nargout>1
            G = [];
        end 

end

function [c,ceq,G,Geq] = myconstraints_pitchup(x_vect,param_1st)
% Non linear contraint function
    
    % Extract variables
    v0 = x_vect(1);     % [m/s]
    gamma0 = x_vect(2); % [rad]
    h0 = x_vect(3);     % [m]
    t_bo = x_vect(4);

    omega0 = 0;         % [rad]
    y0  = [v0;gamma0;h0;omega0];

    %t_burn1 = param_1st.t_burn ;

    % Matching condition
    v_target = param_1st.v_f;
    gamma_target = param_1st.gamma_f;
    h_target = param_1st.h_f;

    % Propagating the trajectory
    options = odeset('reltol', 1e-5, 'abstol', 1e-6);
    [~, xx1] =      ode78(@(t,x) fun_pitchup(t,x,param_1st,t_bo), [0 t_bo], y0,options);

    xx_f = xx1(end,:);
    v_f = xx_f(1);
    gamma_f =xx_f(2);
    h_f = xx_f(3);

    c = [];

    ceq(1) = v_f - v_target;
    ceq(2) = gamma_f - gamma_target;
    ceq(3) = h_f - h_target;

    if nargout>2
        G = [];
        Geq =[];
    end

end

