clear; 
clc;
close all;
format long;

% 3ST0 SOLID ANALYSIS
%  gross liftoff mass = 40957.552014 kg 
% 
%  first stage mass = 43233.481521 Kg 
%  second stage mass = 5742.267535 kg
%  third stage mass = 1076.158916 kg 
% 
%  first stage propellant mass = 31868.561101 kg
%  second stage propellant mass = 4390.060173 kg
%  third stage propellant mass = 847.237434 kg 
% 
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
%  first stage mass = 11343.207282 Kg 
%  second stage mass = 3240.259092 kg
%  third stage mass = 781.862449 kg 
% 
%  first stage propellant mass = 8085.765960 kg
%  second stage propellant mass = 2508.313347 kg
%  third stage propellant mass = 609.997334 kg

t_i = 0;

D = 1.5;                    % [m]
param_1st.T=  60e3;        % [N] T/W=1-->60, lict-->196
param_1st.A = pi*(D/2)^2;   % [m^2]
param_1st.Cd = 0.5;         % [-]
param_1st.Cl = 0;           % [-]
param_1st.m0 = 5742.267535+400;    % [kg]
param_1st.I_sp = 350;       % [s] 
param_1st.delta = 0;        % [rad]
param_1st.m_prop = 4390.060173;
param_1st.n = 1/(-param_1st.m_prop/param_1st.m0 + 1);

param_2st.T=  20e3;        % [N] T/W=1-->20, lict-->36
param_2st.A = pi*(D/2)^2;   % [m^2]
param_2st.Cd = 0.5;         % [-]
param_2st.Cl = 0;           % [-]
param_2st.m0 = 1076.158916+400;     % [kg]
param_2st.I_sp = 350;       % [s] 
param_2st.delta = 0;        % [rad]
param_2st.m_prop = 847.237434;
param_2st.n = 1/(-param_2st.m_prop/param_2st.m0 + 1);

t_burn = param_1st.m_prop*param_1st.I_sp*9.81/param_1st.T;
t_f = t_burn;
param_1st.t_burn = t_burn;
t_burn2 = param_2st.m_prop*param_2st.I_sp*9.81/param_2st.T;
param_2st.t_burn = t_burn2;

v0 = 2000; %400*3.6; % [m/s]
gamma0 = deg2rad(33); % [rad]
omega0 =     0; % [m]
h0 =     10e3; % [m]

y0 = [v0;gamma0;h0;omega0];

% %% Tring to compute the propagation backward in time
% t_i = 300;
% t_f = 0;
% R_E = 6.37813660e+06; %[m]
% v0_back = sqrt(398600./(R_E/1000+400))*1000; % [m/s]
% gamma0_back = 0.1;
% h0_back = 400e3;
% y0 = [v0_back;gamma0_back;h0_back;deg2rad(20)];
% param_1st.m0 = 400; 

% Perform integration
options_STM = odeset('reltol', 1e-6, 'abstol', 1e-8);


[tt, xx] = ode78(@(t,x) fun(t,x,param_1st), [t_i t_f], y0);
[tt2, xx2] = ode78(@(t,x) fun(t,x,param_2st), [0 t_burn2], xx(end,:));

tt_update = [tt;tt2+tt(end)];
xx_update = [xx;xx2];

figure(1);
subplot(2,2,1)
plot(tt_update,xx_update(:,1))
xlabel('t [s]','Interpreter','latex',FontSize=14)
ylabel('v [m/s]','Interpreter','latex',FontSize=14)
title("Variation of velocity in time ",'Interpreter','latex');
grid on;

figure(1);
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
ylabel('$\omega$ [km]','Interpreter','latex',FontSize=14)
title("Angular position as function of time ",'Interpreter','latex');
grid on;

%% Optimization process

delta_t0 = 1;
y0_new = [y0(1:3);delta_t0];

% options = optimoptions("fmincon",...
%     "Algorithm","interior-point",...
%     "EnableFeasibilityMode",true,...
%     "SubproblemAlgorithm","cg");

options = optimoptions("fmincon","Algorithm","sqp",'Display','iter-detailed',...
                      'SpecifyObjectiveGradient',false,'SpecifyConstraintGradient',false);

%options.StepTolerance = 1.000000e-10;
options.ConstraintTolerance = 1e-10;
options.MaxFunctionEvaluations = 5e3;
options.MaxIterations = 800;
options.StepTolerance = 1e-08;
lb = [ 200, 0.001, 7e3, 0];
ub = [ 2000, deg2rad(70), 12e3, 600]; % with 70deg it works perfectly

% Setting the lower and upper constraint on the variables on which fmincon
% will act to serch for the local minimum 

% y0 = 1e4.*[ 0.034668692417506;
%             0.000153546227838;
%             3.740329561789431;
%             0];

[x_min,f,exflag,output] = fmincon(@(x)ObjFun(x,param_1st,param_2st), y0_new,...
                                  [],[],[],[],lb,ub, @(x)myconstraints(x,param_1st,param_2st),...
                                  options);

xx_0_NEW  = [x_min(1:3);0];
delta_t_min = x_min(4);

fprintf([' \n The optimized initial state is :\n initial velocity : %+.5e [m/s] \n initial pitch angle : %+.5e [deg]'  ...
    '\n initial height : %+.5e [km] \n \n delta_t off : %+.5e [s] \n'],xx_0_NEW(1),rad2deg(xx_0_NEW(2)),xx_0_NEW(3)./1e3,delta_t_min)

[tt_bis, xx_bis] = ode78(@(t,x) fun(t,x,param_1st), [t_i t_f], xx_0_NEW);
[tt_off_bis, xx_off_bis] = ode78(@(t,x) fun_off(t,x,param_2st), [0 delta_t_min], xx_bis(end,:));
[tt2_bis, xx2_bis] = ode78(@(t,x) fun(t,x,param_2st), [0 t_burn2], xx_off_bis(end,:));

tt_update_bis = [tt_bis;tt_off_bis + tt_bis(end);tt2_bis + tt_off_bis(end)+tt_bis(end)];
xx_update_bis = [xx_bis;xx_off_bis;xx2_bis];

figure(2);
subplot(2,2,1)
plot(tt_update_bis,xx_update_bis(:,1))
xlabel('t [s]','Interpreter','latex',FontSize=14)
ylabel('v [m/s]','Interpreter','latex',FontSize=14)
title("Variation of velocity in time ",'Interpreter','latex');
grid on;

subplot(2,2,2)
plot(tt_update_bis,rad2deg(xx_update_bis(:,2)) )
xlabel('t [s]','Interpreter','latex',FontSize=14)
ylabel('$\gamma$ [deg]','Interpreter','latex',FontSize=14)
title("Variation of $\gamma$ in time ",'Interpreter','latex');
grid on;

subplot(2,2,3)
plot(tt_update_bis,xx_update_bis(:,3)./1000)
xlabel('t [s]','Interpreter','latex',FontSize=14)
ylabel('h [km]','Interpreter','latex',FontSize=14)
title("Altitude as function of time ",'Interpreter','latex');
grid on;

subplot(2,2,4)
plot(tt_update_bis,xx_update_bis(:,4))
xlabel('t [s]','Interpreter','latex',FontSize=14)
ylabel('$\omega$ [km]','Interpreter','latex',FontSize=14)
title("Angular position as function of time ",'Interpreter','latex');
grid on;

%% ---------------------------FUNCTIONS--------------------------------

function [dydt] = fun(t,y,param)
    v =     y(1); % [m/s]
    gamma = y(2); % [rad]
    h =     y(3); % [m]
    

    T=  param.T;
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

function [dydt] = fun_off(~,y,param)

    v =     y(1); % [m/s]
    gamma = y(2); % [rad]
    h =     y(3); % [m]
    
    A = param.A;
    Cd = param.Cd;
    Cl = param.Cl;
    m0 = param.m0;
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

    m = m0;
    T = 0;

    dydt(1) = T/m*cos(delta) -D/m -g*sin(gamma);
    dydt(2) = v*cos(gamma)/(R_E+h) + T*sin(delta)/(m*v) + L/m -g*cos(gamma)/v ;
    dydt(3) = v*sin(gamma);
    dydt(4) = v*cos(gamma)/(R_E+h);
        
    
end

function [F,G] = ObjFun(x_vect,param_1st,param_2st)
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
    t_f = param_1st.t_burn ;
    t_burn2 =param_2st.t_burn;
    % Propagating the trajectory
    options = odeset('reltol', 1e-5, 'abstol', 1e-6);
    [tt, xx] =      ode78(@(t,x) fun(t,x,param_1st), [0 t_f], y0,options);
    [tt_off, xx_off] =    ode78(@(t,x) fun_off(t,x,param_2st), [0 delta_t], xx(end,:),options);
    [tt2, xx2] =    ode78(@(t,x) fun(t,x,param_2st), [0 t_burn2], xx_off(end,:),options);

    % tt_update = [tt;tt2+tt(end)];
    % xx_update = [xx;xx2];
    
    v_vect1 =       xx(:,1);
    gamma_vect1 =   xx(:,2);
    h_vect1 =       xx(:,3);
    %omega_vect1 =   xx(:,4);

    v_vect_off = xx_off(:,1);
    gamma_vect_off = xx_off(:,2);
    h_vect_off = xx_off(:,3);

    v_vect2 =       xx2(:,1);
    gamma_vect2 =   xx2(:,2);
    h_vect2 =       xx2(:,3);
    %omega_vect2 =   xx2(:,4);

    % Manouvre intensity with sign
    Cb1 = param_1st.m0/(param_1st.A*param_1st.Cd);
    Cb2 = param_2st.m0/(param_2st.A*param_2st.Cd);
    Cb_off = Cb2;
    H_scale = 7500; %[m]
    R_E = 6.37813660e+06; %[m]
    
    rho =   @(h) 1.225*exp(-h./H_scale);
    m1 =    @(t) param_1st.m0 - param_1st.T./(param_1st.I_sp*9.81).*t;
    m2 =    @(t) param_2st.m0 - param_2st.T./(param_2st.I_sp*9.81).*t;
    g =     @(h) 9.81/(1+h./R_E).^2;
    term1_d = -0.5./Cb1.*rho(h_vect1).*v_vect1.^2.*param_1st.m0./m1(tt);
    term2_d = -0.5./Cb2.*rho(h_vect2).*v_vect2.^2.*param_2st.m0./m2(tt2);
    term_off_d = -0.5./Cb_off.*rho(h_vect_off).*v_vect_off.^2.;
    term1_g = -g(h_vect1)'.*sin(gamma_vect1);
    term2_g = -g(h_vect2)'.*sin(gamma_vect2);
    term_off_g = -g(h_vect_off)'.*sin(gamma_vect_off);
    delta_v_drag = abs(trapz(tt,term1_d)) + abs(trapz(tt2,term2_d)) + abs(trapz(tt_off,term_off_d));
    delta_v_grav = abs(trapz(tt,term1_g)) + abs(trapz(tt2,term2_g)) + abs(trapz(tt_off,term_off_g));
    
    F =  delta_v_drag + delta_v_grav;

        if nargout>1
            G = [];
        end 

end

function [c,ceq,G,Geq] = myconstraints(x_vect,param_1st,param_2st)
% Non linear contraint function
    
    % Extract variables
    v0 = x_vect(1);     % [m/s]
    gamma0 = x_vect(2); % [rad]
    h0 = x_vect(3);     % [m]
    delta_t = x_vect(4); %[s]
    omega0 = 0;         % [rad]
    y0  = [v0;gamma0;h0;omega0];
    mu = 398600.435; % gravitational constant of the Earth [km^3/s^2]
    R_E = 6.37813660e+06;   % [m]
    h_f = 400e3;            % [m] target final altitude
    r_f = R_E+h_f;          % [m] final target radius

    % Propagating the trajectory
    options = odeset('reltol', 1e-5, 'abstol', 1e-6);
    t_f = param_1st.t_burn ;
    t_burn2 =param_2st.t_burn;
    [~, xx] =      ode78(@(t,x) fun(t,x,param_1st), [0 t_f], y0,options);
    [~, xx_off] =    ode78(@(t,x) fun_off(t,x,param_2st), [0 delta_t], xx(end,:),options);
    [~, xx2] =    ode78(@(t,x) fun(t,x,param_2st), [0 t_burn2], xx_off(end,:),options);

    xx_f = xx2(end,:);
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

