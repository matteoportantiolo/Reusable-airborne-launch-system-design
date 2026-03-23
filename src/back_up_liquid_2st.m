clear
clc
close all
format long
plotStyle


% TOTAL OPTIMIZATION WITH ENGINES ON/OFF (2 STAGES)


% 2STO LIQUID
% glom	23959	kg
% m1	21016	kg
% mp1	18494	kg
% ms1	2522	kg
% m2	2543	kg
% mp2	2187	kg
% ms2	356  	kg


% structural epsilon
% 1st --> 0.12
% 2st --> 0.14


% TRAJECTORY PHASES PEGASUS XL
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


%% INITIALIZATION

% play with TW and masses

t_i = 0;

% First stage
param_1st.D = 1.6;          % [m]  
param_1st.A = pi*(param_1st.D/2)^2;   % [m^2]
param_1st.Cd = 0.5;         % [-]
param_1st.Cl = 0;           % [-]
param_1st.m0 = 23959;    % [kg]
param_1st.I_sp = 323.16;       % [s] 
param_1st.delta = 0;        % [rad]
param_1st.m_prop = 18494;
param_1st.n = 1/(-param_1st.m_prop/param_1st.m0+1);
param_1st.TW = 1.6;
param_1st.T = param_1st.m0*9.81*param_1st.TW;  % [N] %T/W=1

% Second stage
param_2st.D = 1.3;          % [m] 
param_2st.A = pi*(param_2st.D/2)^2;   % [m^2]
param_2st.Cd = 0;         % [-]
param_2st.Cl = 0;           % [-]
param_2st.m0 = 2543+400;    % [kg]
param_2st.I_sp = 336.58;       % [s] 
param_2st.delta = 0;        % [rad]
param_2st.m_prop = 2187;
param_2st.n = 1/(-param_2st.m_prop/param_2st.m0+1);
param_2st.TW = 1.4;
param_2st.T = param_2st.m0*9.81*param_2st.TW;  % [N] %T/W=1

t_burn1 = param_1st.m_prop*param_1st.I_sp*9.81/param_1st.T;
param_1st.t_burn = t_burn1;

t_burn2 = param_2st.m_prop*param_2st.I_sp*9.81/param_2st.T;
param_2st.t_burn = t_burn2;

param_1st.gamma_rate = deg2rad(30)/param_1st.t_burn;

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

%% Carrier velocity

M_c = 0.82;
h_c = 11.88e3;
gradient = 0.0065; 
T_amb = 15+273.15;
T_c = T_amb - h_c*gradient;
gamma = 1.4;
R = 287;
a_c = sqrt(gamma*R*T_c);
v_c = M_c*a_c;

%% Optimization process

v0 = 220;  % [m/s]
gamma0 = deg2rad(1); % [rad]
omega0 =     0; % [m]
h0 =     11e3; % [m]

y0 = [v0;gamma0;h0;omega0];

t_bo1 = 150;
delta_t0 = 200;
t_bo2 = 170;

y0_new = [y0(1:3);t_bo1;delta_t0;t_bo2];

options = optimoptions("fmincon","Algorithm","sqp",'Display','iter-detailed',...
                      'SpecifyObjectiveGradient',false,'SpecifyConstraintGradient',false);

%options.StepTolerance = 1.000000e-10;
options.ConstraintTolerance = 1e-10;
options.MaxFunctionEvaluations = 5e3;
options.MaxIterations = 800;
options.StepTolerance = 1e-08;
lb = [ 200, 0.001, 9e3, 0, 0, 0];
ub = [ v_c+20, deg2rad(2), 12e3, t_burn1, 500, t_burn2]; 

[x_opt,f] = fmincon(@(x)ObjFun(x,param_1st,param_2st), y0_new,...
                                  [],[],[],[],lb,ub, @(x)myconstraints(x,param_1st,param_2st),...
                                  options);

delta_v_opt = f;

fprintf(' \n The optimized delta_v loss is : %+.5e [m/s] \n',delta_v_opt)

xx_0_NEW  = [x_opt(1:3);0];
t_bo1_opt = x_opt(4);
delta_t_opt = x_opt(5);
t_bo2_opt = x_opt(6);

fprintf([' \n The optimized initial state is :\n \n initial velocity : %+.5e [m/s] \n initial pitch angle : %+.5e [deg]'  ...
    '\n initial height : %+.5e [km] \n burnout time1 : %+.5e [s] \n delta_t off : %+.5e [s] \n burnout time2 : %+.5e [s] \n'], ...
    xx_0_NEW(1),rad2deg(xx_0_NEW(2)),xx_0_NEW(3)./1e3,t_bo1_opt,delta_t_opt,t_bo2_opt)

options_pu = odeset('reltol', 1e-5, 'abstol', 1e-6,'Events',@(t,y) Karman_line(t,y,true));
options = odeset('reltol', 1e-5, 'abstol', 1e-6);
[tt1, xx1] =      ode78(@(t,x) fun_pitchup(t,x,param_1st), [0 t_bo1_opt], xx_0_NEW,options_pu);
[tt_off, xx_off] =      ode78(@(t,x) fun_off(t,x,param_2st), [0 delta_t_opt], xx1(end,:),options);
[tt2, xx2] =      ode78(@(t,x) fun(t,x,param_2st), [0 t_bo2_opt], xx_off(end,:),options);

mu = 398600.435; 
R_E = 6.37813660e+06;  
h_f = 400e3;          
r_f = R_E+h_f;       
v_insertion = sqrt(mu/(r_f/1e3))*1e3;
err = abs(xx2(end,1) - v_insertion);
fprintf(' \n The error in orbit insertion velocity is : %+.5e [m/s] \n',err)

fprintf([' \n The states at separation are: \n velocity : %+.5e [m/s] \n flight path angle : %+.5e [deg] ' ...
    '\n altitude : %+.5e [m] \n'],xx1(end,1),rad2deg(xx1(end,2)),xx1(end,3))

if t_bo1_opt > tt1(end)
    t_bo1_opt_pre = t_bo1_opt;
    fprintf(' \n First stage burn out due to altitude target : delta_t = %+.5e [s] \n',(t_bo1_opt-tt1(end)))
    t_bo1_opt = tt1(end); %burn out due to event function (target h)
end

%% Singular plots

figure(2)
subplot(2,2,1)
plot(tt1,xx1(:,1),'LineWidth',2)
hold on
plot(tt_off + tt1(end),xx_off(:,1),'LineWidth',2)
plot(tt2 + tt_off(end) + tt1(end),xx2(:,1),'LineWidth',2)
xlabel('t [s]','Interpreter','latex',FontSize=14)
ylabel('v [m/s]','Interpreter','latex',FontSize=14)
title("Variation of velocity in time ",'Interpreter','latex')
legend(["First stage pitch-up","Separation and coasting","Second stage"],'Location','best',FontSize=14)
grid on

subplot(2,2,2)
plot(tt1,rad2deg(xx1(:,2)),'LineWidth',2)
hold on
plot(tt_off + tt1(end),rad2deg(xx_off(:,2)),'LineWidth',2)
plot(tt2 + tt_off(end) + tt1(end),rad2deg(xx2(:,2)),'LineWidth',2)
xlabel('t [s]','Interpreter','latex',FontSize=14)
ylabel('$\gamma$ [deg]','Interpreter','latex',FontSize=14)
title("Variation of $\gamma$ in time ",'Interpreter','latex')
legend(["First stage pitch-up","Separation and coasting","Second stage"],'Location','best',FontSize=14)
grid on

subplot(2,2,3)
plot(tt1,xx1(:,3),'LineWidth',2)
hold on
plot(tt_off + tt1(end),xx_off(:,3),'LineWidth',2)
plot(tt2 + tt_off(end) + tt1(end),xx2(:,3),'LineWidth',2)
xlabel('t [s]','Interpreter','latex',FontSize=14)
ylabel('h [km]','Interpreter','latex',FontSize=14)
title("Altitude as function of time ",'Interpreter','latex')
legend(["First stage pitch-up","Separation and coasting","Second stage"],'Location','best',FontSize=14)
grid on

subplot(2,2,4)
plot(tt1,xx1(:,4),'LineWidth',2)
hold on
plot(tt_off + tt1(end),xx_off(:,4),'LineWidth',2)
plot(tt2 + tt_off(end) + tt1(end),xx2(:,4),'LineWidth',2)
xlabel('t [s]','Interpreter','latex',FontSize=14)
ylabel('$\omega$ [deg]','Interpreter','latex',FontSize=14)
title("Angular position as function of time ",'Interpreter','latex')
legend(["First stage pitch-up","Separation and coasting","Second stage"],'Location','best',FontSize=14)
grid on

%% Total plot

tt_update = [tt1;tt_off + tt1(end);tt2 + tt_off(end) + tt1(end)];
xx_update = [xx1;xx_off;xx2];

figure
yyaxis left
plot(tt_update,xx_update(:,1),'LineWidth',2)
ylabel('v [m/s]','Interpreter','latex',FontSize=14)
yyaxis right
plot(tt_update,xx_update(:,2),'LineWidth',2)
ylabel('$\gamma$ [deg]','Interpreter','latex',FontSize=14)
xlabel('t [s]','Interpreter','latex',FontSize=14)
xline(tt1(end),'--','LineWidth',2)
title("Variations in time ",'Interpreter','latex')
legend(["Velocity","Flight path angle"],'Location','best',FontSize=14)
grid on

% % Aggiungi il terzo asse per la terza curva
% ax1 = gca; % Cattura l'asse corrente
% % Crea un nuovo asse nella stessa posizione per la terza curva
% ax2 = axes('Position', ax1.Position, ...
%            'XAxisLocation', 'top', ...
%            'YAxisLocation', 'right', ...
%            'Color', 'none'); % Rendi trasparente il nuovo asse
% % Plotta la terza curva sul nuovo asse
% hold(ax2, 'on');  % Mantiene il grafico corrente
% plot(tt_update,xx_update(:,3),'LineWidth',2)  % Disegna la terza curva
% ylabel('h [km]','Interpreter','latex',FontSize=14) % Etichetta per il terzo asse y
% % Regola la posizione dell'asse aggiuntivo leggermente a destra
% % ax2.YAxis(1).Position(1) = 1.15; % Sposta leggermente il terzo asse y a destra

%% Delta_v losses

v_vect1 = xx1(:,1);
gamma_vect1 = xx1(:,2);
h_vect1 = xx1(:,3);

v_vect_off = xx_off(:,1);
gamma_vect_off = xx_off(:,2);
h_vect_off = xx_off(:,3);

v_vect2 = xx2(:,1);
gamma_vect2 = xx2(:,2);
h_vect2 = xx2(:,3);

% Manouvre intensity with sign
Cb1 = param_1st.m0/(param_1st.A*param_1st.Cd);
Cb_off = param_2st.m0/(param_2st.A*param_2st.Cd);
Cb2 = param_2st.m0/(param_2st.A*param_2st.Cd);

H_scale = 7500; %[m]
R_E = 6.37813660e+06; %[m]

rho =   @(h) 1.225*exp(-h./H_scale);
m1 =    @(t) param_1st.m0 - param_1st.T./(param_1st.I_sp*9.81).*t;
m2 =    @(t) param_2st.m0 - param_2st.T./(param_2st.I_sp*9.81).*t;
g =     @(h) 9.81./((1+h./R_E).^2);

term1_d = -0.5./Cb1.*rho(h_vect1).*v_vect1.^2.*param_1st.m0./m1(tt1);
term_off_d = -0.5./Cb_off.*rho(h_vect_off).*v_vect_off.^2;
term2_d = -0.5./Cb2.*rho(h_vect2).*v_vect2.^2.*param_2st.m0./m2(tt2);

term1_g = -g(h_vect1).*sin(gamma_vect1);
term_off_g = -g(h_vect_off).*sin(gamma_vect_off);
term2_g = -g(h_vect2).*sin(gamma_vect2);

delta_v_drag = abs(trapz(tt1,term1_d)) + abs(trapz(tt_off,term_off_d)) + abs(trapz(tt2,term2_d));
delta_v_grav = abs(trapz(tt1,term1_g)) + abs(trapz(tt_off,term_off_g)) + abs(trapz(tt2,term2_g));

delta_v_tot = delta_v_drag + delta_v_grav;

fprintf(' \n The total delta_v losses are : %+.5e [m/s] \n',delta_v_tot)

%% Masses left

t1_left = t_burn1 - t_bo1_opt;
t2_left = t_burn2 - t_bo2_opt;

m1_left = param_1st.T./(param_1st.I_sp*9.81).*t1_left;
m2_left = param_2st.T./(param_2st.I_sp*9.81).*t2_left;

perc1 = m1_left/param_1st.m_prop*100;
perc2 = m2_left/param_2st.m_prop*100;

fprintf(' \n Propellant masses letf are : \n stage 1 : %+.5e [kg], percentage : %+.5e \n stage 2 : %+.5e [kg], percentage : %+.5e \n',m1_left,perc1,m2_left,perc2)

%% Accelerations

D1 = @(rho,v) 0.5.*rho.*(v.^2).*param_1st.A.*param_1st.Cd;
D_off = @(rho,v) 0.5.*rho.*(v.^2).*param_2st.A.*param_2st.Cd;
D2 = @(rho,v) 0.5.*rho.*(v.^2).*param_2st.A.*param_2st.Cd;

% to do in body axis
% a1x = (param_1st.T./m1(tt1) - D1(rho(h_vect1),xx1(:,1))./m1(tt1) - g(h_vect1).*sin(xx1(:,2)))./9.81.*cos(xx1(:,2));
% a1y = (param_1st.T./m1(tt1) - D1(rho(h_vect1),xx1(:,1))./m1(tt1) - g(h_vect1).*sin(xx1(:,2)))./9.81.*sin(xx1(:,2));
% 
% a_sepx = (D_sep(rho(h_vect_sep),xx_sep(:,1))./(param_1st.m0-param_1st.m_prop) - g(h_vect_sep).*sin(xx_sep(:,2)))./9.81.*cos(xx_sep(:,2));
% a_sepy = (D_sep(rho(h_vect_sep),xx_sep(:,1))./(param_1st.m0-param_1st.m_prop) - g(h_vect_sep).*sin(xx_sep(:,2)))./9.81.*sin(xx_sep(:,2));
% 
% a2x = (param_2st.T./m2(tt2) - D2(rho(h_vect2),xx2(:,1))./m2(tt2) - g(h_vect2).*sin(xx2(:,2)))./9.81.*cos(xx2(:,2));
% a2y = (param_2st.T./m2(tt2) - D2(rho(h_vect2),xx2(:,1))./m2(tt2) - g(h_vect2).*sin(xx2(:,2)))./9.81.*sin(xx2(:,2));
% 
% a_offx = (D_off(rho(h_vect_off),xx_off(:,1))./(param_2st.m0-param_2st.m_prop) - g(h_vect_off).*sin(xx_off(:,2)))./9.81.*cos(xx_off(:,2));
% a_offy = (D_off(rho(h_vect_off),xx_off(:,1))./(param_2st.m0-param_2st.m_prop) - g(h_vect_off).*sin(xx_off(:,2)))./9.81.*sin(xx_off(:,2));
% 
% a3x = (param_3st.T./m3(tt3) - D3(rho(h_vect3),xx3(:,1))./m3(tt3) - g(h_vect3).*sin(xx3(:,2)))./9.81.*cos(xx3(:,2));
% a3y = (param_3st.T./m3(tt3) - D3(rho(h_vect3),xx3(:,1))./m3(tt3) - g(h_vect3).*sin(xx3(:,2)))./9.81.*sin(xx3(:,2));

a1 = (param_1st.T./m1(tt1) - D1(rho(h_vect1),xx1(:,1))./m1(tt1) - g(h_vect1).*sin(xx1(:,2)))./9.81;
a_off = (D_off(rho(h_vect_off),xx_off(:,1))./param_2st.m0 - g(h_vect_off).*sin(xx_off(:,2)))./9.81;
a2 = (param_2st.T./m2(tt2) - D2(rho(h_vect2),xx2(:,1))./m2(tt2) - g(h_vect2).*sin(xx2(:,2)))./9.81;

a_tot = [a1;a_off;a2];

figure
plot(tt_update,a_tot)
xlabel('t [s]','Interpreter','latex',FontSize=14)
ylabel('a [$m/(s^2)$]','Interpreter','latex',FontSize=14)
title("Evolution of acceleration",'Interpreter','latex')
grid on

%% MAX Q 

% PegasusXL --> maxq = 6,923e4 [Pa]

max_q1 = 0.5.*rho(h_vect1).*v_vect1.^2;
max_q_off = 0.5.*rho(h_vect_off).*v_vect_off.^2;
max_q2 = 0.5.*rho(h_vect2).*v_vect2.^2;

figure
plot(tt1,max_q1,'LineWidth',2)
xlabel('t [s]','Interpreter','latex',FontSize=14)
ylabel('q [Pa]','Interpreter','latex',FontSize=14)
title("MAX Q during pitch up",'Interpreter','latex')
grid on

figure
plot(tt1,max_q1,'LineWidth',2)
hold on
plot(tt_off + tt1(end),max_q_off,'LineWidth',2)
plot(tt2 + tt_off(end) + tt1(end),max_q2,'LineWidth',2)
xlabel('t [s]','Interpreter','latex',FontSize=14)
ylabel('q [Pa]','Interpreter','latex',FontSize=14)
title("MAX Q",'Interpreter','latex')
legend(["First stage pitch-up","Separation and coasting","Second stage"],'Location','best',FontSize=14)
grid on

%% ---------------------------FUNCTIONS--------------------------------

function plotStyle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set figure properties for better looking plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% interpreter:
set(0, 'defaultTextInterpreter', 'Latex')
set(0, 'defaultAxesTickLabelInterpreter', 'Latex')
set(0, 'defaultLegendInterpreter', 'Latex')
set(0,'defaultAxesXGrid','on')
set(0,'defaultAxesYGrid','on')
% lines:
set(0,'defaultLineLineWidth', 1.5);
set(0,'defaultLineMarkerSize',6) ;
set(0,'defaultLineMarkerEdgeColor','k')
set(0,'defaultLineMarkerFaceColor','auto')
% legend:
set(0, 'defaultLegendLocation','northoutside');
set(0, 'defaultLegendOrientation','horizontal');
set(0, 'defaultLegendFontSize',12);
% axes:
set(0,'defaultAxesFontSize',16);
end

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
        m = m0 - m_prop0;
    end

    dydt(1) = T/m*cos(delta) -D/m -g*sin(gamma);
    dydt(2) = v*cos(gamma)/(R_E+h) + T*sin(delta)/(m*v) + L/m -g*cos(gamma)/v ;
    dydt(3) = v*sin(gamma);
    dydt(4) = v*cos(gamma)/(R_E+h);
        
end

function [dydt] = fun_pitchup(t,y,param)

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
    gamma_rate = param.gamma_rate; % [rad]
    delta = param.delta; % [rad]
    
    H_scale = 7500; %[m]
    R_E = 6.37813660e+06; %[m]
    g = 9.81/(1+h/R_E)^2;
    rho = 1.225*exp(-h/H_scale);
    D = 0.5*rho*(v^2)*A*Cd;
    %L = 0.5*rho*(v^2)*A*Cl;
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
    % gamma_targ = param.gamma_targ; % [rad]
    %gamma_rate_coast = param.gamma_rate_coast; % [rad]
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
    t_bo1 = x_vect(4);
    delta_t = x_vect(5); %[t]
    t_bo2 = x_vect(6);

    omega0 = 0;         % [rad]
    y0  = [v0;gamma0;h0;omega0];  

    % t_burn1 = param_1st.t_burn;
    % t_burn2 = param_2st.t_burn;
    % t_burn3 =param_3st.t_burn;

    % Propagating the trajectory
    options_pu = odeset('reltol', 1e-5, 'abstol', 1e-6,'Events',@(t,y) Karman_line(t,y,true));
    options = odeset('reltol', 1e-5, 'abstol', 1e-6);
    [tt1, xx1] =      ode78(@(t,x) fun_pitchup(t,x,param_1st), [0 t_bo1], y0,options_pu);
    [tt_off, xx_off] =      ode78(@(t,x) fun_off(t,x,param_2st), [0 delta_t], xx1(end,:),options);
    [tt2, xx2] =      ode78(@(t,x) fun(t,x,param_2st), [0 t_bo2], xx_off(end,:),options);
    
    v_vect1 = xx1(:,1);
    gamma_vect1 = xx1(:,2);
    h_vect1 = xx1(:,3);
    
    v_vect_off = xx_off(:,1);
    gamma_vect_off = xx_off(:,2);
    h_vect_off = xx_off(:,3);

    v_vect2 = xx2(:,1);
    gamma_vect2 = xx2(:,2);
    h_vect2 = xx2(:,3);

    % Manouvre intensity with sign
    Cb1 = param_1st.m0/(param_1st.A*param_1st.Cd);
    Cb_off = param_2st.m0/(param_2st.A*param_2st.Cd);
    Cb2 = param_2st.m0/(param_2st.A*param_2st.Cd);

    H_scale = 7500; %[m]
    R_E = 6.37813660e+06; %[m]
    
    rho =   @(h) 1.225*exp(-h./H_scale);
    m1 =    @(t) param_1st.m0 - param_1st.T./(param_1st.I_sp*9.81).*t;
    m2 =    @(t) param_2st.m0 - param_2st.T./(param_2st.I_sp*9.81).*t;
    g =     @(h) 9.81./((1+h./R_E).^2);

    term1_d = -0.5./Cb1.*rho(h_vect1).*v_vect1.^2.*param_1st.m0./m1(tt1);
    term_off_d = -0.5./Cb_off.*rho(h_vect_off).*v_vect_off.^2;
    term2_d = -0.5./Cb2.*rho(h_vect2).*v_vect2.^2.*param_2st.m0./m2(tt2);

    term1_g = -g(h_vect1).*sin(gamma_vect1);
    term_off_g = -g(h_vect_off).*sin(gamma_vect_off);
    term2_g = -g(h_vect2).*sin(gamma_vect2);

    delta_v_drag = abs(trapz(tt1,term1_d)) + abs(trapz(tt_off,term_off_d)) + abs(trapz(tt2,term2_d));
    delta_v_grav = abs(trapz(tt1,term1_g)) + abs(trapz(tt_off,term_off_g)) + abs(trapz(tt2,term2_g));
    
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
    t_bo1 = x_vect(4);
    delta_t = x_vect(5); %[t]
    t_bo2 = x_vect(6);

    omega0 = 0;         % [rad]
    y0  = [v0;gamma0;h0;omega0];

    % t_burn1 = param_1st.t_burn;
    % t_burn2 = param_2st.t_burn;
    % t_burn3 =param_3st.t_burn;

    mu = 398600.435; % gravitational constant of the Earth [km^3/s^2]
    R_E = 6.37813660e+06;   % [m]
    h_f = 400e3;            % [m] target final altitude
    r_f = R_E+h_f;          % [m] final target radius

    % Propagating the trajectory
    options_pu = odeset('reltol', 1e-5, 'abstol', 1e-6,'Events',@(t,y) Karman_line(t,y,true));
    options = odeset('reltol', 1e-5, 'abstol', 1e-6);
    [~, xx1] =      ode78(@(t,x) fun_pitchup(t,x,param_1st), [0 t_bo1], y0,options_pu);
    [~, xx_off] =      ode78(@(t,x) fun_off(t,x,param_2st), [0 delta_t], xx1(end,:),options);
    [~, xx2] =      ode78(@(t,x) fun(t,x,param_2st), [0 t_bo2], xx_off(end,:),options);

    xx_f = xx2(end,:);
    v_f = xx_f(1);
    gamma_f =xx_f(2);
    h_f = xx_f(3);

    c = [];

    ceq(1) = sqrt(mu/(r_f/1e3))*1e3 - v_f;
    ceq(2) = gamma_f - 0;
    ceq(3) = h_f - 400e3;

    if nargout>2
        G = [];
        Geq =[];
    end

end

function [value, isterminal, direction] = Karman_line(~, y, isTerminal)
    h = y(3);
    value = h - 130e3;    % stop when h = 100e3
    isterminal = isTerminal;       % stop integration
    direction = 0;        % direction
end

