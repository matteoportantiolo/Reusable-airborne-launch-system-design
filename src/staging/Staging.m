clear
close all
clc
format longG
g0   = 9.80665;   
mpay = 400;
DV = 9.17e3;

% reliable assumptions for eps_i and isp_i, i=1,2 have to be made.
% [1] Conceptual Launch Vehicle and Spacecract Design
% RP1-LOX couple literature search + margin to take into account for
% a smaller system wrt the ones in the literature.
eps1 = 0.12;
eps2 = 0.14;

% From CEA analysis
isp1 = 323.16;
isp2 = 336.58; % higher isp in the void with a nozzle with an higher eps  

% Optimization process for a TSTO:
% computation of the effective exaust velocities
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

MASS.m2 = ((n(2)-1)/(1-n(2)*eps2))*mpay;
MASS.m1 = ((n(1)-1)/(1-n(1)*eps1))*(mpay+MASS.m2);

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

% saving computed masses into a struct
MASS.mprop1 = m_prop(1);
MASS.mprop2 = m_prop(2);
MASS.mstru1 = m_struct(1);
MASS.mstru2 = m_struct(2);

% clearing the workspace
clear DV eps eps1 eps2 g0 i l c  m_prop m_struct

%---------------------PROPELLANT TANK SIZING PROCEDURE---------------------

OF = 2.13;
rho_rp1 = 820;
rho_lox = 1140;

% FIRST STAGE NOMINAL
MASS.m_rp1_1 = MASS.mprop1/(1+OF);
MASS.m_lox_1 = MASS.mprop1-MASS.m_rp1_1;
MASS.m_rp1_1 = MASS.m_rp1_1*1; % (no margin considered)
MASS.m_lox_1 = MASS.m_lox_1*1; % (no margin considered)

% FIRST STAGE TAKING INTO ACCOUNT 4 TONS FOR PUSHBACK
MASS.PushBack = 4000; % Kg
LOXPushBack = 4000/(1+OF);
RP1PushBack = 4000-LOXPushBack;
MASS.m_rp1_1 = MASS.m_rp1_1+RP1PushBack;
MASS.m_lox_1 = MASS.m_lox_1+LOXPushBack;


VOLUME.v_rp1_1 = (MASS.m_rp1_1/rho_rp1)*1.03; % (3% margin considered)
VOLUME.v_lox_1 = (MASS.m_lox_1/rho_lox)*1.03; % (3% margin considered)

% SECOND STAGE 
MASS.m_rp1_2 = MASS.mprop2/(1+OF);
MASS.m_lox_2 = MASS.mprop2-MASS.m_rp1_2;
MASS.m_rp1_2 = MASS.m_rp1_2*1; % (no margin considered)
MASS.m_lox_2 = MASS.m_lox_2*1; % (no margin considered)
VOLUME.v_rp1_2 = (MASS.m_rp1_2/rho_rp1)*1.03; % (3% margin considered)
VOLUME.v_lox_2 = (MASS.m_lox_2/rho_lox)*1.03;

clear  isp1 isp2 m n

%-------------------CONFIGURATION AND LENGTH ESTIMATION--------------------
% Tanks shape: cylindrical body + elliptical emispheroids dome
% Tanks configuration: tandem tanks without common bulkhead
% Denser propellant below: LOX at the bottom, RP1 at the top 

% Launcher diameters
D1 = 1.6;
D2 = 1.3;
AR = 2; % aspect ratio of the EH
h_dome = D1/(2*AR);
h_dome2 = D2/(2*AR);

% first stage computations
VOLUME.dome = (2/3*pi*(D1/2)^3/AR);
% RP1 tank length
LENGTH.h_dome = h_dome;
LENGTH.rp1_cyl_1 = (VOLUME.v_rp1_1-2*VOLUME.dome)/(pi*D1^2/4);
LENGTH.rp1_1  = LENGTH.rp1_cyl_1+2*h_dome;
% LOX tank length
LENGTH.lox_cyl_1 = (VOLUME.v_lox_1-2*VOLUME.dome)/(pi*D1^2/4);
LENGTH.lox_1 = LENGTH.lox_cyl_1+2*h_dome;
%-------------------------------------------------------------------------
% second stage computations
VOLUME.dome2 = (2/3*pi*(D2/2)^3/AR);
% RP1 tank length
LENGTH.h_dome2 = h_dome2;
LENGTH.rp1_cyl_2 = (VOLUME.v_rp1_2-2*VOLUME.dome2)/(pi*D2^2/4);
LENGTH.rp1_2  = LENGTH.rp1_cyl_2+2*h_dome2;
% LOX tank length
LENGTH.lox_cyl_2 = (VOLUME.v_lox_2-2*VOLUME.dome2)/(pi*D2^2/4);
LENGTH.lox_2 = LENGTH.lox_cyl_2+2*h_dome2;

% first guess of the system length from Edberg and Costa, 2020
% skirt support strucutre is considered in corrispondance of the EH domes
%--------------------------------------------------------------------------
% thrust structure+first step after skirt ~= D1
% intertank length ~= 1/4*D1+2*h_dome
% interstage length ~= D1
% secnd skirt length ~= 1/3*D2+h_dome2
% intertank length ~= 1/4*D2+2*h_dome2
% foward skirt length ~= 1/3*D2+h_dome2
% payload fairing length ~= 2*D2
%--------------------------------------------------------------------------

LENGTH.first_skirt = D1;
LENGTH.Intertank1 = 1/4*D1+2*h_dome;
LENGTH.Interstage1 = D1;
LENGTH.second_skirt = 1/3*D2+h_dome2;
LENGTH.Intertank2 = 1/4*D2+2*h_dome2;
LENGTH.Foward_skirt = 1/3*D2+h_dome2;
LENGTH.Fairing = 2*D2;

LENGTH.system = LENGTH.first_skirt+LENGTH.lox_cyl_1+LENGTH.Intertank1+...
    LENGTH.rp1_1+LENGTH.Interstage1+LENGTH.second_skirt+LENGTH.lox_cyl_2+...
    LENGTH.Intertank2+LENGTH.rp1_cyl_2+LENGTH.Foward_skirt+LENGTH.Fairing;

%----------SYSTEM MASS ESTIMATION FROM MER AND PRELIMINARY DESIGN----------

TW1 = 1.6;
TW2 = 1.4;
exp1 = 30;
exp2 = 99.056;
A_fairing = 10; % from aerodynamics

% thust structures masses
MASS.ThrustStruct1 = 2.55*1e-4*TW1*(MASS.glom+4000);
MASS.ThrustStruct2 = 2.55*1e-4*TW2*(MASS.m2+mpay);
% LRE weight from generic correlations including nozzle weigth
MASS.LRE1 = 7.81*1E-4*TW1*(MASS.glom+4000)+3.37*1e-5*TW1*(MASS.glom+4000)*sqrt(exp1)+59;
MASS.LRE2 = 7.81*1E-4*TW2*(MASS.m2+mpay)+3.37*1e-5*TW2*(MASS.m2+mpay)*...
    sqrt(exp2)+59;
% avionics, fairing and wiring mass estimations
MASS.avionics = 10*(MASS.glom+4000)^0.361;
MASS.fairing = 4.95*A_fairing^(1.15);
MASS.wiring = 1.058*LENGTH.system^(0.25)*sqrt(MASS.glom+40000);

%-------------mass estimation of the first stage thurbopumps---------------
% [1] from space propulsion slides (i.e. from Sutton book)
P_tank = 2; % bar
P_cc = 150; % bar, from propulsion, eventually to be updated
P_inj_losses = 0.05*P_cc; % bar
P_feed_losses = 0.5; % bar
N = 30000; % pump rotational speed in RPM, rough guess
eta = 0.7; % pump efficiency

Head_bar = P_cc-P_tank+P_inj_losses+P_feed_losses;
% recovering the mass flow rate of LOX and RP1
T1 = TW1*MASS.glom;
m_flow_propellant = T1/c1;
m_flow_RP1 = (1/(1+OF))*m_flow_propellant;
m_flow_LOX = m_flow_propellant-m_flow_RP1;
% recovering volumetric flow rates for LOX and RP1
Q_LOX = m_flow_LOX/rho_lox;
Q_RP1 = m_flow_RP1/rho_rp1;
% Required power
Power_LOX = Q_LOX*Head_bar*100000/eta;
Power_RP1 = Q_RP1*Head_bar*100000/eta;
% Mass estimation
A = 1.95;  % medium value between the expected range of variability
b = 0.634; % medium value between the expected range of variability
MASS.M_TPLOX = A*(Power_LOX/(N*2*pi/60))^b;
MASS.M_TPRP1 = A*(Power_RP1/(N*2*pi/60))^b;
%--------------------------------------------------------------------------
% Cryogenic tank insulation for LOX

% LOX tank surface area computation
% first stage
% Eccentricity of the EH dome
E1 = sqrt(1-h_dome^2/(D1/2)^2);
% is it ln or log?
A_outer_EHdome = pi*(D1/2)^2*(1+1/(2*E1*AR^2)*log((1+E1)/(1-E1)));
A_outer_cyl_1lox = pi*D1*LENGTH.lox_cyl_1;  
A_LOX1 = A_outer_EHdome*2+A_outer_cyl_1lox;
% second stage
% Eccentricity of the EH dome
E2 = sqrt(1-h_dome^2/(D2/2)^2);
% is it ln or log?
A_outer_EHdome2 = pi*(D2/2)^2*(1+1/(2*E2*AR^2)*log((1+E2)/(1-E2)));
A_outer_cyl_2lox = pi*D2*LENGTH.lox_cyl_2; 
A_LOX2 = A_outer_EHdome2*2+A_outer_cyl_2lox;

MASS.LOXinsulation1 = 1.12*A_LOX1;
MASS.LOXinsulation2 = 1.12*A_LOX2;
%--------------------------------------------------------------------------
% Pressurizing system preliminary design for LOX and RP1 tank
% choice of Helium, gamma =1.6, allow for lower temperature at the end of
% the isoentropic expansion, thus it's preferred with respect to Helium
% since we are dealing with LOX

% the aim is to maintained the propellant tanks at 2 bar
T_pg_in = 273.15; % K, initial temperature of the pressurizing gas
P_pg_in = 500*100000; % pa, initial pressure of the press. gas (500 bar)
P_pg_end = (1.75)*100000; % pa, final pressure of the press.gas >=1.75bar

% final temperature of the pressurizing gas computation
% PH: isoentropic expansion from p.g. tank to propellant tank
gamma = 1.66;
T_pg_end = T_pg_in*(P_pg_end/P_pg_in)^((gamma-1)/gamma);

% PG TANK FOR LOX
V_pg_in_lox = (VOLUME.v_lox_1*(P_pg_end/P_pg_in)^(1/gamma))/...
    (1-(P_pg_end/P_pg_in)^(1/gamma));
M_pg_lox = (P_pg_in*V_pg_in_lox)/(T_pg_in*8314/4);
% soupposing five spherical helium spheres to be placed under the LOX tank
% in between the thrust structure
vol = V_pg_in_lox/6;
R_spheres_lox = (vol*3/(4*pi))^(1/3);

% PG TANK FOR RP1
V_pg_in_rp1 = (VOLUME.v_rp1_1*(P_pg_end/P_pg_in)^(1/gamma))/...
    (1-(P_pg_end/P_pg_in)^(1/gamma));
M_pg_rp1 = (P_pg_in*V_pg_in_rp1)/(T_pg_in*8314/4);
% soupposing three spherical helium spheres, to be placed above the RP1
% tank or in between the two tanks
vol = V_pg_in_rp1/6;
R_spheres_rp1 = (vol*3/(4*pi))^(1/3);

MASS.PressurizingGas = M_pg_rp1+M_pg_lox;

% Pressurizing spheres mass estimation
% MATERIAL: Ti-6Al-4V
% ultimate yield = 1170 Mpa
sigma = 1.170; % Gpa
% Density = 4430 Kg/m^3
rho = 4430;
% required thikness for LOX spheres
t_lox = (P_pg_in*2*R_spheres_lox)/(2*sigma*1e9);
% mass of a single tank
MASS.lox_spheres = 5*t_lox*4*pi*R_spheres_lox^2*rho;
% required thikness for RP1 spheres
t_rp1 = (P_pg_in*2*R_spheres_rp1)/(2*sigma*1e9);
% mass of a single tank
MASS.rp1_spheres = 5*t_rp1*4*pi*R_spheres_rp1^2*rho;
%--------------------------------------------------------------------------
% LOX and RP1 tank thikness and mass computation to widstand 2 bar of
% pressure

% required thikness
% cylindrical tank 
% considering 8 g ofmax ammissible acceleration
t_LOX_CYL = ((1.75+rho_lox*8*9.81*LENGTH.lox_1)*2*D1/2)/(sigma*1e9);
% EH dome, first approximation as spherical shape
t_LOX_DOME = (1.75*2*D1/2)/(2*sigma*1e9);
% Selected thikness (common to both tanks)
t = max([1,t_LOX_DOME,t_LOX_CYL]);
t = t*1e-3; % from mm to m

% Recovering the tanks mass for the first stage
% lox
m_cylinder_lox = t*rho*A_outer_cyl_1lox;
m_domes_lox = t*rho*A_outer_EHdome*2;
MASS.LOX_tank = m_cylinder_lox+m_domes_lox;
% rp1
m_cylinder_rp1 = t*rho*pi*D1*LENGTH.rp1_cyl_1;
m_domes_rp1 = t*rho*A_outer_EHdome*2;
MASS.RP1_tank = m_cylinder_rp1+m_domes_rp1;
%--------------------------------------------------------------------------

%---SECOND STAGE TP, PRESSURIZING GAS AND PROPELLANT TANKS COMPUTATIONS----

% recovering the mass flow rate of LOX and RP1 for the second stage
T2 = TW2*(MASS.m2+mpay);
m_flow_propellant = T2/c2;
m_flow_RP1 = (1/(1+OF))*m_flow_propellant;
m_flow_LOX = m_flow_propellant-m_flow_RP1;
% recovering volumetric flow rates for LOX and RP1
Q_LOX = m_flow_LOX/rho_lox;
Q_RP1 = m_flow_RP1/rho_rp1;
% Required power
Power_LOX = Q_LOX*Head_bar*100000/eta;
Power_RP1 = Q_RP1*Head_bar*100000/eta;
% Mass estimation
A = 1.95;  % medium value between the expected range of variability
b = 0.634; % medium value between the expected range of variability
MASS.M_TPLOX_stg2 = A*(Power_LOX/(N*2*pi/60))^b;
MASS.M_TPRP1_stg2 = A*(Power_RP1/(N*2*pi/60))^b;
%--------------------------------------------------------------------------
% Pressurizing gas estimation
% PG TANK FOR LOX
V_pg_in_lox = (VOLUME.v_lox_2*(P_pg_end/P_pg_in)^(1/gamma))/...
    (1-(P_pg_end/P_pg_in)^(1/gamma));
M_pg_lox = (P_pg_in*V_pg_in_lox)/(T_pg_in*8314/4);
% soupposing five spherical helium spheres to be placed under the LOX tank
% in between the thrust structure
vol = V_pg_in_lox/3;
R_spheres_lox = (vol*3/(4*pi))^(1/3)

% PG TANK FOR RP1
V_pg_in_rp1 = (VOLUME.v_rp1_2*(P_pg_end/P_pg_in)^(1/gamma))/...
    (1-(P_pg_end/P_pg_in)^(1/gamma));
M_pg_rp1 = (P_pg_in*V_pg_in_rp1)/(T_pg_in*8314/4);
% soupposing three spherical helium spheres, to be placed above the RP1
% tank or in between the two tanks
vol = V_pg_in_rp1/3;
R_spheres_rp1 = (vol*3/(4*pi))^(1/3)
MASS.PressurizingGas_2stg = M_pg_rp1+M_pg_lox;

% required thikness for LOX spheres
t_lox = (P_pg_in*2*R_spheres_lox)/(2*sigma*1e9);
% mass of a single tank
MASS.lox_spheres_2stg = 5*t_lox*4*pi*R_spheres_lox^2*rho;
% required thikness for RP1 spheres
t_rp1 = (P_pg_in*2*R_spheres_rp1)/(2*sigma*1e9);
% mass of a single tank
MASS.rp1_spheres_2stg = 5*t_rp1*4*pi*R_spheres_rp1^2*rho;
%--------------------------------------------------------------------------
% LOX and RP1 tank thikness and mass computation to widstand 2 bar of
% pressure

% required thikness
% cylindrical tank 
t_LOX_CYL = ((1.75+rho_lox*LENGTH.lox_cyl_2)*2*D2/2)/(sigma*1e9);
% EH dome, first approximation as spherical shape
t_LOX_DOME = (1.75*2*D2/2)/(2*sigma*1e9);
% Selected thikness (common to both tanks)
t = max([1,t_LOX_DOME,t_LOX_CYL]);
t = t*1e-3; % from mm to m

% Recovering the tanks mass for the first stage
% lox
m_cylinder_lox = t*rho*A_outer_cyl_2lox;
m_domes_lox = t*rho*A_outer_EHdome2*2;
MASS.LOX_tank_2stg = m_cylinder_lox+m_domes_lox;
% rp1
m_cylinder_rp1 = t*rho*pi*D2*LENGTH.rp1_cyl_2;
m_domes_rp1 = t*rho*A_outer_EHdome2*2;
MASS.RP1_tank_2stg = m_cylinder_rp1+m_domes_rp1;
%--------------------------------------------------------------------------
% CASE OF A PRESSURE FED SYSTEM
% pressurizant needed
% Pressurizing gas estimation
% Final pressure of the pg == final pressure of the tanks == pcc+losses
P_pg_end = (P_cc+P_inj_losses+P_feed_losses)*100000;

% PG TANK FOR LOX
V_pg_in_lox = (VOLUME.v_lox_2*(P_pg_end/P_pg_in)^(1/gamma))/...
    (1-(P_pg_end/P_pg_in)^(1/gamma));
MASS.M_pg_lox_PF = (P_pg_in*V_pg_in_lox)/(T_pg_in*8314/4);

% PG TANK FOR RP1
V_pg_in_rp1 = (VOLUME.v_rp1_2*(P_pg_end/P_pg_in)^(1/gamma))/...
    (1-(P_pg_end/P_pg_in)^(1/gamma));
MASS.M_pg_rp1_PF = (P_pg_in*V_pg_in_rp1)/(T_pg_in*8314/4);
%--------------------------------------------------------------------------


