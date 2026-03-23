%% OPTIMAL STAGING ANALYSIS

clear
close all
clc

% ----------------------------3ST0 SOLID-----------------------------------
% Recovering the specific problem data
g0   = 9.80665;   

% conservative estimate on isp
isp1 = 250;     
isp2 = 250;     
isp3 = 250;     
eps1 = 0.08;
eps2 = 0.11;
eps3 = 0.14;
% from mission requirements
mpay = 400;
% see MA Staging Propulsion rationale document on  OneDrive for reasoning
% behind the first DV budget estimation:
% ITERATA 1: DV LOSS = 1.3 KM/S DV TOT = 8.92e3;
DV   = 8.92e3-0.6e3; %m/s

% computation of the effective exaust velocities
c1 = isp1*g0;
c2 = isp2*g0;
c3 = isp3*g0;

% Solving for lambda:=l
f = @(l) -(DV-c1*log((c1*l-1)/(eps1*c1*l))-c2*log((c2*l-1)/(eps2*c2*l))...
    -c3*log((c3*l-1)/(eps3*c3*l)));

x_vect = linspace(-0.2,0.2,100000);
f_eval = zeros(length(x_vect),1);
for j = 1:length(x_vect)
    f_eval(j) = f(x_vect(j));
end

figure
plot(x_vect,f_eval,'o',MarkerSize=1);
grid on
yline(0)
title('COST FUNCTION f(l) BEHAVIOUR','Interpreter','latex',FontSize=12);
xlabel('l',FontSize=15,Interpreter='latex');
ylabel('fval',FontSize=15,Interpreter='latex');

% solution
options = optimoptions("fsolve","Display",'off');
[l,~,exitflag] =fsolve(f,1,options);

% recovering n1,n2,n3
% defining eps and c vectors
eps = [eps1; eps2; eps3];
c   = [c1; c2; c3];
n = zeros(3,1);
for i=1:3
    n(i) = (c(i)*l-1)/(eps(i)*c(i)*l);
end

% once n1, n2, n3 have been computed, it's possible to recover the stage
% masses according to the following equations:

m3 = ((n(3)-1)/(1-n(3)*eps3))*mpay;
m2 = ((n(2)-1)/(1-n(2)*eps2))*(mpay+m3);
m1 = ((n(1)-1)/(1-n(1)*eps1))*(mpay+m2+m3);
m = [m1;m2;m3];

% GROSS LIFT OFF MASS 
glom = m1+m2+m3+mpay; % [Kg]

% Knowing the mass of the stages and the structural mass indexces
% m_structure  = m_stage*eps_stage
% m_propellant = m_stage*(1-eps_stage)

m_struct = zeros(3,1);
m_prop   = zeros(3,1);
for i=1:3
    m_struct(i) = m(i)*eps(i);
    m_prop(i) = m(i)*(1-eps(i));
end

clc
fprintf('\n3ST0 SOLID ANALYSIS')
fprintf('\n gross liftoff mass = %f kg \n',glom );
fprintf('\n first stage mass = %f Kg ',m1);
fprintf('\n second stage mass = %f kg',m2);
fprintf('\n third stage mass = %f kg \n',m3);
fprintf('\n first stage propellant mass = %f kg',m_prop(1));
fprintf('\n second stage propellant mass = %f kg',m_prop(2));
fprintf('\n third stage propellant mass = %f kg \n',m_prop(3));


%%  Sensitiviy of GLOM with respect to 15% variation of eps1, esp2, esp3

% Given the nominal values of the structural indexes, the GLOM computation
% is now perform considering a 10-15% change in the structural mass indexes 

eps1_plus  = linspace(eps1,0.15*eps1+eps1,11)';
eps1_minus = linspace(-0.15*eps1+eps1,eps1,10)';
Eps_1 = [eps1_minus;eps1_plus(2:end)];

eps2_plus  = linspace(eps2,0.15*eps2+eps2,11)';
eps2_minus = linspace(-0.15*eps2+eps2,eps2,10)';
Eps_2 = [eps2_minus;eps2_plus(2:end)];

eps3_plus  = linspace(eps3,0.15*eps3+eps3,11)';
eps3_minus = linspace(-0.15*eps3+eps3,eps3,10)';
Eps_3 = [eps3_minus;eps3_plus(2:end)];

% Computation of the GLOM for all the possible combinations of variability 
% of the structural mass indexces:

GLOM = zeros(20,20,20);
for i=1:20
    for j= 1:20
        for k = 1:20
            f = @(l) -(DV-c1*log((c1*l-1)/(Eps_1(i)*c1*l))-c2*log((c2*l-1)/...
                (Eps_2(j)*c2*l))-c3*log((c3*l-1)/(Eps_3(k)*c3*l)));
            l =fsolve(f,1,options);
            
            n1 = (c1*l-1)/(c1*Eps_1(i)*l);
            n2 = (c2*l-1)/(c2*Eps_2(j)*l);
            n3 = (c3*l-1)/(c3*Eps_3(k)*l);

            m3 = ((n3-1)/(1-n3*Eps_3(k)))*mpay;
            m2 = ((n2-1)/(1-n2*Eps_2(j)))*(mpay+m3);
            m1 = ((n1-1)/(1-n1*Eps_1(i)))*(mpay+m2+m3);
            GLOM(i,j,k) = m1+m2+m3+mpay;
        end
    end
end

%% Scatter plot of the variability region
[x, y, z] = ndgrid(Eps_1, Eps_2, Eps_3);
x = x(:); 
y = y(:);
z = z(:);
glom_flat = GLOM(:);

figure;
scatter3(x, y, z, 20, glom_flat, 'filled');
xlabel('\epsilon_1',FontSize=15);
ylabel('\epsilon_2',FontSize=15);
zlabel('\epsilon_3',FontSize=15);
title('3STO solid GLOM expected variability','Interpreter','latex',...
    FontSize=15);
colorbar;
clim([min(glom_flat) max(glom_flat)]);
hold on
scatter3(eps1,eps2,eps3,glom,'filled',MarkerFaceColor='red',SizeData=30);

%% SENSITIVITY ANALYSIS

% Extraxting the submatrix GLOM(1:20;10,10) to see how the variability of
% the first stack structural mass index influence the GLOM with respect to
% the nominal case.

figure
% INFLUENCE OF THE FIRST STRUCTURAL MASS INDEX, WITH CONSTANT EPS_2, EPS_3
TradeOffMatrix1 = GLOM(1:20,10,10);
plot(((Eps_1-Eps_1(10).*ones(20,1))./Eps_1(10))*100,((TradeOffMatrix1-...
    glom*ones(20,1))./glom)*100,'o',MarkerFaceColor='blue')
grid on
hold on
xlabel('Change in structural mass index [\%] ',...
   'Interpreter','latex', 'FontSize',15);
ylabel('Change in GLOM [\%]', 'Interpreter','latex',FontSize=15);
title('Optimization strategy analysis', 'Interpreter','latex',FontSize=15);

% INFLUENCE OF THE SECOND STRUCTURAL MASS INDEX, WITH CONSTANT EPS_1, EPS_3
TradeOffMatrix2 = GLOM(10,1:20,10)';
plot(((Eps_2-Eps_2(10).*ones(20,1))./Eps_2(10))*100,((TradeOffMatrix2-...
    glom*ones(20,1))./glom)*100,'o',MarkerFaceColor='green')

% INFLUENCE OF THE THIRD STRUCTURAL MASS INDEX, WITH CONSTANT EPS_1, EPS_2
TradeOffMatrix3 = squeeze(GLOM(10,10,1:20));
plot(((Eps_3-Eps_3(10).*ones(20,1))./Eps_3(10))*100,((TradeOffMatrix3-...
    glom*ones(20,1))./glom)*100,'o',MarkerFaceColor='yellow')

% reference point
plot((Eps_1(10)-Eps_1(10))./Eps_1(10),(TradeOffMatrix1(10)-glom)/...
    TradeOffMatrix1(10),'o',MarkerFaceColor='red')

legend('eps1 variable','eps2 variable','eps3 variable','ref','FontSize',...
    14,Interpreter = 'latex');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOME CONSIDERATIONS: IT'S MUCH BETTER OPTIMIZING THE FIRST STAGE
% STRUCTURAL MASS INDEXCES, AS FOR THE SAME % CHANGE, THE OVERALL EFFECT ON
% THE GLOM IS MUCH STRONGER WITH RESPECT TO THE OPTIMIZATION OF THE
% STRUCTURAL MASS INDEX OF THE UPPER STAGES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%-----------------------------2STO LIQUID----------------------------------
% reliable assumptions for eps_i and isp_i, i=1,2 have to be made.
% [1] Conceptual Launch Vehicle and Spacecract Design for Risk Assessment
eps1 = 0.15;
eps2 = 0.17;

% [2] Atlas V and Falcon 9 both exploit LOX-RP1 for the first stage
isp1 = 320;
isp2 = 320;  

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

m2 = ((n(2)-1)/(1-n(2)*eps2))*mpay;
m1 = ((n(1)-1)/(1-n(1)*eps1))*(mpay+m2);

% GROSS LIFT OFF MASS 
glom_TSTO = m1+m2+mpay; % [Kg]

m_struct_TSTO = zeros(2,1);
m_prop_TSTO   = zeros(2,1);
m = [m1;m2];
for i=1:2
    m_struct_TSTO(i) = m(i)*eps(i);
    m_prop_TSTO(i) = m(i)*(1-eps(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The mission seems possible from the energetic point of view, but some
% constraints may be present on the acceleration profile, thus a TSTO may
% not be a feasible option, depending on the loads that the structure
% and/or the crew can soupport.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n2STO LIQUID')
fprintf('\n gross liftoff mass = %f kg \n',glom_TSTO );
fprintf('\n first stage mass = %f Kg ',m1);
fprintf('\n second stage mass = %f kg \n',m2);
fprintf('\n first stage propellant mass = %f kg',m_prop_TSTO(1));
fprintf('\n second stage propellant mass = %f kg \n',m_prop_TSTO(2));

% Sensitivity analysis
% Given the nominal values of the structural indexes, the GLOM computation
% is now perform considering a 10-15% change in the structural mass indexes 

eps1_plus  = linspace(eps1,0.15*eps1+eps1,11)';
eps1_minus = linspace(-0.15*eps1+eps1,eps1,10)';
Eps_1 = [eps1_minus;eps1_plus(2:end)];

eps2_plus  = linspace(eps2,0.15*eps2+eps2,11)';
eps2_minus = linspace(-0.15*eps2+eps2,eps2,10)';
Eps_2 = [eps2_minus;eps2_plus(2:end)];

% Computation of the GLOM for all the possible combinations of variability 
% of the structural mass indexces:

GLOM = zeros(20,20);
for i=1:20
    for j= 1:20
            f = @(l) DV-c1*log((c1*l-1)/(Eps_1(i)*c1*l))-c2*log((c2*l-1)...
                /(Eps_2(j)*c2*l));
            l =fsolve(f,1,options);
            
            n1 = (c1*l-1)/(c1*Eps_1(i)*l);
            n2 = (c2*l-1)/(c2*Eps_2(j)*l);

            m2 = ((n2-1)/(1-n2*Eps_2(j)))*(mpay);
            m1 = ((n1-1)/(1-n1*Eps_1(i)))*(mpay+m2);
            GLOM(i,j) = m1+m2+mpay;     
    end
end

% Scatter plot of the variability region
[x, y] = ndgrid(Eps_1, Eps_2);
x = x(:); 
y = y(:);
glom_flat = GLOM(:);

figure;
scatter(x, y, 20, glom_flat, 'filled');
xlabel('\epsilon_1',FontSize=15);
ylabel('\epsilon_2',FontSize=15);
title('2STO liquid GLOM expected variability','Interpreter','latex',...
    FontSize=15);
colorbar;
clim([min(glom_flat) max(glom_flat)]);
hold on
scatter3(eps1,eps2,eps3,glom,'filled',MarkerFaceColor='red',SizeData=30);
%%
%------------------------------3STO LIQUID---------------------------------
DV = 8.92e3-0.6e3;
% conservative estimate on isp
isp1 = 320; % LOX-RP1     
isp2 = 320; % LOX-RP1    
isp3 = 320; % LOX-RP1    
% to be refine
eps1 = 0.15;
eps2 = 0.16;
eps3 = 0.17;

% computation of the effective exaust velocities
c1 = isp1*g0;
c2 = isp2*g0;
c3 = isp3*g0;

% Solving for lambda:=l
f = @(l) -(DV-c1*log((c1*l-1)/(eps1*c1*l))-c2*log((c2*l-1)/(eps2*c2*l))...
    -c3*log((c3*l-1)/(eps3*c3*l)));

% solution
options = optimoptions("fsolve","Display",'off');
[l,~,exitflag] =fsolve(f,1,options);

% recovering n1,n2,n3
% defining eps and c vectors
eps = [eps1; eps2; eps3];
c   = [c1; c2; c3];
n = zeros(3,1);
for i=1:3
    n(i) = (c(i)*l-1)/(eps(i)*c(i)*l);
end

% once n1, n2, n3 have been computed, it's possible to recover the stage
% masses according to the following equations:

m3 = ((n(3)-1)/(1-n(3)*eps3))*mpay;
m2 = ((n(2)-1)/(1-n(2)*eps2))*(mpay+m3);
m1 = ((n(1)-1)/(1-n(1)*eps1))*(mpay+m2+m3);
m = [m1;m2;m3];

% GROSS LIFT OFF MASS 
glom = m1+m2+m3+mpay; % [Kg]

% Knowing the mass of the stages and the structural mass indexces
% m_structure  = m_stage*eps_stage
% m_propellant = m_stage*(1-eps_stage)

m_struct = zeros(3,1);
m_prop   = zeros(3,1);
for i=1:3
    m_struct(i) = m(i)*eps(i);
    m_prop(i) = m(i)*(1-eps(i));
end

fprintf('\n3ST0 LIQUID ANALYSIS')
fprintf('\n gross liftoff mass = %f kg \n',glom );
fprintf('\n first stage mass = %f Kg ',m1);
fprintf('\n second stage mass = %f kg',m2);
fprintf('\n third stage mass = %f kg \n',m3);
fprintf('\n first stage propellant mass = %f kg',m_prop(1));
fprintf('\n second stage propellant mass = %f kg',m_prop(2));
fprintf('\n third stage propellant mass = %f kg \n',m_prop(3));

%--------------------CHOICE OF PRESSURIZING THECNOLOGY---------------------
% First choice of specific impulse for LOX-RP1 230s
% Second iteration for specific impulse: from CEA code analysis 230s with
% O/F = 2.13. From literature: rho_rp1 = 820 kg/m^3; rho_lox = 1140 kg/m^3

OF = 2.13;
rho_rp1 = 820;
rho_lox = 1140;
% FIRST STAGE
m_rp1_1 = m_prop(1)/(1+OF);
m_lox_1 = m_prop(1)-m_rp1_1;
v_rp1_1 = m_rp1_1/rho_rp1;
v_lox_1 = m_lox_1/rho_lox;
% SECOND STAGE
m_rp1_2 = m_prop(2)/(1+OF);
m_lox_2 = m_prop(2)-m_rp1_2;
v_rp1_2 = m_rp1_2/rho_rp1;
v_lox_2 = m_lox_2/rho_lox;
% THIRD STAGE
m_rp1_3 = m_prop(3)/(1+OF);
m_lox_3 = m_prop(3)-m_rp1_3;
v_rp1_3 = m_rp1_3/rho_rp1;
v_lox_3 = m_lox_3/rho_lox;

% CHOICE OF TANK CONFIGURATION: [slide 11. 06-01]
% tandem with common bulk head, internal piping, seems the most convinient
% May be critical a common head for semi criogenic couples
% fuel tank always above the oxidizer tank 

% tank shape: from literature
% cylindrical body + elliptical emispheroids dome --> Expected tank axpect
% ratio AR>1

% tank cyclindrical diameter ~= launcher dyameter
diameter_launcher = 1.6; % ext
diameter_tanks = diameter_launcher -0.2; % inner 

% tank length computation

% dome volume computation
AR = 2; % Aspect ratio, to be confirmed by the literature
vol_domes = 2*(2/3*pi*(diameter_tanks/2)^3/AR);
% first stage length of cylindrical part
Lenghts.one.rp1 = (v_rp1_1-vol_domes)/(pi*diameter_tanks^2/4); 
Lenghts.one.lox = (v_lox_1)/(pi*diameter_tanks^2/4);
% second stage lenght of cyndrical part
Lenghts.two.rp1 = (v_rp1_2-vol_domes)/(pi*diameter_tanks^2/4); 
Lenghts.two.lox = (v_lox_2)/(pi*diameter_tanks^2/4);

% third stage lenght --> a cylinder will be too short --> sort of sphere
% shape, with the borders touching the external walls of the launcher
% third stage domes
diameter_three = 1.2; % to be updated
AR_three = 4;
vol_domes_three = 2*(2/3*pi*(diameter_three/2)^3/AR_three);
Lenghts.three.rp1 = (v_rp1_3-vol_domes_three)/(pi*diameter_three^2/4); 
Lenghts.three.lox = (v_lox_3-vol_domes_three)/(pi*diameter_three^2/4);

h_dome_three = (diameter_three/2)/AR_three;
h_dome_one = (diameter_launcher/2)/AR;

% old sphere assumption: diameter is too small
% Lenghts.three.rp1 = ((v_rp1_3)*3/(4*pi))^(1/3); 
% Lenghts.three.lox = ((v_lox_3)*3/(4*pi))^(1/3); 

% tanks total lenght:
Tank.one.length = Lenghts.one.lox+Lenghts.one.rp1+2*(diameter_tanks/2)/AR;
Tank.two.length = Lenghts.two.lox+Lenghts.two.rp1+2*(diameter_tanks/2)/AR;
Tank.three.length = Lenghts.three.lox*2+Lenghts.three.rp1*2;

% Launcher total lenght estimate:  [slide 10 06-01]
L_TOT = diameter_launcher+h_dome_three...
    +1/3*diameter_three+Lenghts.three.lox*2+...
    1/4*diameter_three+Lenghts.three.rp1*2+1/3*diameter_three+...
    2*diameter_three+ h_dome_three+2*h_dome_three+Lenghts.three.rp1+...
    Lenghts.three.lox+h_dome_one+diameter_launcher+1/3*diameter_launcher+...
    diameter_launcher+Lenghts.one.rp1+Lenghts.one.lox+diameter_launcher

L_STACK2 = diameter_launcher+h_dome_three...
    +1/3*diameter_three+Lenghts.three.lox*2+...
    1/4*diameter_three+Lenghts.three.rp1*2+1/3*diameter_three+...
    2*diameter_three+ h_dome_three+2*h_dome_three+Lenghts.two.rp1+...
    Lenghts.two.lox+h_dome_one+diameter_launcher+1/3*diameter_launcher+...
    diameter_launcher-h_dome_one

% The lenght of the first stack is too big with the tank diameter constant.
% As first guess: diameter_three = 1.2 m
% check feasibility, lenght includes payload fairing ~=2.5 m

L_STACK3 = diameter_launcher+h_dome_three...
    +1/3*diameter_three+Lenghts.three.lox*2+...
    1/4*diameter_three+Lenghts.three.rp1*2+1/3*diameter_three+...
    2*diameter_three+ h_dome_three+2*h_dome_three-h_dome_one

% fairing shape may be too thin --> potrebbe doversi allargare per
% accomodare il payload, per poi richiudersi sulla forma del naso.

% ESTIMATION OF SINGLE COMPONENTS MASSES
MASS.wiring = 1.058*L_TOT^0.25*sqrt(glom);
MASS.Avionics = 10*glom^0.361;
% MASS FAIRING == REQUIRES EXTERNAL AREA OF THE FAIRING -> REQUIRE SHAPE

% THUST STRUCUTRE, LRE CORRELATIONS REQUIRES THE DEFINITION OF THE THRUST
% LEVEL, THUS THE DEFINITION OF THE T/W RATIO FOR ALL THE STAGES
T_W1 =  1;
T_W2 = 0.8;
T_W3 = 0.8;
% from graphical correlations
T1 = T_W1*g0*glom;
T2 = T_W2*g0*(mpay+m2+m3);
T3 = T_W3*g0*(mpay+m3);
% thrust structure weight
MASS.Tstruct1 = 2.55*T1*1e-4;
MASS.Tstruct2 = 2.55*T2*1e-4;
MASS.Tstruct3 = 2.55*T3*1e-4;
% eps nozzle definition to exploit empirical correlations
exp1 = 32.58;
exp2 = 40;
exp3 = 40;
% LRE mass estimation
MASS.LRE1 = 7.81*1e-4*T1+3.37*1e-5*T1*sqrt(exp1)+59;
MASS.LRE2 = 7.81*1e-4*T2+3.37*1e-5*T2*sqrt(exp1)+59;
MASS.LRE3 = 7.81*1e-4*T3+3.37*1e-5*T3*sqrt(exp1)+59;
% Cryogenic tanks insulation for LOX
% required computation of external tank area [slide 14 06-01]
% considered the tank as cylinder+EH dome
% eccentricity of the dome
E = sqrt(1-h_dome_one^2/(diameter_launcher/2)^2);
R_launher = diameter_launcher/2;
% is it ln or log?
A_outer_EHdome = pi*R_launher^2*(1+1/(2*E*AR^2)*log((1+E)/(1-E)));
A_outer_cyl_1lox = pi*diameter_launcher*Lenghts.one.lox;  
AOUT.one.loxtank = A_outer_EHdome*2+A_outer_cyl_1lox;
A_outer_cyl_2lox = pi*diameter_launcher*Lenghts.two.lox; 
AOUT.two.loxtank = A_outer_EHdome*2+A_outer_cyl_2lox;
% computation for the third stage
E_three = sqrt(1-h_dome_three^2/(diameter_three/2)^2);
R3 =diameter_three/2;
A_outer_EHdome3 = pi*R3^2*(1+1/(2*E*AR_three^2)*log((1+E_three)/(1-E_three)));
A_outer_cyl_3lox = pi*diameter_launcher*Lenghts.three.lox; 
AOUT.three.loxtank = A_outer_EHdome3*2+A_outer_cyl_3lox;

% Expected cryogenic insulation mass for LOX tanks
MASS.insulation1 = 1.12*AOUT.one.loxtank;
MASS.insulation2 = 1.12*AOUT.two.loxtank;
MASS.insulation3 = 1.12*AOUT.three.loxtank;

% Tanks mass computations --> burst pressur sizing

close all
%%  Sensitiviy of GLOM with respect to 15% variation of eps1, esp2, esp3

% Given the nominal values of the structural indexes, the GLOM computation
% is now perform considering a 10-15% change in the structural mass indexes 

eps1_plus  = linspace(eps1,0.15*eps1+eps1,11)';
eps1_minus = linspace(-0.15*eps1+eps1,eps1,10)';
Eps_1 = [eps1_minus;eps1_plus(2:end)];

eps2_plus  = linspace(eps2,0.15*eps2+eps2,11)';
eps2_minus = linspace(-0.15*eps2+eps2,eps2,10)';
Eps_2 = [eps2_minus;eps2_plus(2:end)];

eps3_plus  = linspace(eps3,0.15*eps3+eps3,11)';
eps3_minus = linspace(-0.15*eps3+eps3,eps3,10)';
Eps_3 = [eps3_minus;eps3_plus(2:end)];

% Computation of the GLOM for all the possible combinations of variability 
% of the structural mass indexces:

GLOM = zeros(20,20,20);
for i=1:20
    for j= 1:20
        for k = 1:20
            f = @(l) -(DV-c1*log((c1*l-1)/(Eps_1(i)*c1*l))-c2*log((c2*l-1)/...
                (Eps_2(j)*c2*l))-c3*log((c3*l-1)/(Eps_3(k)*c3*l)));
            l =fsolve(f,1,options);
            
            n1 = (c1*l-1)/(c1*Eps_1(i)*l);
            n2 = (c2*l-1)/(c2*Eps_2(j)*l);
            n3 = (c3*l-1)/(c3*Eps_3(k)*l);

            m3 = ((n3-1)/(1-n3*Eps_3(k)))*mpay;
            m2 = ((n2-1)/(1-n2*Eps_2(j)))*(mpay+m3);
            m1 = ((n1-1)/(1-n1*Eps_1(i)))*(mpay+m2+m3);
            GLOM(i,j,k) = m1+m2+m3+mpay;
        end
    end
end

%% Scatter plot of the variability region
[x, y, z] = ndgrid(Eps_1, Eps_2, Eps_3);
x = x(:); 
y = y(:);
z = z(:);
glom_flat = GLOM(:);

figure;
scatter3(x, y, z, 20, glom_flat, 'filled');
xlabel('\epsilon_1',FontSize=15);
ylabel('\epsilon_2',FontSize=15);
zlabel('\epsilon_3',FontSize=15);
title('3STO liquid GLOM expected variability','Interpreter','latex',...
    FontSize=15);
colorbar;
clim([min(glom_flat) max(glom_flat)]);
hold on
scatter3(eps1,eps2,eps3,glom,'filled',MarkerFaceColor='red',SizeData=30);









