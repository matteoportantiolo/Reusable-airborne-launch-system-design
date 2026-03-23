% Clear workspace and setup
clear; clc; close all;
format long;
plotStyle;

% ========================================================================
% OPTIMIZATION OF 2-STAGE ROCKET WITH ENGINE ON/OFF
% ========================================================================

% 2STO LIQUID
% Total mass (glom): 23959 kg
% Stage 1 mass (m1): 21016 kg
% Stage 1 propellant mass (mp1): 18494 kg
% Stage 1 structural mass (ms1): 2522 kg
% Stage 2 mass (m2): 2543 kg
% Stage 2 propellant mass (mp2): 2187 kg
% Stage 2 structural mass (ms2): 356 kg

% Structural efficiency (epsilon)
% Stage 1 (ε₁): 0.12
% Stage 2 (ε₂): 0.14

% TRAJECTORY PHASES FOR PEGASUS XL
%
% 1) **Release Phase**:
%    - Altitude (h) = 39,000 feet (~11,887 meters)
%    - Mach number (M) = 0.82
%    - Duration = 5 seconds
%
% 2) **Stage 1 Ignition**:
%    - Max thrust (Tmax) = 163,000 pounds (~726,000 N)
%    - Maximum altitude during Stage 1 = 53.9 km
%
% 3) **Stage 1 Separation**:
%    - Occurs after 16 seconds of flight
%
% 4) **Stage 2 Ignition**:
%    - Ignition occurs 1 second after Stage 1 separation
%    - Altitude during Stage 2 ignition = 72.1 km
%
% 5) **Fairing Separation**:
%    - Occurs at 116 km altitude
%    - Thrust during this phase = 44,000 pounds (~196,000 N)
%
% 6) **Coasting Phase**:
%    - Duration = 4 minutes (no thrust)
%
% 7) **Stage 2 Separation**:
%    - Occurs after the coasting phase
%
% 8) **Stage 3 Ignition**:
%    - Occurs 69 seconds after Stage 2 separation
%
% **Total Mission Duration**:
%    - Total flight time (t_tot) = 11 minutes 20 seconds

% ========================================================================
%% INITIALIZATION
% ========================================================================
t_start = 0; % Initial time

% Stage 1 parameters
stage1.D = 1.6;            % Diameter [m]
stage1.A = pi * (stage1.D / 2)^2; % Cross-sectional area [m^2]
stage1.Cd = 0.5;           % Drag coefficient
stage1.Cl = 0;             % Lift coefficient
stage1.m0 = 23959;         % Initial mass [kg]
stage1.I_sp = 323.16;      % Specific impulse [s]
stage1.m_prop = 18494;     % Propellant mass [kg]
stage1.TW = 1.8;           % Thrust-to-weight ratio
stage1.delta = 0;          % Thrust vector angle [rad]
stage1.T = stage1.m0 * 9.81 * stage1.TW; % Thrust [N]
stage1.t_burn = stage1.m_prop * stage1.I_sp * 9.81 / stage1.T;
stage1.gamma_rate = deg2rad(50)/stage1.t_burn; % Imposed gamma rate [rad/s]

% Stage 2 parameters
stage2.D = 1.3;            % Diameter [m]
stage2.A = pi * (stage2.D / 2)^2; % Cross-sectional area [m^2]
stage2.Cd = 0.5;           % Drag coefficient
stage2.Cl = 0;             % Lift coefficient
stage2.m0 = 2543 + 400;    % Initial mass [kg] including payload
stage2.I_sp = 336.58;      % Specific impulse [s]
stage2.m_prop = 2187;      % Propellant mass [kg]
stage2.TW = 1.8;           % Thrust-to-weight ratio
stage2.delta = 0;          % Thrust vector angle [rad]
stage2.T = stage2.m0 * 9.81 * stage2.TW; % Thrust [N]
stage2.t_burn = stage2.m_prop * stage2.I_sp * 9.81 / stage2.T;

% Carrier velocity calculation
M_c = 0.82;               % Mach number, carrier velocity relative to the speed of sound
h_c = 11.88e3;            % Carrier altitude (11.88 km)
gradient = 0.0065;        % Standard atmospheric temperature lapse rate
T_amb = 15 + 273.15;      % Ambient temperature at sea level in Kelvin (15°C)
T_c = T_amb - h_c * gradient;  % Temperature at the carrier's altitude, calculated using the lapse rate
gamma = 1.4;              % Specific heat ratio (adiabatic index)
R = 287;                  % Specific gas constant for air in J/(kg·K)
a_c = sqrt(gamma * R * T_c);   % Speed of sound at the carrier's altitude, using the equation for sound speed
v_c = M_c * a_c;          % Carrier's velocity, calculated as Mach number times speed of sound

% Trajectory optimization bounds
lb = [200, deg2rad(0.001), 9e3, 0, 0, 0]; % Lower bounds
ub = [2200, deg2rad(2), 12e3, stage1.t_burn, 500, stage2.t_burn]; % Upper bounds

% ========================================================================
%% OPTIMIZATION
% ========================================================================
initial_guess = [220, deg2rad(1), 11e3, stage1.t_burn-10, 350, stage2.t_burn-10]; % [v0, gamma0, h0, t_bo1, delta_t, t_bo2]

options = optimoptions("fmincon", ...
    "Algorithm", "sqp", ...
    'Display', 'iter-detailed', ...
    'MaxFunctionEvaluations', 5e3, ...
    'MaxIterations', 800, ...
    'ConstraintTolerance', 1e-10, ...
    'StepTolerance', 1e-8);

% Run optimization
[x_opt, delta_v_opt] = fmincon(@(x) ObjFun(x, stage1, stage2), ...
    initial_guess, [], [], [], [], lb, ub, ...
    @(x) myconstraints(x, stage1, stage2), options);

% ========================================================================
%% RESULTS DISPLAY
% ========================================================================
fprintf(' \nOptimized delta-v loss: %+.5e [m/s]\n', delta_v_opt);

% Extract optimized parameters
v_opt = x_opt(1);            % Optimized velocity [m/s]
gamma_opt = rad2deg(x_opt(2)); % Flight path angle [deg]
h_opt = x_opt(3) / 1e3;      % Optimized altitude [km]
t_bo1_opt = x_opt(4);        % Burnout time 1 [s]
delta_t_opt = x_opt(5);      % Coast phase duration [s]
t_bo2_opt = x_opt(6);        % Burnout time 2 [s]

fprintf([' \nOptimized initial state:\n' ...
    'Initial velocity: %+.5e [m/s]\n' ...
    'Initial pitch angle: %+.5e [deg]\n' ...
    'Initial altitude: %+.5e [km]\n' ...
    'Burnout time 1: %+.5e [s]\n' ...
    'Coast phase duration: %+.5e [s]\n' ...
    'Burnout time 2: %+.5e [s]\n'], ...
    v_opt, gamma_opt, h_opt, t_bo1_opt, delta_t_opt, t_bo2_opt);

%% PROPAGATION

% New state to propagate: [initial velocity, flight path angle, altitude, initial angular velocity]
xx_0_NEW = [x_opt(1:3), 0];

% Set the options for ODE solvers
options_pu = odeset('reltol', 1e-5, 'abstol', 1e-6, 'Events', @(t,y) Karman_line(t, y, true)); 
options = odeset('reltol', 1e-5, 'abstol', 1e-6);  % Options for the other phases without event function

% Step 1: Integrate the trajectory for the pitch-up phase
[tt1, xx1] = ode78(@(t,x) fun_pitchup(t,x,stage1), [0 t_bo1_opt], xx_0_NEW, options_pu); 

% Display the results of the pitch-up phase
fprintf('\nResults at the end of the pitch-up phase:\n')
fprintf('Final velocity: %+.5e [m/s]\n', xx1(end, 1))
fprintf('Final flight path angle: %+.5e [deg]\n', rad2deg(xx1(end, 2)))
fprintf('Final altitude: %+.5e [m]\n', xx1(end, 3))

% Step 2: Integrate the coasting phase after stage separation
[tt_off, xx_off] = ode78(@(t,x) fun_off(t,x,stage2), [0 delta_t_opt], xx1(end,:), options);

% Display the results of the coasting phase
fprintf('\nResults at the end of the coasting phase:\n')
fprintf('Final velocity: %+.5e [m/s]\n', xx_off(end, 1))
fprintf('Final flight path angle: %+.5e [deg]\n', rad2deg(xx_off(end, 2)))
fprintf('Final altitude: %+.5e [m]\n', xx_off(end, 3))

% Step 3: Integrate the final phase with stage 2 propulsion
[tt2, xx2] = ode78(@(t,x) fun(t,x,stage2), [0 t_bo2_opt], xx_off(end,:), options);

% Display the results of the stage 2 phase
fprintf('\nResults at the end of the stage 2 phase:\n')
fprintf('Final velocity: %+.5e [m/s]\n', xx2(end, 1))
fprintf('Final flight path angle: %+.5e [deg]\n', rad2deg(xx2(end, 2)))
fprintf('Final altitude: %+.5e [m]\n', xx2(end, 3))

% Constants for orbital insertion calculation
mu = 398600.435;  % Gravitational parameter for Earth in km^3/s^2
R_E = 6.37813660e+06;  % Earth's radius in meters
h_f = 400e3;  % Final altitude (in meters) for orbit insertion
r_f = R_E + h_f;  % Final orbital radius (from Earth's center to the target altitude)
gamma_f = 0;  % Final flight path angle

% Calculate the required velocity for orbital insertion using the vis-viva equation
v_insertion = sqrt(mu / (r_f / 1e3)) * 1e3;  % Convert from km/s to m/s

% Calculate the errors 
err_v = abs(xx2(end, 1) - v_insertion);  % The final velocity error (in m/s)
err_gamma = abs(xx2(end, 2) - gamma_f);  % The final velocity error (in m/s)
err_h = abs(xx2(end, 3) - h_f);  % The final velocity error (in m/s)

% Output the error in velocity to the console
fprintf(' \n The error in orbit insertion velocity is : %+.5e [m/s] \n', err_v)
fprintf(' \n The error in final flight path angle is : %+.5e [rad] \n', err_gamma)
fprintf(' \n The error in final altitude is : %+.5e [m] \n', err_h)

% Check if the burn-out time exceeds the time from the ODE solution
if t_bo1_opt > tt1(end)
    % If the burnout time is greater than the time of the ODE integration (i.e., target altitude reached),
    % adjust the burn-out time to match the event condition (altitude).
    t_bo1_opt_pre = t_bo1_opt;
    fprintf(' \n First stage burn out due to altitude target : delta_t = %+.5e [s] \n', (t_bo1_opt - tt1(end)))
    t_bo1_opt = tt1(end);  % Update the burn-out time based on the event function
end

%% Singular Plots: Stage-wise Dynamics
figure('Name', 'Stage Dynamics', 'NumberTitle', 'off')

% Velocity Plot
subplot(2, 2, 1)
plot(tt1, xx1(:, 1), 'LineWidth', 2); hold on;
plot(tt_off + tt1(end), xx_off(:, 1), 'LineWidth', 2);
plot(tt2 + tt_off(end) + tt1(end), xx2(:, 1), 'LineWidth', 2);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Velocity [m/s]', 'Interpreter', 'latex', 'FontSize', 14);
title('Velocity vs. Time', 'Interpreter', 'latex');
legend(["First stage pitch-up", "Separation and coasting", "Second stage"], ...
    'Location', 'best', 'FontSize', 12);
grid on;

% Flight Path Angle
subplot(2, 2, 2)
plot(tt1, rad2deg(xx1(:, 2)), 'LineWidth', 2); hold on;
plot(tt_off + tt1(end), rad2deg(xx_off(:, 2)), 'LineWidth', 2);
plot(tt2 + tt_off(end) + tt1(end), rad2deg(xx2(:, 2)), 'LineWidth', 2);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Flight Path Angle [deg]', 'Interpreter', 'latex', 'FontSize', 14);
title('Flight Path Angle vs. Time', 'Interpreter', 'latex');
grid on;

% Altitude Plot
subplot(2, 2, 3)
plot(tt1, xx1(:, 3), 'LineWidth', 2); hold on;
plot(tt_off + tt1(end), xx_off(:, 3), 'LineWidth', 2);
plot(tt2 + tt_off(end) + tt1(end), xx2(:, 3), 'LineWidth', 2);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Altitude [km]', 'Interpreter', 'latex', 'FontSize', 14);
title('Altitude vs. Time', 'Interpreter', 'latex');
grid on;

% Angular Velocity
subplot(2, 2, 4)
plot(tt1, rad2deg(xx1(:, 4)) , 'LineWidth', 2); hold on;
plot(tt_off + tt1(end), rad2deg(xx_off(:, 4)) , 'LineWidth', 2);
plot(tt2 + tt_off(end) + tt1(end), rad2deg(xx2(:, 4)) , 'LineWidth', 2);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Angular Velocity [deg/s]', 'Interpreter', 'latex', 'FontSize', 14);
title('Angular Velocity vs. Time', 'Interpreter', 'latex');
grid on;

%% Combined Plot: Velocity and Flight Path Angle
figure('Name', 'Combined Dynamics', 'NumberTitle', 'off');
tt_combined = [tt1; tt_off + tt1(end); tt2 + tt_off(end) + tt1(end)];
xx_combined = [xx1; xx_off; xx2];

yyaxis left
plot(tt_combined, xx_combined(:, 1), 'LineWidth', 2);
ylabel('Velocity [m/s]', 'Interpreter', 'latex', 'FontSize', 14);

yyaxis right
plot(tt_combined, rad2deg(xx_combined(:, 2)), 'LineWidth', 2);
ylabel('Flight Path Angle [deg]', 'Interpreter', 'latex', 'FontSize', 14);

xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
xline(tt1(end), '--', 'LineWidth', 1.5, 'DisplayName', 'Stage Separation');
title('Combined Dynamics: Velocity and Flight Path Angle', 'Interpreter', 'latex');
legend(["Velocity", "Flight Path Angle"], 'Location', 'best', 'FontSize', 12);
grid on;

%% Delta-V Losses

% Extracting the state vectors for each phase (pitch-up, coasting, stage 2)
v_vect1 = xx1(:,1);        % Velocity during the pitch-up phase
gamma_vect1 = xx1(:,2);    % Flight path angle during pitch-up
h_vect1 = xx1(:,3);        % Altitude during pitch-up

v_vect_off = xx_off(:,1);  % Velocity during the coasting phase
gamma_vect_off = xx_off(:,2);  % Flight path angle during coasting
h_vect_off = xx_off(:,3);  % Altitude during coasting

v_vect2 = xx2(:,1);        % Velocity during stage 2 propulsion
gamma_vect2 = xx2(:,2);    % Flight path angle during stage 2
h_vect2 = xx2(:,3);        % Altitude during stage 2

% Coefficient for drag calculation (Manoeuvre intensity with sign)
Cb1 = stage1.m0 / (stage1.A * stage1.Cd);  % For stage 1
Cb_off = stage2.m0 / (stage2.A * stage2.Cd);  % For coasting phase
Cb2 = stage2.m0 / (stage2.A * stage2.Cd);  % For stage 2

% Constants for atmospheric density and gravitational acceleration
H_scale = 7500;  % Scale height for atmospheric density model [m]
R_E = 6.37813660e+06;  % Earth's radius [m]

% Function to calculate air density as a function of altitude
rho = @(h) 1.225 * exp(-h ./ H_scale);

% Function to calculate the mass as a function of time for both stages
m1 = @(t) stage1.m0 - stage1.T ./ (stage1.I_sp * 9.81) .* t;  % Stage 1 mass
m2 = @(t) stage2.m0 - stage2.T ./ (stage2.I_sp * 9.81) .* t;  % Stage 2 mass

% Function to calculate gravitational acceleration as a function of altitude
g = @(h) 9.81 ./ ((1 + h ./ R_E) .^ 2);

% Drag losses for each phase (pitch-up, coasting, stage 2)
term1_d = -0.5 ./ Cb1 .* rho(h_vect1) .* v_vect1 .^ 2 .* stage1.m0 ./ m1(tt1);  % Drag during pitch-up
term_off_d = -0.5 ./ Cb_off .* rho(h_vect_off) .* v_vect_off .^ 2;  % Drag during coasting
term2_d = -0.5 ./ Cb2 .* rho(h_vect2) .* v_vect2 .^ 2 .* stage2.m0 ./ m2(tt2);  % Drag during stage 2

% Gravitational losses for each phase
term1_g = -g(h_vect1) .* sin(gamma_vect1);  % Gravitational loss during pitch-up
term_off_g = -g(h_vect_off) .* sin(gamma_vect_off);  % Gravitational loss during coasting
term2_g = -g(h_vect2) .* sin(gamma_vect2);  % Gravitational loss during stage 2

% Calculate the total drag losses by integrating over the time for each phase
delta_v_drag = abs(trapz(tt1, term1_d)) + abs(trapz(tt_off, term_off_d)) + abs(trapz(tt2, term2_d));

% Calculate the total gravitational losses by integrating over the time for each phase
delta_v_grav = abs(trapz(tt1, term1_g)) + abs(trapz(tt_off, term_off_g)) + abs(trapz(tt2, term2_g));

% Total delta_v losses (sum of drag and gravitational losses)
delta_v_tot = delta_v_drag + delta_v_grav;

% Output the total delta_v losses to the console
fprintf(' \nThe total delta_v losses are : %+.5e [m/s] \n', delta_v_tot)
fprintf(' The drag delta_v losses are : %+.5e [m/s] \n', delta_v_drag)
fprintf(' The gravity delta_v losses are : %+.5e [m/s] \n', delta_v_grav)

%% Remaining Propellant Mass

m1_left = stage1.T / (stage1.I_sp * 9.81) * (stage1.t_burn - t_bo1_opt);
m2_left = stage2.T / (stage2.I_sp * 9.81) * (stage2.t_burn - t_bo2_opt);

perc1 = m1_left / stage1.m_prop * 100;
perc2 = m2_left / stage2.m_prop * 100;

fprintf('\nRemaining Propellant Mass:\n');
fprintf('  - Stage 1: %.2f [kg] (%.2f%%)\n', m1_left, perc1);
fprintf('  - Stage 2: %.2f [kg] (%.2f%%)\n\n', m2_left, perc2);

%% Acceleration Plot

D1 = @(rho,v) 0.5.*rho.*(v.^2).*stage1.A.*stage1.Cd;
D_off = @(rho,v) 0.5.*rho.*(v.^2).*stage2.A.*stage2.Cd;
D2 = @(rho,v) 0.5.*rho.*(v.^2).*stage2.A.*stage2.Cd;

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

a1 = (stage1.T./m1(tt1) - D1(rho(h_vect1),xx1(:,1))./m1(tt1) - g(h_vect1).*sin(xx1(:,2)))./9.81;
a_off = (D_off(rho(h_vect_off),xx_off(:,1))./stage2.m0 - g(h_vect_off).*sin(xx_off(:,2)))./9.81;
a2 = (stage2.T./m2(tt2) - D2(rho(h_vect2),xx2(:,1))./m2(tt2) - g(h_vect2).*sin(xx2(:,2)))./9.81;

figure('Name', 'Acceleration', 'NumberTitle', 'off');
a_combined = [a1; a_off; a2];
plot(tt_combined, a_combined, 'LineWidth', 2);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Acceleration [$m/s^2$]', 'Interpreter', 'latex', 'FontSize', 14);
title('Acceleration Evolution', 'Interpreter', 'latex');
grid on;

max_a = 0;

if max(a1) > max_a
    max_a = max(a1);
    if max(a_off) > max_a
        max_a = max(a_off);
        if max(a2) > max_a
            max_a = max(a2);
        end
    end
end

fprintf(' \nMax acceleration is : %+.5e [g] \n', max_a)

%% Max-Q Analysis

% PegasusXL --> maxq = 6,923e4 [Pa]

max_q1 = 0.5.*rho(h_vect1).*v_vect1.^2;
max_q_off = 0.5.*rho(h_vect_off).*v_vect_off.^2;
max_q2 = 0.5.*rho(h_vect2).*v_vect2.^2;

figure('Name', 'Max-Q Analysis', 'NumberTitle', 'off');
q_max = [max_q1; max_q_off; max_q2];
plot(tt_combined, q_max, 'LineWidth', 2);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Dynamic Pressure [Pa]', 'Interpreter', 'latex', 'FontSize', 14);
title('Dynamic Pressure (Max-Q)', 'Interpreter', 'latex');
grid on;

max_q = 0;

if max(max_q1) > max_q
    max_q = max(max_q1);
    if max(max_q_off) > max_q
        max_q = max(max_q_off);
        if max(max_q2) > max_q
            max_q = max(max_q2);
        end
    end
end

fprintf(' \nMax q is : %+.5e [Pa] \n', max_q)

%% ---------------------------FUNCTIONS--------------------------------

% Function to set plot styles for better visualization
function plotStyle
    % General interpreter settings for LaTeX rendering
    set(0, 'defaultTextInterpreter', 'Latex');
    set(0, 'defaultAxesTickLabelInterpreter', 'Latex');
    set(0, 'defaultLegendInterpreter', 'Latex');
    
    % Grid settings
    set(0, 'defaultAxesXGrid', 'on');
    set(0, 'defaultAxesYGrid', 'on');
    
    % Line appearance
    set(0, 'defaultLineLineWidth', 1.5);
    set(0, 'defaultLineMarkerSize', 6);
    set(0, 'defaultLineMarkerEdgeColor', 'k');
    set(0, 'defaultLineMarkerFaceColor', 'auto');
    
    % Legend appearance
    set(0, 'defaultLegendLocation', 'northoutside');
    set(0, 'defaultLegendOrientation', 'horizontal');
    set(0, 'defaultLegendFontSize', 12);
    
    % Axes settings
    set(0, 'defaultAxesFontSize', 16);
end

% Dynamics function for the main propulsion phase
function [dydt] = fun(t, y, param)
    % Unpack state variables
    v = y(1); % Velocity [m/s]
    gamma = y(2); % Flight path angle [rad]
    h = y(3); % Altitude [m]
    
    % Unpack parameters
    T = param.T; % Thrust [N]
    A = param.A; % Cross-sectional area [m^2]
    Cd = param.Cd; % Drag coefficient
    Cl = param.Cl; % Lift coefficient
    m0 = param.m0; % Initial mass [kg]
    I_sp = param.I_sp; % Specific impulse [s]
    m_prop0 = param.m_prop; % Initial propellant mass [kg]
    delta = param.delta; % Thrust angle (unused since T=0) [rad]
    
    % Constants
    H_scale = 7500; % Scale height [m]
    R_E = 6.37813660e+06; % Earth radius [m]
    g = 9.81 / (1 + h / R_E)^2; % Gravity at altitude
    rho = 1.225 * exp(-h / H_scale); % Air density
    D = 0.5 * rho * (v^2) * A * Cd; % Drag force
    L = 0.5 * rho * (v^2) * A * Cl; % Lift force
    t_burn = m_prop0 * I_sp * 9.81 / T; % Burn time

    % Compute mass during and after burn
    if t < t_burn
        m = m0 - T / (I_sp * 9.81) * t;
    else
        T = 0; % Thrust off
        m = m0 - m_prop0; % Remaining mass
    end

    % Compute derivatives
    dydt = zeros(4, 1);
    dydt(1) = T/m * cos(delta) - D/m - g * sin(gamma); % dv/dt
    dydt(2) = v * cos(gamma) / (R_E + h) + T * sin(delta) / (m * v) + L/m - g * cos(gamma) / v; % dgamma/dt
    dydt(3) = v * sin(gamma); % dh/dt
    dydt(4) = v * cos(gamma) / (R_E + h); % domega/dt (angular velocity)

    % % Debugging output for key variables
    % if mod(t, 1) < 1e-2 % Print every second approximately
    %     fprintf('Time: %.2f s, Velocity: %.2f m/s, Altitude: %.2f m\n', t, v, h);
    % end
end

% Dynamics function for pitch-up phase
function [dydt] = fun_pitchup(t, y, param)
    % Function to model dynamics during the pitch-up phase
    % INPUTS:
    %   t      - Time [s]
    %   y      - State vector [v, gamma, h] where
    %              v: velocity [m/s]
    %              gamma: flight path angle [rad]
    %              h: altitude [m]
    %   param  - Structure containing physical and aerodynamic parameters
    % OUTPUT:
    %   dydt   - Derivatives of the state vector
    
    % Extract state variables
    v = y(1); % Velocity [m/s]
    gamma = y(2); % Flight path angle [rad]
    h = y(3); % Altitude [m]

    % Unpack parameters
    T = param.T; % Thrust [N]
    A = param.A; % Cross-sectional area [m^2]
    Cd = param.Cd; % Drag coefficient
    m0 = param.m0; % Initial mass [kg]
    I_sp = param.I_sp; % Specific impulse [s]
    gamma_rate = param.gamma_rate; % Rate of change of flight path angle [rad/s]
    m_prop0 = param.m_prop; % Initial propellant mass [kg]
    delta = param.delta; % Thrust angle (unused since T=0) [rad]

    % Constants
    H_scale = 7500; % Atmospheric scale height [m]
    R_E = 6.37813660e+06; % Earth radius [m]
    g = 9.81 / (1 + h / R_E)^2; % Gravity at current altitude [m/s^2]
    rho = 1.225 * exp(-h / H_scale); % Air density [kg/m^3]
    D = 0.5 * rho * (v^2) * A * Cd; % Drag force [N]
    t_burn = m_prop0 * I_sp * 9.81 / T; % Burn duration [s]

    % Initialize dydt
    dydt = NaN(length(y), 1);

    % Determine mass during and after burn
    if t < t_burn
        m = m0 - T / (I_sp * 9.81) * t; % Mass during burn [kg]
    else
        T = 0; % Thrust is off
        m = m0 - m_prop0; % Mass after burn [kg]
    end

    % Compute derivatives
    dydt(1) = T/m * cos(delta) - D/m - g * sin(gamma); % dv/dt
    dydt(2) = gamma_rate; % dgamma/dt (constant pitch rate)
    dydt(3) = v * sin(gamma); % dh/dt
    dydt(4) = v * cos(gamma) / (R_E + h); % domega/dt (angular velocity)

    % % Debugging output for critical variables
    % if mod(t, 1) < 1e-2 % Output every second approximately
    %     fprintf('Pitch-up Phase | Time: %.2f s | Velocity: %.2f m/s | Altitude: %.2f m\n', t, v, h);
    %     fprintf('Drag Force: %.2f N | Mass: %.2f kg | Thrust: %.2f N\n', D, m, T);
    % end
end

% Dynamics function for coasting phase
function [dydt] = fun_off(~, y, param)
    % Function to model dynamics during the coasting phase (no thrust)
    % INPUTS:
    %   y      - State vector [v, gamma, h] where
    %              v: velocity [m/s]
    %              gamma: flight path angle [rad]
    %              h: altitude [m]
    %   param  - Structure containing physical and aerodynamic parameters
    % OUTPUT:
    %   dydt   - Derivatives of the state vector

    % Extract state variables
    v = y(1); % Velocity [m/s]
    gamma = y(2); % Flight path angle [rad]
    h = y(3); % Altitude [m]

    % Unpack parameters
    A = param.A; % Cross-sectional area [m^2]
    Cd = param.Cd; % Drag coefficient
    Cl = param.Cl; % Lift coefficient
    m0 = param.m0; % Total mass [kg]
    delta = param.delta; % Thrust angle (unused since T=0) [rad]

    % Constants
    H_scale = 7500; % Atmospheric scale height [m]
    R_E = 6.37813660e+06; % Earth radius [m]
    g = 9.81 / (1 + h / R_E)^2; % Gravity at current altitude [m/s^2]
    rho = 1.225 * exp(-h / H_scale); % Air density [kg/m^3]
    D = 0.5 * rho * (v^2) * A * Cd; % Drag force [N]
    L = 0.5 * rho * (v^2) * A * Cl; % Lift force [N]

    % Initialize dydt
    dydt = NaN(length(y), 1);

    % No thrust during coasting phase
    T = 0;

    % Compute derivatives
    dydt(1) = T/m0 * cos(delta) - D/m0 - g * sin(gamma); % dv/dt
    dydt(2) = v * cos(gamma) / (R_E + h) + T * sin(delta) / (m0 * v) + L/m0 - g * cos(gamma) / v; % dgamma/dt
    dydt(3) = v * sin(gamma); % dh/dt
    dydt(4) = v * cos(gamma) / (R_E + h); % domega/dt (angular velocity)
    
    % % Debugging output for critical variables
    % if mod(h, 1000) < 1e-2 % Output for every 1000 m altitude approximately
    %     fprintf('Coasting Phase | Altitude: %.2f m | Velocity: %.2f m/s | Flight Path Angle: %.2f rad\n', h, v, gamma);
    %     fprintf('Drag Force: %.2f N | Lift Force: %.2f N | Gravity: %.2f m/s^2\n', D, L, g);
    % end
end

% Objective function for optimization
function [F, G] = ObjFun(x_vect, param_1st, param_2st)
    % Objective function for constrained minimization (fmincon)
    % This function calculates the total deltaV losses due to drag and gravity,
    % and optionally computes the gradient of the function.
    % INPUTS:
    %   x_vect   - Vector of variables to optimize:
    %               x_vect(1) = initial velocity (v0) [m/s]
    %               x_vect(2) = initial flight path angle (gamma0) [rad]
    %               x_vect(3) = initial altitude (h0) [m]
    %               x_vect(4) = burn time for phase 1 (t_bo1) [s]
    %               x_vect(5) = coasting time (delta_t) [s]
    %               x_vect(6) = burn time for phase 2 (t_bo2) [s]
    %   param_1st - Structure containing the first set of parameters for phase 1
    %   param_2st - Structure containing the second set of parameters for phase 2
    %
    % OUTPUTS:
    %   F - Scalar value of the total deltaV losses (drag + gravity)
    %   G - Gradient of the objective function (optional)
    
    % Extract variables from the optimization vector
    v0 = x_vect(1);     % Initial velocity [m/s]
    gamma0 = x_vect(2); % Initial flight path angle [rad]
    h0 = x_vect(3);     % Initial altitude [m]
    t_bo1 = x_vect(4);  % Burn time phase 1 [s]
    delta_t = x_vect(5); % Coasting time [s]
    t_bo2 = x_vect(6);  % Burn time phase 2 [s]
    
    % Initial state vector for propagation
    omega0 = 0;         % Initial angular velocity [rad]
    y0 = [v0; gamma0; h0; omega0];
    
    % Set ODE options
    options_pu = odeset('reltol', 1e-5, 'abstol', 1e-6, 'Events', @(t, y) Karman_line(t, y, true));
    options = odeset('reltol', 1e-5, 'abstol', 1e-6);
    
    % Propagate the trajectory for each phase
    [tt1, xx1] = ode78(@(t, x) fun_pitchup(t, x, param_1st), [0 t_bo1], y0, options_pu);
    [tt_off, xx_off] = ode78(@(t, x) fun_off(t, x, param_2st), [0 delta_t], xx1(end, :), options);
    [tt2, xx2] = ode78(@(t, x) fun(t, x, param_2st), [0 t_bo2], xx_off(end, :), options);
    
    % Extract states from the trajectory
    v_vect1 = xx1(:, 1);
    gamma_vect1 = xx1(:, 2);
    h_vect1 = xx1(:, 3);
    
    v_vect_off = xx_off(:, 1);
    gamma_vect_off = xx_off(:, 2);
    h_vect_off = xx_off(:, 3);

    v_vect2 = xx2(:, 1);
    gamma_vect2 = xx2(:, 2);
    h_vect2 = xx2(:, 3);
    
    % Constants for drag and mass calculations
    Cb1 = param_1st.m0 / (param_1st.A * param_1st.Cd);
    Cb_off = param_2st.m0 / (param_2st.A * param_2st.Cd);
    Cb2 = param_2st.m0 / (param_2st.A * param_2st.Cd);

    % Constants for atmospheric calculations
    H_scale = 7500; % [m]
    R_E = 6.37813660e+06; % [m]
    
    % Define functions for atmospheric density, gravity, and mass decay
    rho = @(h) 1.225 * exp(-h ./ H_scale); % Air density [kg/m^3]
    m1 = @(t) param_1st.m0 - param_1st.T ./ (param_1st.I_sp * 9.81) .* t; % Mass during phase 1
    m2 = @(t) param_2st.m0 - param_2st.T ./ (param_2st.I_sp * 9.81) .* t; % Mass during phase 2
    g = @(h) 9.81 ./ ((1 + h ./ R_E) .^ 2); % Gravity at altitude [m/s^2]

    % Drag and gravity terms for the objective function
    term1_d = -0.5 ./ Cb1 .* rho(h_vect1) .* v_vect1 .^ 2 .* param_1st.m0 ./ m1(tt1);
    term_off_d = -0.5 ./ Cb_off .* rho(h_vect_off) .* v_vect_off .^ 2;
    term2_d = -0.5 ./ Cb2 .* rho(h_vect2) .* v_vect2 .^ 2 .* param_2st.m0 ./ m2(tt2);

    term1_g = -g(h_vect1) .* sin(gamma_vect1);
    term_off_g = -g(h_vect_off) .* sin(gamma_vect_off);
    term2_g = -g(h_vect2) .* sin(gamma_vect2);
    
    % Total deltaV losses due to drag and gravity
    delta_v_drag = abs(trapz(tt1, term1_d)) + abs(trapz(tt_off, term_off_d)) + abs(trapz(tt2, term2_d));
    delta_v_grav = abs(trapz(tt1, term1_g)) + abs(trapz(tt_off, term_off_g)) + abs(trapz(tt2, term2_g));
    
    % Total objective function value (deltaV losses)
    F = delta_v_drag + delta_v_grav;

    % Gradient of the objective function (optional)
    if nargout > 1
        G = []; % Currently no gradient is computed. Add analytical or numerical gradient here if required.
    end
end

% Non-linear constraints for optimization
function [c, ceq, G, Geq] = myconstraints(x_vect, param_1st, param_2st)
    % Nonlinear constraint function for trajectory optimization
    % INPUTS:
    %   x_vect   - Optimization variables [v0, gamma0, h0, t_bo1, delta_t, t_bo2]
    %   param_1st - Parameters for the first phase (e.g., propulsion, drag, etc.)
    %   param_2st - Parameters for the second phase (e.g., coasting, drag, etc.)
    %
    % OUTPUTS:
    %   c     - Inequality constraints (empty if none)
    %   ceq   - Equality constraints [final velocity, final flight path angle, final altitude]
    %   G     - Gradient of inequality constraints (empty if none)
    %   Geq   - Gradient of equality constraints (empty if none)

    % Extract variables from the optimization vector
    v0 = x_vect(1);     % Initial velocity [m/s]
    gamma0 = x_vect(2); % Initial flight path angle [rad]
    h0 = x_vect(3);     % Initial altitude [m]
    t_bo1 = x_vect(4);  % Burn time for phase 1 [s]
    delta_t = x_vect(5); % Coasting time [s]
    t_bo2 = x_vect(6);  % Burn time for phase 2 [s]

    omega0 = 0;         % Initial angular velocity [rad]
    y0 = [v0; gamma0; h0; omega0];  % Initial state vector

    % Constants
    mu = 398600.435;    % Gravitational parameter of Earth [km^3/s^2]
    R_E = 6.37813660e+06; % Earth radius [m]
    h_f_target = 400e3;  % Target final altitude [m]
    r_f_target = R_E + h_f_target; % Target final radius [m]

    % Propagate the trajectory using the ODE solver
    options_pu = odeset('reltol', 1e-5, 'abstol', 1e-6, 'Events', @(t, y) Karman_line(t, y, true));
    options = odeset('reltol', 1e-5, 'abstol', 1e-6);
    [~, xx1] = ode78(@(t, x) fun_pitchup(t, x, param_1st), [0 t_bo1], y0, options_pu);
    [~, xx_off] = ode78(@(t, x) fun_off(t, x, param_2st), [0 delta_t], xx1(end, :), options);
    [~, xx2] = ode78(@(t, x) fun(t, x, param_2st), [0 t_bo2], xx_off(end, :), options);

    % Final state after phase 2
    xx_f = xx2(end, :);
    v_f = xx_f(1);      % Final velocity [m/s]
    gamma_f = xx_f(2);  % Final flight path angle [rad]
    h_f = xx_f(3);      % Final altitude [m]

    % Inequality constraints (currently none)
    c = [];

    % Equality constraints:
    % 1. Final velocity should match the circular orbital velocity at target radius
    ceq(1) = sqrt(mu / (r_f_target / 1e3)) * 1e3 - v_f;  % Circular orbital velocity [m/s]
    
    % 2. Final flight path angle should be zero (horizontal orbit)
    ceq(2) = gamma_f;  % Final flight path angle [rad]
    
    % 3. Final altitude should be 400 km
    ceq(3) = h_f - 400e3;  % Target altitude [m]

    % Gradient of constraints (optional, but currently not computed)
    if nargout > 2
        G = [];   % Gradient of inequality constraints (none)
        Geq = []; % Gradient of equality constraints (none)
    end
end

% Event function to stop integration
function [value, isterminal, direction] = Karman_line(~, y, isTerminal)
    % Event function for stopping the ODE solver when the altitude reaches the Karman line (130 km)
    % INPUTS:
    %   ~         - Time (not used in this function)
    %   y         - State vector [v; gamma; h; omega] (h is the altitude)
    %   isTerminal - Boolean that controls whether to stop the integration
    %
    % OUTPUTS:
    %   value     - The event condition (stop when h = 130 km)
    %   isterminal - Boolean flag indicating if the integration should stop
    %   direction  - The direction of the event trigger (0 means trigger in any direction)
    
    h = y(3); % Extract altitude from the state vector
    value = h - 80e3; % Event condition: trigger when altitude exceeds 130 km
    isterminal = isTerminal; % Stop integration if event is triggered (1 = stop)
    direction = 0; % Event triggers in either direction (h increasing or decreasing)
end


