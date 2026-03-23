function [C_L, C_D, center_of_pressure, total_moment] = fin(M, alpha, d_main, chord_root, chord_tip, span, x_CoM, x_le, h)
    % Inputs

    %[~,a,~,rho_inf,~,~]  = atmosisa(h,extended=true);   % for MATLAB 2024
    [~,a,~,rho_inf,~]  = atmosisa(h,extended=true);      % for MATLAB 2023

    alpha = deg2rad(alpha); % Convert angle of attack to radian
    V_inf = M*a;       % Freestream velocity (m/s, assuming Mach=1 at 343 m/s)
    gamma = 1.4;
    C_D_base = 0.1; 
    A_ref = (d_main/2)^2*pi;
    
    % Fin gemetrical properties
    area = 0.5 * (chord_root + chord_tip) * span; % Trapezoidal area
    AR = span^2 / area;                           % Compute Aspect Ratio (AR)
    
    % Oswald efficiency factor for supersonic flow (assume e = 0.9)
    e = 0.9;

    % Interference factor for fins
    Kf_D = 1.04;        % Drag interference factor (Governale, 2021)

    if M <= 1
        % Subsonic case
        C_L = 2 * pi * alpha; % Thin airfoil theory
        C_D = C_D_base + (C_L^2) / (pi * AR * e); % Viscous + induced drag
    else
        % Supersonic case
        beta = sqrt(M^2 - 1);
        C_L = 4 * alpha / beta; % Supersonic thin airfoil theory
        C_D_wave = (gamma / 2) * (alpha^2) / beta; % Wave drag term
        C_D = C_D_base + (C_L^2) / (pi * AR * e) + C_D_wave; % Total drag
    end

    % Interference factor correction
    C_D = C_D*Kf_D;

    % Compute dynamic pressure (q_inf)
    q_inf = 0.5 * rho_inf * V_inf^2;

    % Compute Lift and Drag Forces
    fin_lift = C_L * q_inf * area*2;
    fin_drag = C_D * q_inf * area*2;

    % Rescaling the aerodynamic coefficients onto the launcher reference
    % area
    C_L = C_L * area/A_ref;
    C_D = C_D * area/A_ref;

    % Center of Pressure (corrected for trapezoidal fin)
    % Step 1: Mean Aerodynamic Chord (MAC)
    MAC = (2/3) * (chord_root + chord_tip - (chord_root * chord_tip) / (chord_root + chord_tip));
    
    % Step 2: CoP is at 25% of the MAC from the leading edge
    center_of_pressure = x_le + 0.25 * MAC;

    % Moment Arms
    x_L = center_of_pressure - x_CoM; % Distance from CoM to lift force

    % Compute Pitching Moment
    total_moment = fin_lift * x_L;

end

