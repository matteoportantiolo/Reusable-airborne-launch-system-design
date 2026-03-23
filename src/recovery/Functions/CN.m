% this function computes preliminary aerodynamic normal coefficients
% according to the concept of component build-up and Jorgensen method

% Inputs: 
% a,b semimajor and semi-minor axes of the cross section of the body in meters 
% x is the coordinate at which cross section (a,b) are provided 
% d_nose nose diameter 
% l_nose length of the nose 

function [Cn, Cnf_alpha] = CN(a, b, x, l_nose, d_nose, Ap_nose, d_main, M, alpha, h)  
    
    if h > 84000
        Cn = 0;     % limit for atmospheric effects 
    else
        
        alpha   = deg2rad(alpha);     % [rad]
        l       = x(end);

        fn_nose = l_nose/d_nose;    % fineness of the nose 
        
        % Deal with AoA 
        
        if alpha <= deg2rad(90) && alpha >= 0 
                alpha = + alpha;
        elseif alpha <= deg2rad(180)  && alpha >= deg2rad(90)
                alpha = deg2rad(180-alpha);
        end
        
        % areas definition 
        Ab = d_main^2*pi/4;
        Ap = d_main*(l-l_nose)+ Ap_nose;
        Ar = d_main^2*pi/4;
        
        %% Definizione di Cd estrapolato male da un grafico (rifare meglio se necessario)  
        C_sub = 1.2;       % Minimum drag coefficient, keeping it at 1.2 for subsonic
        C_peak = 2.2;      % Maximum drag coefficient at the transonic peak
        M0 = 1.0;          % Mach number where the peak occurs
        sigma = 0.15;      % Controls the width of the peak in the transonic region
        C_sup = 1.2;       % Supersonic asymptotic drag coefficient, set to match the subsonic minimum
        a = 0.4;

        %Cd = 1.2 not true in transonic and for certain values of crossflow Re 
            if M< 0.8
                % Subsonic region
                Cd = C_sub;
            elseif M >= 0.8 && M < 1.0
                % Smooth transition region from subsonic to transonic peak
                % Using a smooth exponential function to gradually increase
                Cd = C_sub + (C_peak - C_sub) * (1 - exp(-((M - 0.8) / 0.2)));
            elseif M >= 1.0 && M <= 1.5
                % Transonic peak and smooth transition to avoid downward peak
                % Value at M=1 is directly at the peak
                if M == 1.0
                    Cd = C_peak;
                else
                    % Interpolate smoothly between peak at M=1 and supersonic start at M=1.5
                    C_start = C_peak;
                    C_end = C_sup + a / (1.5^2);
                    Cd = C_start + (C_end - C_start) * (M - 1) / 0.5; % Linear interpolation
                end
            else
                % Supersonic region
                Cd = C_sup + a / (M^2);
            end
        
        if M >= 1
           eta = 1; 
        elseif M < 1 
           eta = 0.05*l/d_main + 0.52;
           beta = sqrt(1-M^2);
        end
        
        C_sb = 1;
        C_newt = 1;
        % They change for winged configuration
        
        %% ----- RESULTANT -----
        
        Cn = Ab/Ar*sin(2*alpha)*cos(alpha/2)*C_sb +eta*Cd*Ap/Ar*sin(alpha)^2*C_newt;

    end
end

