% This function computes preliminry aerodynamic axial coefficients
% according to the concept of component build-up 

% Inputs: 
% a,b semimajor and semi-minor axes of the cross section of the body in meters 
% x is the coordinate at which cross section (a,b) are provided 
% d_nose nose diameter 
% l_nose length of the nose 
% nose_type: 'char' - ogive, cone...

function [Ca, Ca_w, Ca_sf, Ca_b, y] = CA(a, b, x, l_nose, d_nose, S_nose, nose_type, d_main, A_ex, M, alpha, h, Pw) 
    
    Ca    = 0;
    Ca_w  = 0;
    Ca_sf = 0;
    Cf_i  = 0;
    Ca_b  = 0;
 
    if h > 84000    % Altitude check
        Ca = 0;     % limit for atmospheric effects 
    else
        
        % General parameters
        alpha   = deg2rad(alpha);
        A_ref = (d_main/2)^2*pi;    % reference area
        fn_nose = l_nose/d_nose;    % fineness ratio of nose 
        
        %[~,a,~,~,ni,~]  = atmosisa(h,extended=true);  % for MATLAB 2024
        [~,a,~,~,ni]  = atmosisa(h,extended=true);     % for MATLAB 2023

        v = M*a;
        Re = v*x(end)/ni;           % Reynolds number
        
        % Deal with AoA 
        if alpha <= deg2rad(90) && alpha >= 0 
                alpha = + alpha;
        elseif alpha <= deg2rad(180)  && alpha >= deg2rad(90)
                alpha = deg2rad(180-alpha);
        end
        
        if alpha <= deg2rad(90)
            theta = atan(1/2/fn_nose);  % cone half-angle
        elseif alpha > deg2rad(90)
            theta = pi/2;
        end
        
        % Nose type definition 
        switch nose_type

            case 'C'    % conical 
                if M <= 1  
                    Ca_w = 0.8*sin(theta)^2;    
                elseif M > 1
                    beta = sqrt(M^2 - 1);
                    Ca_w = (4*sin(theta)^2*(2.5+8*beta*sin(theta)))/(1+16*beta*sin(theta)); % Linnell-Bailey 
                    %epsilon = atan(1/(2*fn_nose));                         
                    %Ca_w = 2.1*sin(epsilon)^2 + 0.5*sin(epsilon)/beta;
                    %alternative valid for M>1.3
                end
        
            case 'TO'   % tangent ogive
                if M <= 1
                    Ca_w = 0.8*sin(theta)^2;             % in transonic interpolate
                elseif M > 1
                    Ca_w = wvdrogive(d_nose,l_nose,M);       % Miles method 
                end
        
            case 'HVK'  % von karman 
                % C3 extrapolation (See VK_fitting.m)
                load("VK_poly_coeff.mat", 'p1', 'p2', 'p3')
                if M<0.9
                    C3 = 0.005;
                elseif M>=0.9 && M<1.1
                    C3 = polyval(p1,M);
                elseif M>=1.1 && M<1.25
                    C3 = 0.075;
                elseif M>=1.25 && M<1.45
                    C3 = polyval(p2,M);
                elseif M>=1.45 && M<1.95
                    C3 = polyval(p3,M);
                else
                    C3 = 0.08;
                end

                if M<1
                     q_ratio = 1 + M^2/4+M^4/40;
                else
                    q_ratio = 1.84 -0.76/M^2 +0.166/M^4 + 0.065/M^6;
                end
                C0 = 0.85*q_ratio;
                log4_fun = @(x) log(x) / log(4);
                ex = log4_fun(fn_nose+1);
                Ca_w = C0*(C3/C0)^(ex);
        
            case 'CYL'  
                if M<1          % source: OpenRocket technical documentation 
                    q_ratio = 1+ M^2/4 + M^4/40;
                else 
                    q_ratio = 1.84 - 0.76/M^2 + 0.166/M^4 + 0.035/M^6;
                end 
                Ca_w = 0.85*q_ratio;
        
         end
        
        %% ----- 2. BASE PRESSURE ------ 
        
        gamma = 1.4;
        if M >= 1
            Cp_b = 2/(gamma*M^2)*((2/(gamma+1))^(1.4)*(1/M)^2.8*((2*gamma*M^2-gamma+1)/(gamma+1))-1);          % slide 57
            Ca_b = - Cp_b;
        elseif M < 1 
            Ca_b = .12 + .13*M^2;
        end
        
        if Pw == true
            Ca_b = Ca_b*(1-A_ex/A_ref);
        end
        
        %% ----- 3. SKIN-FRICTION ------
        R_s = 5; % Assume a uniform roughness of 5
        Rec = 51 * (R_s/x(end))^(-1.039);
        % Reynolds check
        if Re<10^4
            Cf_i = 1.48*10^(-2);
        elseif Re>=10^4 && Re<Rec
            Cf_i = 1/(1.5*log(Re)-5.6)^2;
        else
            Cf_i = 0.0032*(R_s/x(end))^0.2;
        end
        
        % Compressibility correction
        if M<1
            C_sf = Cf_i*(1-0.1*M^2);
        else
            C_sf = max([Cf_i/(1+0.18*M^2) ,Cf_i/(1+0.15*M^2)^0.58]) ;
        end
        
        A_wet_body = pi*(x(end)-l_nose)*d_nose + S_nose;
        
        Ca_sf = C_sf*(1+1/2/(x(end)/d_main)*A_wet_body)/A_ref;
        
        %% ----- 4. RESULTANT -----
        
        Ca = Ca_w + Ca_sf + Ca_b;
        Ca = Ca*1.1;        % take into account parasitic drag and other drag sources 
        
        Ca = Ca*cos(alpha)^2;
    end
end