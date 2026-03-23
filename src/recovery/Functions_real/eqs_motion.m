function [dydt,L,D,rho,throttle] = eqs_motion(t,y,t_push)

global S_ref g0 Thr Isp mf mbb1 tbb1                                       %Call global variable

v = y(1); gam = y(2); h = y(3); x =  y(4); m =  y(5);                      %extract variable

% GRAVITATIONAL ACCELERATION
G = 6.67430e-11;                                                           %gravitational constant [Nm^2/kg^2]
Re = 6371e3;                                                               %earth radius [m]
M = 5.972e24;                                                              %earth mass [kg]
g = G * M / (Re+h)^2;                                                      %gravitational acceleration at given altitude [m/s^2]

% AERODYNAMIC FORCES
rho = 0; Ma=0; CL=0; CD=0;
if h<84852 && h>0                                                          %extract density
    [~,a,~,rho] = atmosisa(h,'extended','on');
    Ma = v/a;
    [CL,CD] = AeroCoeff(0,Ma,h);
end



L = 0.5*rho*v^2 * S_ref*CL;                                                %aerodynamic forces
D = 0.5*rho*v^2 * S_ref*CD;

%THRUST
T = 0; m_dot = 0; flag_push = 0; delt = +pi; throttle = 0;                 %FIRST BOOSTBACK
k1 = 0.85+0.00; g_tresh = 5.00-0;
if t>tbb1   
    if m>mbb1                                                              %If true, push
        T = Thr*6.00*0.70;
        m_dot     = -T/(Isp*g0);
    end
end

if t>t_push                                                                                                        
    if m>mf                                                                %SECOND BOOSTBACK
        throttle = k1 * (g_tresh-sqrt(L^2+D^2)/(m*g0))/g_tresh;

        if throttle<0.5
            throttle = 0.5;
        elseif throttle>1
            throttle = 1;
        end

        T = throttle * Thr*3;
        m_dot     = -T/(Isp*g0);
    end
end


%EQUATION OF MOTION IN CANONICAL FORM
dydt    = zeros(4,1);

dydt(1) = T/m*cos(delt)-D/m-g*sin(gam);                                      
dydt(2) = ( (v^2/(Re+h)-g)*cos(gam) +L/m )/v+T/(m*v)*sin(delt);
dydt(3) = v*sin(gam);
dydt(4) = Re/(Re+h) * v*cos(gam);
dydt(5) = m_dot;

D = D + T;