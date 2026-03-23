function [CL,CD] = AeroCoeff(alpha,Mach,h)
%
%function [Cl,Cd] = AeroCoeff(alpha,Mach,h)
%
%INPUT:
%   -alpha          -->         Angle of attack in °    [1x1]
%   -Mach           -->         Mach number             [1x1]
%   -h              -->         Altitude in m           [1x1]
%
%OUTPUT:
%   -Cl             -->         Lift coefficient        [1x1]
%   -Cd             -->         Drag coefficient        [1x1]
%

% Body geometry
a = [0 0.8 0.8]; % Semi-major axis for each station
b = a;           % Semi- minor axis
x_CoM = 13.56;
x_le = 20.75;
d_main = 1.6;    % Main body diameter

% Nose properties
d_nose = 1.6;
x = [0 0.01 15]; % x station
l_nose = x(2);

nose_type = 'CYL';

% Cylinder geometry
l_cyl = x(3)-l_nose;

% Fins geometry (Clipped delta shape - trapezoidal)
AR = 0.6;               % Aspect ratio: typical values from 0.5 to 3; TBD 
chord_root = 1.5;       % Root chord length (m) TBD
chord_tip = 0.5;        % Tip chord length (m) TBD
span = 0.6;             % Fin span (m) TBD
gamma = atand((chord_root-chord_tip)/span);    % Sweep angle: typical values from 0 to 30°; TBD 

x_fin = x(end) + [0  0 -chord_tip -chord_root];
y_fin = d_main/2 + [0 span span 0];

% Geometry shape
x_nose = linspace(0,l_nose);
C = 0;

theta = acos(1-2*x_nose./l_nose);
y_nose = a(2)/sqrt(pi)*sqrt(theta-sin(2*theta)/2+C*sin(theta).^3);

% Nose surface area and volume computation (Chad O’Brien, 2014)
K = sqrt(1+diff(y_nose).^2);
S_nose  = trapz(x_nose(2:end),2*pi*y_nose(2:end).*K);  % Nose wet stuface
Ap_nose = 2*trapz(x_nose,y_nose);                       % Nose planar surface

x_cyl = linspace(l_nose, l_nose + l_cyl, 50); 
y_cyl = ones(size(x_cyl)) * a(2);

% Nozzle exit area
d_ex = 0.2;         % Nozzle diameter
N = 9;              % Number of nozzles
A_ex = N*pi*(d_ex/2)^2;           % Overall exit nozzle area


% Set powered or unpowered phase
Pw = true;


% Compute axial and normal coefficient of the launcher body
[Ca, Ca_w, Ca_sf, Ca_b] = CA(a, b , x, l_nose, d_nose, S_nose, nose_type, d_main, A_ex, Mach, alpha, h, Pw);
[Cn] = CN(a, b, x, l_nose, d_nose, Ap_nose, d_main, Mach, alpha, h);

% Compute lif and drag coefficient of fins scaled on the launcher
% reference area
[CL_fin, CD_fin, XP, moment] = fin(Mach, alpha, d_main, chord_root, chord_tip, span, x_CoM, x_le, h);
% Compute overall lift and drag coefficient
CL = Cn*cos(deg2rad(alpha))-Ca*sin(deg2rad(alpha)) + CL_fin;
CD = Cn*sin(deg2rad(alpha))+Ca*cos(deg2rad(alpha)) + CD_fin;        


