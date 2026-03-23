% Von Karman HAACK ogive
clc, clear all, close all

R_base = 0.65;
finess_ratio = 2;
L_nose = R_base*2*finess_ratio;
C = 0;

x = linspace(0,L_nose);

theta = acos(1-2*x./L_nose);
y = R_base/sqrt(pi)*sqrt(theta-sin(2*theta)/2+C*sin(theta).^3);

% Nose surface area and volume computation (Chad O’Brien, 2014)
K = sqrt(1+diff(y).^2);
S = trapz(x(2:end),2*pi*y(2:end).*K);

% Nose planar surface
Ap_nose = 2*trapz(x,y);

% Nose internal volume (formule per solidi di rotazione)
V = trapz(x,2*pi*y.^2);


color = [0.8500 0.3250 0.0980];
plot(x,y,'Color',color)
grid on
hold on
plot(x,-y,'Color',color)
axis("equal")

z = zeros(1,length(x));

Points = [z' y' -x']*1000;

writematrix(Points,'Points_SW.txt')