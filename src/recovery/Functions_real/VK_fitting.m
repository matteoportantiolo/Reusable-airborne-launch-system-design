% Function to shape the pressure drag coefficient of Von Karman Nose shape
% with fr = 3
clc, clear all, close all

M = 0.7:0.01:2;
C3 = zeros(length(M), 1);

% Vectors of point between 0.9 and 1.1
x1 = [0.9 0.95 0.975 1 1.025 1.05 1.075 1.1];
y1 = [0.005 0.01 0.015 0.025 0.05 0.065 0.07 0.075];

p1 = polyfit(x1, y1, 4);

% Costant between 1.1 and 1.25

% Linear interpolation between 1.25 and 1.45
x2 = [1.25 1.45];
y2 = [0.075 0.09];

p2 = polyfit(x2, y2, 1);

% Linear interpolation between 1.45 and 1.95
x3 = [1.45 1.95];
y3 = [0.09 0.08];

p3 = polyfit(x3, y3, 1);

for i=1:length(M)
    if M(i)<0.9
        C3(i) = 0.005;
    elseif M(i)>=0.9 && M(i)<1.1
        C3(i) = polyval(p1,M(i));
    elseif M(i)>=1.1 && M(i)<1.25
        C3(i) = 0.075;
    elseif M(i)>=1.25 && M(i)<1.45
        C3(i) = polyval(p2,M(i));
    elseif M(i)>=1.45 && M(i)<1.95
        C3(i) = polyval(p3,M(i));
    else
        C3(i) = 0.08;
    end
end

plot(M, C3)

save('VK_poly_coeff.mat', 'p1', 'p2', 'p3')