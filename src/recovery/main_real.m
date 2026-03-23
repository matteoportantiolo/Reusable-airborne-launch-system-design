clear all
close all
clc

set(0,'defaultFigureRenderer','painters');
set(0,'defaultTextInterpreter','latex');  
set(0,'defaultAxesTickLabelInterpreter','latex');  
set(0,'defaultLegendInterpreter','latex');
addpath('.\Functions_real');
size_font = 22;
size_tit  = 30;

%% DATA DEFINITION
global d S_ref g0 h_stop Thr Isp mf mbb1 tbb1

%Release at [V0 , gamma0 , h0 , x0 , m0]                                   %Initial condition
y0 = [2.99954e+03 , deg2rad(3.54840e+01) , 88.95611e+04 , 0 , 3331.22];
Thr = 2.899729849050720e5/9;                                                      %Thrust               [N]
Isp = 323;                                                                 %Isp                  [s]
mf = 1.996190000000000e+03;                                                                 %Final mass           [Kg]



tbb1 = 30;                                                                 %First boostback
mbb1 = mf + (y0(end)-mf)*0.6;                                              %Mass after first boostback
d = 1.6;                                                                   %Diameter             [m]
S_ref = pi*d^2/4;                                                          %Reference surface    [m^2]
g0 = 9.8067;                                                               %Gravity acceleration [m/s^2]
h_stop = 15000;                                                            %Stop integration     [m]

% MINIMIZE g
options = odeset('abstol',1e-4,'reltol',1e-4, 'MaxStep',1,'Events',@(t,y) h_zero(t,y) );
[t,y_out] = ode45(@(t,y) eqs_motion(t,y,[5000 ; 0 ; 0]), [0 5000] , y0 , options );  %Determine max time
lb = [     50 ];
ub = [ t(end) ];
t_push_guess = mean([lb(1) ub(1)])+100;

options = optimoptions('fmincon','Display','iter','Algorithm','active-set','FiniteDifferenceType','central','FiniteDifferenceStepSize',[1e-2],'MaxIter',10);
% [t_push_opt,g_max] = fmincon(@(y) g_min(y,y0),t_push_guess,[],[],[],[],lb,ub,@(y) v_term_constr(y,y0),options);
[t_push_opt,g_max] = fmincon(@(y) g_min(y,y0),t_push_guess,[],[],[],[],lb,ub,'',options);
t_push_opt

%% LOAD OPTIMAL
close all
clear y_out L D rho throt
% y0 = [+3.9e+03 , deg2rad(+3.1e+01) , 109*1e3 , 0 , 3100];
% mbb1 = mf + (y0(end)-mf)*0.6;  
% % 
% t_push_opt = 382+3;

options = odeset('abstol',1e-4,'reltol',1e-4, 'MaxStep',0.1,'Events',@(t,y) h_zero(t,y) );
[t,y_out] = ode45(@(t,y) eqs_motion(t,y,t_push_opt), [0 5000] , y0 , options );

v = y_out(:,1); gam = rad2deg(y_out(:,2)); h = y_out(:,3); x = y_out(:,4); m = y_out(:,5);

L = zeros(length(t),1); D = zeros(length(t),1); rho = zeros(length(t),1);  %Initialize vectors
for i=1:length(t)                                                          %Retrieve L,D, and rho
    [~,L(i),D(i),rho(i),throt(i)] = eqs_motion(t(i), y_out(i,:),t_push_opt);
end

load_fact = sqrt(L.^2+D.^2)./(m*g0);                                      %Compute load factor           [g]
q_stag    = 1.85*(0.5.*rho.*v.^2)*1e-3;                                    %Compute stagnation point load [KPa]
heat_flux = 1.7416*1e-4.*sqrt(rho/0.5).*v.^3.*1e-3;                        %Compute heat flux             [kW/m2]

load_max  = 6;                                                             %Maximum load factor      [g]
q_max     = 70;                                                            %Maximum dynamic pressure [Kpa]
heat_max  = 600;                                                           %Maximum heat flux        [kW/m2]


% PLOT
figure                                                                     %Plot results
    subplot(2,2,1)
        plot(t,v,'LineWidth',2);
        tit1=title('VELOCITY');
        ylabel('V [m/s]');
        xlabel('t [s]');
        grid on;
        box on;
        xlim([min(t) max(t)]);
    subplot(2,2,2)
        plot(t,gam,'LineWidth',2);
        tit2=title('GAMMA');
        ylabel('$\gamma [deg]$');
        xlabel('t [s]');
        grid on;
        box on;
        xlim([min(t) max(t)]);
    subplot(2,2,3)
        plot(t,h,'LineWidth',2);
        tit3=title('ALTITUDE');
        ylabel('h [m]');
        xlabel('t [s]');
        grid on;
        box on;
        xlim([min(t) max(t)]);
    subplot(2,2,4)
        plot(t,x,'LineWidth',2);
        tit4=title('DOWNRANGE');
        ylabel('r [m]');
        xlabel('t [s]');
        grid on;
        box on;
        xlim([min(t) max(t)]);
    set(findall(gcf,'-property','FontSize'),'FontSize',size_font);
    tit1.FontSize = size_tit;
    tit2.FontSize = size_tit;
    tit3.FontSize = size_tit;
    tit4.FontSize = size_tit;

figure                                                                     %Plot loads
    subplot(1,3,1)
        hold on
        plot(t,load_fact,'LineWidth',2);
        yline(load_max,'--r','LineWidth',2)
        tit1=title('LOAD FACTOR');
        ylabel('g');
        xlabel('t [s]');
        grid on;
        box on;
        xlim([min(t) max(t)]);
    subplot(1,3,2)
        hold on
        plot(t,q_stag,'LineWidth',2);
        yline(q_max,'--r','LineWidth',2)
        tit2=title('STAGNATION POINT q');
        ylabel('kPa');
        xlabel('t [s]');
        grid on;
        box on;
        xlim([min(t) max(t)]);
    subplot(1,3,3)
        hold on
        plot(t,heat_flux,'LineWidth',2);
        yline(heat_max,'--r','LineWidth',2)
        tit3=title('HEAT FLUX');
        ylabel('$kW/m^2$');
        xlabel('t [s]');
        grid on;
        box on;
        xlim([min(t) max(t)]);
    set(findall(gcf,'-property','FontSize'),'FontSize',size_font);
    tit1.FontSize = size_tit;
    tit2.FontSize = size_tit;
    tit3.FontSize = size_tit;

figure
    subplot(1,2,1)
        plot(t,throt .* 100,'LineWidth',2);
        tit1=title('THROTTLE');
        ylabel('$\%$');
        xlabel('t [s]');
        grid on;
        box on;
        xlim([min(t) max(t)]);
    subplot(1,2,2)
        plot(t,m,'LineWidth',2);
        tit2=title('MASS');
        ylabel('Kg');
        xlabel('t [s]');
        grid on;
        box on;
        xlim([min(t) max(t)]);
    set(findall(gcf,'-property','FontSize'),'FontSize',size_font);
    tit1.FontSize = size_tit;
    tit2.FontSize = size_tit;

% PRINT SECTION

clc

fprintf('COMPUTATION STOPPED AT %dm WITH %.0fm/s',h_stop,v(end));
% 
% fprintf('\n\n');
% 
% fprintf('\t\t\tV0 [m/s] \t h0 [m]    x0 [m] \t gam0 [°] \t alph0 [°] \t MAX LOAD [g]\n');
% fprintf('INITIAL -->\t  %.0f \t\t %.0f \t %.0f \t      %.1f \t\t   %.1f \t    %.1f\n',y0(1),y0(3),y0(4),rad2deg(y0(2)),y0(5),max(load_first));
% fprintf('FINAL \t-->\t  %.0f \t\t %.0f \t %.0f \t      %.1f \t\t   %.1f \t    %.1f\n',y0_opt(1),y0_opt(3),y0_opt(4),rad2deg(y0_opt(2)),y0_opt(5),g_max(end));
