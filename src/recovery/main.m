clear all
close all
clc

set(0,'defaultFigureRenderer','painters');
set(0,'defaultTextInterpreter','latex');  
set(0,'defaultAxesTickLabelInterpreter','latex');  
set(0,'defaultLegendInterpreter','latex');
addpath('.\Functions');
size_font = 22;
size_tit  = 30;

%% DATA DEFINITION
global d S_ref g0 h_stop Thr Isp mf mbb1 tbb1                                                               

y0 = [2.98369e+03 , deg2rad(3.52433e+01) , 9.00180e+04 , 0 , 3416.51];
mf = 1997; 

Thr    = 2.89972*1e+5/9;                                                   %Thrust                      [N]
Isp    = 323.16;                                                           %Isp                         [s]
mf     = 1997;                                                             %Final mass                  [Kg]
tbb1   = 30;                                                               %First boostback             [s]
mbb1   = mf + (y0(end)-mf)*0.75;                                           %Mass after first boostback  [Kg]
d      = 1.6;                                                              %Diameter                    [m]
S_ref  = pi*d^2/4;                                                         %Reference surface           [m^2]
g0     = 9.8067;                                                           %Gravity acceleration        [m/s^2]
h_stop = 15000;                                                            %Stop integration            [m]

%% INTEGRATION

t_push_opt = 352+12;

options = odeset('abstol',1e-4,'reltol',1e-4, 'MaxStep',0.1,'Events',@(t,y) h_zero(t,y) );
[t,y_out] = ode45(@(t,y) eqs_motion(t,y,t_push_opt), [0 5000] , y0 , options );


%% LOAD COMPUTATION
v = y_out(:,1); gam = rad2deg(y_out(:,2)); h = y_out(:,3); x = y_out(:,4); m = y_out(:,5);

L = zeros(length(t),1); D = zeros(length(t),1); rho = zeros(length(t),1);  %Initialize vectors
for i=1:length(t)                                                          %Retrieve L,D, and rho
    [~,L(i),D(i),rho(i),T(i),g(i)] = eqs_motion(t(i), y_out(i,:),t_push_opt);
end

load_fact = sqrt(L.^2+D.^2)./(m*g0);                                      %Compute load factor           [g]
q_stag    = 1.85*(0.5.*rho.*v.^2)*1e-3;                                    %Compute stagnation point load [KPa]
heat_flux = 1.7416*1e-4.*sqrt(rho/0.5).*v.^3.*1e-3;                        %Compute heat flux             [kW/m2]

load_max  = 6;                                                             %Maximum load factor      [g]
q_max     = 70;                                                            %Maximum dynamic pressure [Kpa]
heat_max  = 600;                                                           %Maximum heat flux        [kW/m2]

%% LOSSES COMPUTADION
mpbb1 = y0(end)-mbb1;                                                      %First boostback
dv1   = g0*Isp*log(y0(end)/mbb1);

mpbb2 = mbb1-mf;                                                           %Second boostback
dv2   = g0*Isp*log((mf+mpbb2)/mf);


dv_tot  = v(end)-y0(1);

dv_prop = -(dv1+dv2);
dv_drag = -trapz(t,(D-T')./m);
dv_grav =  trapz(t,-g'.*sin(deg2rad(gam)));

%% PLOT AND PRINT

thrust_flag = T>0;

unpowered_v = v;
powered_v   = v;
unpowered_v(thrust_flag) = NaN;
powered_v(~thrust_flag)  = NaN;

unpowered_gam = gam;
powered_gam   = gam;
unpowered_gam(thrust_flag) = NaN;
powered_gam(~thrust_flag)  = NaN;

unpowered_h = h;
powered_h   = h;
unpowered_h(thrust_flag) = NaN;
powered_h(~thrust_flag)  = NaN;

unpowered_x = x;
powered_x   = x;
unpowered_x(thrust_flag) = NaN;
powered_x(~thrust_flag)  = NaN;

unpowered_load = load_fact;
powered_load   = load_fact;
unpowered_load(thrust_flag) = NaN;
powered_load(~thrust_flag)  = NaN;

unpowered_q_stag = q_stag;
powered_q_stag   = q_stag;
unpowered_q_stag(thrust_flag) = NaN;
powered_q_stag(~thrust_flag)  = NaN;

unpowered_thrust = T;
powered_thrust   = T;
unpowered_thrust(thrust_flag) = NaN;
powered_thrust(~thrust_flag)  = NaN;

unpowered_m = m;
powered_m   = m;
unpowered_m(thrust_flag) = NaN;
powered_m(~thrust_flag)  = NaN;

for i=1:length(t)-1                                                        %Fix jump
    if (thrust_flag(i) == 0 & thrust_flag(i+1)==1)
        powered_thrust(i)   = unpowered_thrust(i);
        unpowered_thrust(i) = NaN;

        powered_load(i)   = unpowered_load(i);
        unpowered_load(i) = NaN;
    end

    if (thrust_flag(i) == 1 & thrust_flag(i+1)==0)
        powered_thrust(i+1)   = unpowered_thrust(i+1);
        unpowered_thrust(i+1) = NaN;

        powered_load(i+1)   = unpowered_load(i+1);
        unpowered_load(i+1) = NaN;
    end
end

figure                                                                     %Plot results
    subplot(2,2,1)
        plot(t,unpowered_v,t,powered_v,'LineWidth',2);
        tit1=title('VELOCITY');
        ylabel('V [m/s]');
        xlabel('t [s]');
        grid on;
        box on;
        xlim([min(t) max(t)]);
    subplot(2,2,2)
        plot(t,unpowered_gam,t,powered_gam,'LineWidth',2);
        tit2=title('GAMMA');
        ylabel('$\gamma [deg]$');
        xlabel('t [s]');
        grid on;
        box on;
        xlim([min(t) max(t)]);
        legend('UNPOWERED','POWERED'); 
    subplot(2,2,3)
        plot(t,unpowered_h,t,powered_h,'LineWidth',2);
        tit3=title('ALTITUDE');
        ylabel('h [km]');
        xlabel('t [s]');
        grid on;
        box on;
        xlim([min(t) max(t)]);
    subplot(2,2,4)
        plot(t,unpowered_x,t,powered_x,'LineWidth',2);
        tit4=title('DOWNRANGE');
        ylabel('r [km]');
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
    subplot(1,2,1)
        hold on
        plot(t,unpowered_load,t,powered_load,'LineWidth',2);
        yline(load_max,'--r','LineWidth',2)
        tit1=title('LOAD FACTOR');
        ylabel('g');
        xlabel('t [s]');
        grid on;
        box on;
        xlim([min(t) max(t)]);
        ylim([0 6.5]);
    subplot(1,2,2)
        hold on
        plot(t,unpowered_q_stag,t,powered_q_stag,'LineWidth',2);
        yline(q_max,'--r','LineWidth',2)
        tit2=title('STAGNATION POINT q');
        ylabel('kPa');
        xlabel('t [s]');
        grid on;
        box on;
        xlim([min(t) max(t)]);
        ylim([0 q_max+5]);
        legend('UNPOWERED','POWERED'); 
    set(findall(gcf,'-property','FontSize'),'FontSize',size_font);
    tit1.FontSize = size_tit;
    tit2.FontSize = size_tit;

figure
    subplot(1,2,1)
        plot(t,unpowered_thrust*1e-3,t,powered_thrust*1e-3,'LineWidth',2);
        tit1=title('THRUST');
        ylabel('$T [kN]$');
        xlabel('t [s]');
        grid on;
        box on;
        xlim([min(t) max(t)]);
    subplot(1,2,2)
        plot(t,unpowered_m,t,powered_m,'LineWidth',2);
        tit2=title('MASS');
        ylabel('Kg');
        xlabel('t [s]');
        grid on;
        box on;
        xlim([min(t) max(t)]);
        legend('UNPOWERED','POWERED'); 
    set(findall(gcf,'-property','FontSize'),'FontSize',size_font);
    tit1.FontSize = size_tit;
    tit2.FontSize = size_tit;

%% PRINT SECTION

clc

fprintf('COMPUTATION STOPPED AT %dm WITH %.0fm/s',h_stop,v(end));

fprintf('\n\n');

fprintf('MAXIMUM LOAD FACTOR %.2fg\n',max(load_fact));
fprintf('MAXIMUM DYNAMIC PRESSURE AT THE STAGNATION POINT %.2fkPa',max(q_stag));

fprintf('\n\n');

fprintf('\t\t\t\tINITIAL TIME [s] \t\t\t PROPELLANT CONSUMED [Kg] \t\t\t CORRESPONDING DV [m/s]\n');
fprintf('1st BOOST --> \t\t %.1f \t\t\t\t\t\t\t %.1f \t\t\t\t\t\t\t\t %.1f\n',tbb1,mpbb1,dv1);
fprintf('1st BOOST --> \t\t%.1f \t\t\t\t\t\t\t%.1f \t\t\t\t\t\t\t\t%.1f\n',t_push_opt,mpbb2,dv2);

fprintf('\n\n');

fprintf('DV_tot [m/s] \t\t DV_prop [m/s] \t\t DV_drag [m/s] \t\tDV_grav [m/s]\n');
fprintf('  %.1f \t\t\t\t%.1f \t\t\t%.1f \t\t   %.1f',dv_tot,dv_prop, dv_drag,dv_grav)

fprintf('\n\n');
