function load_fact = g_min(thr_control,y0)
%
% Compute the maximum load factor for a given initial condition; to be used
% in order to minimize load;
%
global g0

options = odeset('abstol',1e-4,'reltol',1e-4, 'MaxStep',1,'Events',@(t,y) h_zero(t,y) ); %Propagate initial guess
[t,y_out] = ode45(@(t,y) eqs_motion(t,y,thr_control), [0 5000] , y0 , options );

L = zeros(length(t),1); D = zeros(length(t),1); rho = zeros(length(t),1);  %Initialize vectors
for i=1:length(t)                                                          %Retrieve L,D, and rho
    [~,L(i),D(i),rho(i)] = eqs_motion(t(i), y_out(i,:),thr_control);
end

load_fact = sqrt(L.^2+D.^2)./(y_out(:,end)*g0);                            %Compute load factor
load_fact = max(load_fact(round(length(load_fact)*2/3):end));