function [c , ceq ] =  v_term_constr(t_push,y0)

global g0

options = odeset('abstol',1e-4,'reltol',1e-4, 'MaxStep',1,'Events',@(t,y) h_zero(t,y) ); %Propagate initial guess
[t,y_out] = ode45(@(t,y) eqs_motion(t,y,t_push), [0 5000] , y0 , options );

c   = y_out(end,1)-550;
ceq = [];

% c   = [];
% ceq = y_out(end,1)-500;