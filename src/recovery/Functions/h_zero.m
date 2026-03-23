function [position,isterminal,direction] = h_zero(t,y)
%
% Stop the integration when h=0
%
global h_stop,

    position   = y(3)-0.001-h_stop; 
    isterminal = 1;
    direction  = 0;
end