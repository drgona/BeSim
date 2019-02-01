function [u, heat] = BeRBC(T,R,heat,TSup,ctrl)
% T = actual temperatures in all zones



% THE ONLY DIFFERENCE from on-off control
T_diff = TSup - T; % difference of room temperatures from supply water in current time step
% gain factor: gf = max(TSup)-min(TLow)
% FIXME:  tunung factor should be 60, but there is poor performance
gf = 40;  % manual tuning 

% TOOD: temporary manual fix
% G = ctrl.umax/gf; %radiator gains
G = ctrl.umax(1:length(T_diff))/gf; %radiator gains

PW = G.*T_diff;       %radiator powers


        % logic of the central on-off thermostat
        above = T(ctrl.RBC.zone) >= R + ctrl.RBC.w;
        below = T(ctrl.RBC.zone) <= R - ctrl.RBC.w;
        doheat = (heat & ~above) | (~heat & below);
        %  on off ctrl action
        heat = doheat;
%         u = double(heat)*ones(model.sim.nu,1).*PW;   % ctrl action
%         without P controller
         
% P controller
        dTNom = 0.5;           
        P = max(0, min(1, (R-T)/dTNom ));
        P(ctrl.RBC.zone) = 1; %Remove proportional controller of thermostat zone.
        u = P.*double(heat).*PW; 
        
%         TODO temporal fix
        u(length(T_diff)+1:length(ctrl.umax)) = 0;
        
        u = min(u,ctrl.umax);
          
end
