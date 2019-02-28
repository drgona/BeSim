function controller = BeCtrl(model, CtrlParam)

if nargin == 0
   buildingType = 'Infrax';  
   ModelOrders.range = [100, 200, 600]; % add any reduced order model you wish to have
   ModelOrders.choice = 200;            % insert model order or 'full' for full order SSM 
   ModelOrders.off_free = 0;            %  augmented model
   reload = 0;
%    construct the model
   model = BeModel(buildingType, ModelOrders, reload); 
end
if nargin < 2
   CtrlParam.use = 1;
%    CtrlParam.precomputed = 1;
   CtrlParam.MPC.use = 0;
   CtrlParam.MPC.Condensing = 1;
   CtrlParam.RBC.use = 0;
   CtrlParam.PID.use = 0;
   CtrlParam.MLagent.use = 0;
end

% controller parameters
controller.use = CtrlParam.use;
% controller.precomputed.use = CtrlParam.precomputed;
controller.MPC.use =    CtrlParam.MPC.use;
controller.MPC.Condensing =    CtrlParam.MPC.Condensing;
controller.RBC.use =    CtrlParam.RBC.use;
controller.PID.use =    CtrlParam.PID.use;
controller.MLagent.use =    CtrlParam.MLagent.use;

% %  % input constraints  [W]
% % controller.umax = ; % 
% % controller.umin = ; % 
% if strcmp(model.buildingType,'Reno')
%     controller.umax = [1680, 685, 154, 1000, 320, 232]'; 
% elseif strcmp(model.buildingType,'Old')
%     controller.umax = [2940, 960, 300, 1400, 460, 253]';
% elseif strcmp(model.buildingType,'RenoLight')
%     controller.umax = [1680, 685, 154, 1000, 320, 232]'/2; 
% % else
% %     disp('no input constraints');
% end

fprintf('\n------------------ Controller -----------------------\n');

if not(controller.use)    % precomputed inputs and outputs or real measurements
    fprintf('*** Load pre-computed controls ... \n')
    path = ['../buildings/', model.buildingType];  
    load([path '/preComputed_matlab/preComputedControls.mat']);
    controller.precomputed.U = U;  
    controller.precomputed.Y = Y;    
    fprintf('*** Done.\n') 
    
        %    CTRL DESIGN 
% RBC, MPC, PID, ML, etc
    
    
elseif CtrlParam.MPC.use  
    fprintf('*** Create MPC controller ... \n')

   
if  strcmp(model.buildingType,'HollandschHuys')    
     % horizons
    controller.MPC.N = 32;
    controller.MPC.Nc = 32;
    controller.MPC.Nrp = 32;
    controller.MPC.Ndp = 32;
    % weight diagonal matrices 
    controller.MPC.Qsb = 1e10*eye(model.pred.ny);
    controller.MPC.Qsa = 1e10*eye(model.pred.ny);
    controller.MPC.Qu = 1e0*eye(model.pred.nu);
else 
     % horizons
    controller.MPC.N = 22;
    controller.MPC.Nc = 22;
    controller.MPC.Nrp = 22;
    controller.MPC.Ndp = 22;
    % weight diagonal matrices 
    controller.MPC.Qsb = 1e8*eye(model.pred.ny);
    controller.MPC.Qsa = 1e8*eye(model.pred.ny);
    controller.MPC.Qu = 1e0*eye(model.pred.nu);
end

   
    
    %  MPC optimizer synthesis   
    controller.MPC.optimizer = BeMPCdesign(model, controller.MPC);
       
    fprintf('*** Done.\n')
    
    
elseif CtrlParam.PID.use  
    fprintf('*** Create PID controller ... \n')
    
    
    
    fprintf('*** Done.\n')

    
elseif CtrlParam.RBC.use      %% RBC heat curve controller
    fprintf('*** Create RBC controller ... \n')
    
    controller.RBC.w = 0.5; %  on off thermostat width of the switching zone zone
    controller.RBC.zone = 2; % zone = choose location of the on-off thermostat (output)   
    
    fprintf('*** Done.\n')

    
% elseif CtrlParam.ML.use  
    
    
end





end