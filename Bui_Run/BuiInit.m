
%% BuiSim 
% Matlab toolbox for fast developlent, simulation and deployment of
% advanced building climate controllers


% functionality intended for automatic construction of controls and
% estimation for a given linear building model


% TODO long term:
% API for end user
% communication with BMS integrate
% code generation integrate, so implementation is matlab free and matlab is
% used only during design phase


yalmip('clear');

addpath('../Bui_Modeling/')
addpath('../Bui_Disturbances/')
addpath('../Bui_References/')
addpath('../Bui_Estimation/')
addpath('../Bui_Control/')
addpath('../Bui_Simulation/')


%% MODEL   emulator + predictor
% available buildings  'Infrax',  'HollandschHuys', 'Reno', 'Old', 'RenoLight'
buildingType = 'RenoLight';  ModelOrders.range = [4, 7, 10, 15, 20, 30, 40, 100];
% buildingType = 'Infrax'; ModelOrders.range = [100, 200, 600]; 
% buildingType = 'HollandschHuys'; ModelOrders.range = [100, 200, 600]; 
% buildingType = 'Borehole';  ModelOrders.range = [10, 15, 20, 40, 100];  % orderds for borehole
% ModelOrders.choice = 100; 
ModelOrders.choice = 'full';
ModelOrders.off_free = 0;    %  augmented model
reload = 0; 

model = BuiModel(buildingType, ModelOrders, reload);

% TODO:
% HVAC time delays modeling and incporporation into prediction model

% TODO: issue with low order??

%% disturbacnes 

% TODO: HH buidling mismatch in Ed matrix and d0
dist = BuiDist(buildingType, reload);

%% References 
% TODO:  fix to be universal for all models, add multiple options
% TODO: error with the path to common files - fix that
% TODO: works only with Reno, adapt to Infrax
refs = BuiRefs(model);

%%  estimator 
EstimParam.LOPP.use = 0;      %  Luenberger observer via pole placement - Not implemented
EstimParam.SKF.use = 0;    % stationary KF
EstimParam.TVKF.use = 1;   % time varying KF
EstimParam.MHE.use = 0;   % moving horizon estimation via yalmip
EstimParam.MHE.Condensing = 1;   % state condensing 

% TODO: flag for using estimation, if not use perfect state update from plant model
EstimParam.use = 1;

estim = BuiEstim(model, EstimParam);


%% controller 
% CtrlParam.precomputed = 0; % TODO: error with the path to common files - fix that

% TODO: flag for controller use, if not perform open loop simulations 
CtrlParam.use = 1;   % 0 for precomputed u,y    1 for closed loop control
CtrlParam.MPC.use = 1;
CtrlParam.MPC.Condensing = 1;
CtrlParam.RBC.use = 0;
CtrlParam.PID.use = 0;

% TODO: finish implementation of RBC, MPC, PID, ML


ctrl = BuiCtrl(model, CtrlParam);

% TODO higher level tuning params
% PID coefficients
% MPC - just magnutide of weight on comfort and energy 

% if CtrlParam.RBC.use
%     % TODO: mismatch with ny and nu for RBC!!!
%     % TODO: implement this sythematically
%     ctrl.umax = 2000*ones(model.pred.nu,1);
%     ctrl.umax(1:4) = 10;
%     ctrl.umax(25:28) = 0;
%     ctrl.umax([17,16,8,5,11,12]) = 0;
% 
%     ctrl.umin = -2000*ones(model.pred.nu,1);
%     ctrl.umax(1:4) = -5;
%     ctrl.uin(25:28) = 0;
%     ctrl.umin([17,16,8,5,11,12]) = 0;
% end

% TODO:
% implement different levels control
% primary slow control RBC - floors TABS inputs
% secondary fast control RBC - room temperatures direct feedback
% VAVs

% TODO RBC:
% implement heating curve for supply water temp for infrax for RBC
% ctrl and move it to the control block
% identify ambient temp signal from outside and based on that compute T sup

% TODO
% implement original RBC for each model as a part of the model itself
% instead of that implement generic Model structures for development


%% Simulate
SimParam.run.start = 1;
SimParam.run.end = 13; 
% SimParam.run.start = 1;
% SimParam.run.end = 183; 
% sim.run.end = 2; 
SimParam.verbose = 1;
SimParam.flagSave = 0;
SimParam.comfortTol = 1e-1;
% flag distinguishing emulation and real measurements
%  0 - measurements
%  1 - emulation
SimParam.emulate = 1;

%TODO: add open loop (fast results), closed loop (RHC) - done 

% TODO: modify plot option in BuiSim init
PlotParam.flagPlot = 1;     % plot 0 - no 1 - yes
PlotParam.plotStates = 0;        % plot states
PlotParam.plotDist = 1;        % plot disturbances
PlotParam.plotEstim = 1;        % plot estimation
PlotParam.plotCtrl = 1;        % plot control
% PlotParam.Transitions = 1;      % pot dynamic transitions of Ax matrix

% PlotParam.reduced = 0;   %  reduced paper plots formats 0 - no 1 - yes
% PlotParam.zone = 2;     % choose zone if reduced
% PlotParam.only_zone = 0;    %  plot only zone temperatures 0 - no 1 - yes  

profile on 
% %  simulation file with embedded plotting file
outdata = BuiSim(model, estim, ctrl, dist, refs, SimParam, PlotParam);

p = profile('info')
profile viewer


 
 
 