
%% BuiSim 
% Matlab toolbox for fast developlent, simulation and deployment of
% advanced building climate controllers 
% Model Predictive Control (MPC)
% functionality intended for automatic construction of controls and
% estimation for a given linear building model

yalmip('clear');
close all

addpath('../Bui_Modeling/')
addpath('../Bui_Disturbances/')
addpath('../Bui_References/')
addpath('../Bui_Estimation/')
addpath('../Bui_Control/')
addpath('../Bui_Simulation/')
addpath('../Bui_Learn/')
addpath('../Bui_RealTime/')
addpath('../Bui_RealTime/MervisApp')

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

%% Disturbacnes 
% ambient temperature, solar radiation, internal heat gains
dist = BuiDist(buildingType, reload);

%% References 
% comfort constraints, price profiles
RefsParam.Price.variable = 1;       %1 =  variable price profile, 0 = fixed to 1
refs = BuiRefs(model, RefsParam);

%%  estimator 
EstimParam.LOPP.use = 0;      %  Luenberger observer via pole placement - Not implemented
EstimParam.SKF.use = 0;    % stationary KF
EstimParam.TVKF.use = 1;   % time varying KF
EstimParam.MHE.use = 0;   % moving horizon estimation via yalmip
EstimParam.MHE.Condensing = 1;   % state condensing 
EstimParam.use = 1;

estim = BuiEstim(model, EstimParam);

%% controller 
CtrlParam.use = 1;   % 0 for precomputed u,y    1 for closed loop control
CtrlParam.MPC.use = 1;
CtrlParam.MPC.Condensing = 1;
CtrlParam.RBC.use = 0;
CtrlParam.PID.use = 0;

ctrl = BuiCtrl(model, CtrlParam);

%% Simulate
SimParam.run.start = 1;
SimParam.run.end = 13; 
SimParam.verbose = 1;
SimParam.flagSave = 0;
SimParam.comfortTol = 1e-1;
% flag distinguishing emulation and real measurements
%  0 - measurements
%  1 - emulation
SimParam.emulate = 1;
SimParam.profile = 0;  % profiler function for CPU evaluation

PlotParam.flagPlot = 1;     % plot 0 - no 1 - yes
PlotParam.plotStates = 0;        % plot states
PlotParam.plotDist = 0;        % plot disturbances
PlotParam.plotEstim = 0;        % plot estimation
PlotParam.plotCtrl = 1;        % plot control
PlotParam.plotPrice = 1;        % plot price signal
% PlotParam.Transitions = 1;      % pot dynamic transitions of Ax matrix
% PlotParam.reduced = 0;   %  reduced paper plots formats 0 - no 1 - yes
% PlotParam.zone = 2;     % choose zone if reduced
% PlotParam.only_zone = 0;    %  plot only zone temperatures 0 - no 1 - yes  

% %  simulation file with embedded plotting file
outdata = BuiSim(model, estim, ctrl, dist, refs, SimParam, PlotParam);



 
 
 