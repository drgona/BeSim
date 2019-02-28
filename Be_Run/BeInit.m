
%% BeSim 
% Matlab toolbox for fast developlent, simulation and deployment of
% advanced building climate controllers 
% Model Predictive Control (MPC)
% functionality intended for automatic construction of controls and
% estimation for a given linear building model

clear
yalmip('clear');
% close all

addpath('../Be_Modeling/')
addpath('../Be_Disturbances/')
addpath('../Be_References/')
addpath('../Be_Estimation/')
addpath('../Be_Control/')
addpath('../Be_Simulation/')
addpath('../Be_Learn/')

%% Model: emulator + predictor (controller)

% =========== 1, choose building model =================
% == Option 1: select from library of available models 
% buildingType = ModelIdentifier 
% ModelIdentifier for residential houses with radiators:   'Reno', 'Old', 'RenoLight'
% ModelIdentifier for office buildings with TABS:          'Infrax', 'HollandschHuys'
% ModelIdentifier for borehole:                            'Borehole' 

% =========== 2, choose model order =================
% ModelParam.Orders.range  = [4, 10, 20, 40, 100, ... ]   % vector of model orders 
% ModelParam.Orders.choice = 100                          % particular model order selection  
% ModelParam.Orders.choice = 'full'                       % full model order selection  
% ModelParam.Orders.off_free = 0 or 1                     % augmented model

% =========== 3, construct model structue =================
% model = BeModel(buildingType, ModelParam);

% =========== particular example =================
% buildingType = 'Reno';  ModelParam.Orders.range = [4, 7, 10, 15, 20, 30, 40, 100];
buildingType = 'HollandschHuys'; ModelParam.Orders.range = [100, 200, 600]; 
ModelParam.Orders.choice = 'full';
ModelParam.Orders.off_free = 0;    
ModelParam.reload = 0; 

model = BeModel(buildingType, ModelParam);      % construct a model object   


%% Disturbacnes 
% ambient temperature, solar radiation, internal heat gains
DistParam.reload = 0;

dist = BeDist(model, DistParam);        % construct a disturbances object  

%% References 
% comfort constraints, price profiles
RefsParam.Price.variable = 0;       %1 =  variable price profile, 0 = fixed to 1

refs = BeRefs(model, dist, RefsParam);     % construct a references object  

%%  estimator 
EstimParam.SKF.use = 0;          % stationary KF
EstimParam.TVKF.use = 1;         % time varying KF
EstimParam.MHE.use = 0;          % moving horizon estimation via yalmip
EstimParam.MHE.Condensing = 1;   % state condensing 
EstimParam.use = 1;

estim = BeEstim(model, EstimParam);      % construct an estimator object  

%% controller 
CtrlParam.use = 1;   % 0 for precomputed u,y    1 for closed loop control
CtrlParam.MPC.use = 1;
CtrlParam.MPC.Condensing = 1;
CtrlParam.RBC.use = 0;
CtrlParam.PID.use = 0;
CtrlParam.MLagent.use = 0;

ctrl = BeCtrl(model, CtrlParam);       % construct a controller object  

%% Simulate
SimParam.run.start = 1;
SimParam.run.end = 7; 
SimParam.verbose = 1;
SimParam.flagSave = 0;
SimParam.comfortTol = 1e-1;
SimParam.emulate = 1;  % emulation or real measurements:  0 = measurements,  1 = emulation
SimParam.profile = 0;  % profiler function for CPU evaluation

% %  simulation file with embedded plotting file
outdata = BeSim(model, estim, ctrl, dist, refs, SimParam);


%% Plot Results
PlotParam.flagPlot = 1;     % plot 0 - no 1 - yes
PlotParam.plotStates = 1;        % plot states
PlotParam.plotDist = 1;        % plot disturbances
PlotParam.plotEstim = 1;        % plot estimation
PlotParam.plotCtrl = 1;        % plot control
PlotParam.plotPrice = 1;        % plot price signal
% PlotParam.Transitions = 1;      % pot dynamic transitions of Ax matrix
% PlotParam.reduced = 0;   %  reduced paper plots formats 0 - no 1 - yes
% PlotParam.zone = 2;     % choose zone if reduced
% PlotParam.only_zone = 0;    %  plot only zone temperatures 0 - no 1 - yes  

if PlotParam.flagPlot
    BePlot(outdata,PlotParam)
end
 
 
 