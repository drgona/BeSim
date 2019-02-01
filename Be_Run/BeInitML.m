
%% BeSim 
% Matlab toolbox for fast developlent, simulation and deployment of
% advanced building climate controllers 

% functionality intended for automatic construction of controls and
% estimation for a given linear building model

% Main control strategies
% 1, Model Predictive Control (MPC) 
% 2, deep learning control supervised by MPC


addpath('../Be_Modeling/')
addpath('../Be_Disturbances/')
addpath('../Be_References/')
addpath('../Be_Estimation/')
addpath('../Be_Control/')
addpath('../Be_Simulation/')
addpath('../Be_Learn/')
addpath('../Be_RealTime/')
addpath('../Be_RealTime/Mervis')

%% Model: emulator + predictor (controller)
% =========== 1, choose building model =================
% == Option 1: load custom model     %%%% TODO  %%%%
% buildingType = 'Load'
% == Option 2: select from library of available models 
% buildingType = ModelIdentifier 
% ModelIdentifier for residential houses with radiators:   'Reno', 'Old', 'RenoLight'
% ModelIdentifier for office buildings with TABS:          'Infrax', 'HollandschHuys'
% ModelIdentifier for borehole:                            'Borehole' 
% TODO: missing disturbance file for borehole large file on github
% =========== 2, choose model order =================
% ModelParam.Orders.range  = [4, 10, 20, 40, 100, ... ]   % vector of model orders 
% ModelParam.Orders.choice = 100                          % particular model order selection  
% ModelParam.Orders.choice = 'full'                       % full model order selection  
% ModelParam.Orders.off_free = 0 or 1                     % augmented model
% =========== 3, construct model structue =================
% model = BeModel(buildingType, ModelParam);

buildingType = 'Reno';  ModelParam.Orders.range = [4, 7, 10, 15, 20, 30, 40, 100];
% buildingType = 'Infrax'; ModelParam.Orders.range = [100, 200, 600]; 
% buildingType = 'HollandschHuys'; ModelParam.Orders.range = [100, 200, 600]; 
% buildingType = 'Borehole';  ModelParam.Orders.range = [10, 15, 20, 40, 100];  
ModelParam.Orders.choice = 'full';
ModelParam.Orders.off_free = 0;    
ModelParam.reload = 0; 

model = BeModel(buildingType, ModelParam);      % construct a model object   

%% Constraints
% TODO: state, input, algebraic equations...
% const = BuiConstraints(model,ConstrParam)

%% Disturbacnes 
% ambient temperature, solar radiation, internal heat gains
% =========== 1, choose disturbances =================
% == Option 1: load custom data    %%%% TODO  %%%%
% == Option 2: select from library of available models 

DistParam.reload = 0;

dist = BeDist(model, DistParam);        % construct a disturbances object  

%% References 
% comfort constraints, price profiles
RefsParam.Price.variable = 0;       %1 =  variable price profile, 0 = fixed to 1

refs = BeRefs(model, RefsParam);     % construct a references object  

%%  estimator 
EstimParam.LOPP.use = 0;      %  Luenberger observer via pole placement - Not implemented
EstimParam.SKF.use = 0;    % stationary KF
EstimParam.TVKF.use = 1;   % time varying KF
EstimParam.MHE.use = 0;   % moving horizon estimation via yalmip
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
SimParam.run.end = 13; 
SimParam.verbose = 1;
SimParam.flagSave = 0;
SimParam.comfortTol = 1e-1;
SimParam.emulate = 1;  % emulation or real measurements:  0 = measurements,  1 = emulation
SimParam.profile = 0;  % profiler function for CPU evaluation

% %  simulation file with embedded plotting file
outdata = BeSim(model, estim, ctrl, dist, refs, SimParam);


%% Plot Results
PlotParam.flagPlot = 1;     % plot 0 - no 1 - yes
PlotParam.plotStates = 0;        % plot states
PlotParam.plotDist = 0;        % plot disturbances
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

%% ========================================================================
%% BuiInitML
% machine learning approximations of MPC

%% ====== Machine Learning Agent ======

% pre-defined ML agent models
AgentParam.RT.use = 0;       % regression tree with orthogonal splits - TODO
AgentParam.regNN.use = 0;    % function fitting regression neural network - TODO
AgentParam.TDNN.use = 1;     % time delayed neural network for time series approximation - OK
AgentParam.custom.use = 0;   % custom designed ML model - TODO

% initialize MLagent type and structure
MLagent = BeMLAgent(outdata, AgentParam);


%% ====== Features ======

FeaturesParam.time_transform = 0;  % time transformations suitable for R
% feature refuction parameters
FeaturesParam.reduce.PCA.use = 1;
FeaturesParam.reduce.PCA.component = 0.999;   % principal component weight threshold
FeaturesParam.reduce.PCA.feature = 0.95;      % PCA features weight threshold
FeaturesParam.reduce.D_model.use = 1;
FeaturesParam.reduce.D_model.feature = 0.99;   % model features weight threshold
FeaturesParam.reduce.lincols.use = 1;
FeaturesParam.reduce.flagPlot = 1;

% generate training data {traindata} for a given agent from simulation data {outdata}
[traindata, MLagent] = BeFeatures(outdata, dist, MLagent, FeaturesParam);


%% ====== MLagent training ======

% train ML agent
MLagent = BeLearn(MLagent,traindata);


%% ====== MLagent Simulation ======

ctrl.use = 1;
ctrl.MLagent = MLagent;
ctrl.MLagent.use = 1;
ctrl.MPC.use = 0;
ctrl.RBC.use = 0;
ctrl.PID.use = 0;

MLoutdata = BeSim(model, estim, ctrl, dist, refs, SimParam);

if PlotParam.flagPlot
    BePlot(MLoutdata,PlotParam)
end









 
 
 