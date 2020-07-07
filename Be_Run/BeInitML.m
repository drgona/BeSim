
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

%% Model: emulator + prediction

% =========== 1, choose building model =================
% select from a library of available models 
% buildingType = ModelIdentifier 
% ModelIdentifier for residential houses with radiators:   'Reno', 'Old', 'RenoLight'
% ModelIdentifier for office buildings with TABS:          'Infrax', 'HollandschHuys'
% ModelIdentifier for borehole:                            'Borehole' 

buildingType = 'HollandschHuys';  
% buildingType = 'HollandschHuys';  


% =========== 2, choose model order =================
% ModelParam.Orders.range = [4, 7, 10, 15, 20, 30, 40, 100];    % suggested = model orders for 'Reno', 'Old', 'RenoLight'
ModelParam.Orders.range = [100, 200, 600];                  % suggested model orders for 'Infrax', 'HollandschHuys'
ModelParam.Orders.choice = 100;                            % model order selection for prediction
ModelParam.off_free = 0;                                      % augmented model with unmeasured disturbances
ModelParam.reload = 0;                                        % if 1 reload ROM, if 0 load saved ROM

% =========== 4, choose model analysis =================
ModelParam.analyze.SimSteps = 2*672; % Number of simulation steps (Ts = 900 s),  672 = one week
ModelParam.analyze.openLoop.use = false;             %  open loop simulation   - TODO
ModelParam.analyze.openLoop.start = 1;              % starting day of the analysis
ModelParam.analyze.openLoop.end = 7;                % ending day of the analysis
ModelParam.analyze.nStepAhead.use = false;           % n-step ahead predicion error  - TODO
ModelParam.analyze.nStepAhead.steps = [1, 10, 40];  % x*Ts  
ModelParam.analyze.HSV = false;                      %  hankel singular values of ROM
ModelParam.analyze.frequency = false;                % frequency analysis - TODO

% =========== 4, construct model structue =================
model = BeModel(buildingType, ModelParam);      % construct a model object   


%% Disturbacnes 
% ambient temperature, solar radiation, internal heat gains
DistParam.reload = 0;

dist = BeDist(model, DistParam);        % construct a disturbances object  

%% References 
% comfort constraints, price profiles
RefsParam.Price.variable = 1;       %1 =  variable price profile, 0 = fixed to 1

refs = BeRefs(model, dist, RefsParam);     % construct a references object  

%% Estimator 
EstimParam.SKF.use = 0;          % stationary KF
EstimParam.TVKF.use = 1;         % time varying KF
EstimParam.MHE.use = 0;          % moving horizon estimation via yalmip
EstimParam.MHE.Condensing = 1;   % state condensing 
EstimParam.use = 1;

estim = BeEstim(model, EstimParam);      % construct an estimator object  

%% Controller 
CtrlParam.use = 1;   % 0 for precomputed u,y    1 for closed loop control
CtrlParam.MPC.use = 1;
CtrlParam.MPC.Condensing = 1;
CtrlParam.LaserMPC.use = 0;
CtrlParam.LaserMPC.Condensing = 1;
CtrlParam.RBC.use = 0;
CtrlParam.PID.use = 0;
CtrlParam.MLagent.use = 0;

ctrl = BeCtrl(model, CtrlParam);       % construct a controller object  

%% Simulate - MPC train set
SimParam.run.start = 1;
SimParam.run.end = 2; 
SimParam.verbose = 1;
SimParam.flagSave = 0;
SimParam.comfortTol = 1e-1;
SimParam.emulate = 1;  % emulation or real measurements:  0 = measurements,  1 = emulation
SimParam.profile = 0;  % profiler function for CPU evaluation

% %  simulation file with embedded plotting file
outdata = BeSim(model, estim, ctrl, dist, refs, SimParam);



%% Diagnose the MPC problem via Yalmip optimize
DiagnoseParam.diagnoseFlag = 0;
DiagnoseParam.Duals.plotCheck = 0;
DiagnoseParam.Reduce.lincols.use = 1;
DiagnoseParam.Reduce.PCA.use = 1;
DiagnoseParam.Reduce.PCA.normalize = 0;             % normalize constraints based on types
DiagnoseParam.Reduce.PCA.component = 0.999999999;   % principal component weight threshold
DiagnoseParam.Reduce.PCA.feature = 0.999999999;     % PCA features weight threshold
if DiagnoseParam.diagnoseFlag
    % solve single instance of the MPC problem via Yalmip optimize
    [diagnostics, con, obj, outdata.con_info] = BeMPC_DualCheck(outdata, model, DiagnoseParam);
end

%% Plot Results
PlotParam.flagPlot = true;          % plot 0 - no 1 - yes
PlotParam.plotStates = 1;        % plot states
PlotParam.plotStates3D = 0;      % ribbon plot states
PlotParam.plotDist = 0;          % plot disturbances
PlotParam.plotDist3D = 0;        % ribbon plot disturbances
PlotParam.plotEstim = 0;         % plot estimation
PlotParam.plotEstim3D = 0;       % ribbon plot estimation
PlotParam.plotCtrl = 1;          % plot control
PlotParam.plotPrice = 1;          % plot price signal
if DiagnoseParam.diagnoseFlag
    PlotParam.plotPrimalDual = 1;          % plot primal and dual varibles
    PlotParam.plotPrimalDual3D = 1;        % ribbon plot primal and dual varibles
    PlotParam.plotDualActive = 1;     % activation of the dual varibles
    PlotParam.plotPCA_Dual = 1;        % principal components of PCA reduced dual variables   
    PlotParam.plotActiveSet = 1;         % plot active sets
else
    PlotParam.plotPrimalDual = 0;          % plot primal and dual varibles
    PlotParam.plotPrimalDual3D = 0;        % ribbon plot primal and dual varibles
    PlotParam.plotDualActive = 0;     % activation of the dual varibles
    PlotParam.plotPCA_Dual = 0;        % principal components of PCA reduced dual variables 
    PlotParam.plotActiveSet = 0;     % plot active sets
end
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

%% Extending MPC training dataset 
ExtraSampling = false;
if ExtraSampling
    % increasing state inititalization range for generalization purposes
    SampleParam.use = 1;   % use dataset extension 
    % Methods:
    SampleParam.StateRandom = 0; % 1, random sampling along the optimal trajectory - TODO
    SampleParam.StateInit = 1; % 2,  range of state initialization per day
    PlotParam.flagPlot = 0;     % plot 0 - no 1 - yes
    SampleParam.continue.use = false;

    % extended dataset
    outdataSampled = BeSample(outdata,SampleParam,PlotParam);

    if true
        % continue sampling of outdataSampled from a given day i and sample k
        SampleParam.use = 1;   % use dataset extension 
        SampleParam.StateRandom = 0; % 1, random sampling along the optimal trajectory - TODO
        SampleParam.StateInit = 1; % 2,  range of state initialization per day
        PlotParam.flagPlot = 0;     % plot 0 - no 1 - yes
        SampleParam.continue.use = true;
        SampleParam.continue.i = 123;
        SampleParam.continue.k = 8;
        outdataSampled.ctrl = ctrl;
        outdataSampled = BeSample(outdataSampled,SampleParam,PlotParam);
    end

else
    outdataSampled = outdata;
end

%% ====== Machine Learning Agent ======

% pre-defined ML agent models
AgentParam.RT.use = 0;       % regression tree with orthogonal splits - TODO
AgentParam.regNN.use = 0;    % function fitting regression neural network - TODO
AgentParam.TDNN.use = 1;     % time delayed neural network for time series approximation - OK
AgentParam.custom.use = 0;   % custom designed ML model - TODO

% initialize MLagent type and structure
MLagent = BeMLAgent(outdataSampled, AgentParam);

%% ====== Features ======
FeaturesParam.time_transform = 0;  % time transformations suitable for R
% feature refuction parameters
FeaturesParam.reduce.PCA.use = 1;
FeaturesParam.reduce.PCA.component = 0.999;   % principal component weight threshold
FeaturesParam.reduce.PCA.feature = 0.999;      % PCA features weight threshold
FeaturesParam.reduce.D_model.use = 1;
FeaturesParam.reduce.D_model.feature = 0.999;   % model features weight threshold
FeaturesParam.reduce.lincols.use = 1;
FeaturesParam.reduce.flagPlot = 1;
% generate training data {traindata} for a given agent from simulation data {outdata}
[traindata, MLagent] = BeFeatures(outdataSampled, dist, MLagent, FeaturesParam);


%% ====== MLagent training ======
% train ML agent
Meta_Epochs = 5;
for jj = 1:Meta_Epochs
    MLagent = BeLearn(MLagent,traindata);
%     save MLagent_TDNN2.mat MLagent
%     pause(900)
end


%% ====== MLagent Test Set Simulation ======
ctrl.use = 1;
ctrl.MLagent = MLagent;
ctrl.MLagent.use = 1;
ctrl.MPC.use = 0;
ctrl.RBC.use = 0;
ctrl.PID.use = 0;
SimParam.run.start = 1;
% SimParam.run.end = 30; 
SimParam.run.end = 330; 
% SimParam.run.start = 331;
% SimParam.run.end = 360; 
MLoutdata = BeSim(model, estim, ctrl, dist, refs, SimParam);

PlotParam.flagPlot = 1;     % plot 0 - no 1 - yes
if PlotParam.flagPlot
    BePlot(MLoutdata,PlotParam)
end

%% ====== MPC Test Set Simulation ======
ctrl.use = 1;
ctrl.MLagent.use = 0;
ctrl.MPC.use = 1;
MPC_test_outdata = BeSim(model, estim, ctrl, dist, refs, SimParam);
if PlotParam.flagPlot
    BePlot(MPC_test_outdata,PlotParam)
end



 
 
 