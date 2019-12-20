
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

%% Model: emulator + prediction

% =========== 1, choose building model =================
% select from a library of available models 
% buildingType = ModelIdentifier 
% ModelIdentifier for residential houses with radiators:   'Reno', 'Old', 'RenoLight'
% ModelIdentifier for office buildings with TABS:          'Infrax', 'HollandschHuys'
% ModelIdentifier for borehole:                            'Borehole'  - % TODO: missing disturbances precomputed file for borehole 
buildingType = 'Reno';  

% =========== 2, choose model order =================
ModelParam.Orders.range = [4, 7, 10, 15, 20, 30, 40, 100];    % suggested = model orders for 'Reno', 'Old', 'RenoLight'
% ModelParam.Orders.range = [100, 200, 600];                  % suggested model orders for 'Infrax', 'HollandschHuys'
ModelParam.Orders.choice = 'full';                            % model order selection for prediction
ModelParam.off_free = 1;                                      % augmented model with unmeasured disturbances
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
RefsParam.Price.variable = 0;       %1 =  variable price profile, 0 = fixed to 1

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
CtrlParam.RBC.use = 0;
CtrlParam.PID.use = 0;
CtrlParam.MLagent.use = 0;

ctrl = BeCtrl(model, CtrlParam);       % construct a controller object  

%% Simulate
SimParam.run.start = 11;
SimParam.run.end = 17; 
SimParam.verbose = 1;
SimParam.flagSave = 0;
SimParam.comfortTol = 1e-1;
SimParam.emulate = 1;  % emulation or real measurements:  0 = measurements,  1 = emulation
SimParam.profile = 0;  % profiler function for CPU evaluation

% %  simulation file with embedded plotting file
outdata = BeSim(model, estim, ctrl, dist, refs, SimParam);


%% Diagnose the MPC problem via Yalmip optimize

diagnoseFlag = true;
if diagnoseFlag
    % solve single instance of the MPC problem via Yalmip optimize
    [diagnostics, con, obj, outdata.con_info] = BeMPC_DualCheck(outdata, model)
end

%% Plot Results
PlotParam.flagPlot = 1;          % plot 0 - no 1 - yes
PlotParam.plotStates = 0;        % plot states
PlotParam.plotStates3D = 0;      % ribbon plot states
PlotParam.plotDist = 0;          % plot disturbances
PlotParam.plotDist3D = 0;        % ribbon plot disturbances
PlotParam.plotEstim = 0;         % plot estimation
PlotParam.plotEstim3D = 0;       % ribbon plot estimation
PlotParam.plotCtrl = 1;          % plot control
PlotParam.plotPrimalDual = 1;          % plot primal and dual varibles
PlotParam.plotPrimalDual3D = 1;        % ribbon plot primal and dual varibles
PlotParam.plotDualActive = 1;     % activation of the dual varibles
PlotParam.plotPrice = 0;          % plot price signal
% PlotParam.Transitions = 1;      % pot dynamic transitions of Ax matrix
% PlotParam.reduced = 0;   %  reduced paper plots formats 0 - no 1 - yes
% PlotParam.zone = 2;     % choose zone if reduced
% PlotParam.only_zone = 0;    %  plot only zone temperatures 0 - no 1 - yes  

if PlotParam.flagPlot
    BePlot(outdata,PlotParam)
end

%% Save Results
SaveParam.path = ['../Data/Simulations/',buildingType]; % savepath
SaveParam.save = false;                     % save or not
SaveParam.data.states = false;              % X      
SaveParam.data.outputs = true;              % Y    
SaveParam.data.inputs = true;               % U   
SaveParam.data.disturbances = true;         % D 
SaveParam.data.references = false;          % WA, WB 
SaveParam.solver.objective = true;          % objective values of the QP optimization problem
SaveParam.solver.duals = true;              % dual variables of the QP optimization problem
SaveParam.solver.primals = false;           % primal variables of the QP optimization problem
SaveParam.solver.SolverTime = true;         % solvertime
SaveParam.solver.iters = true;              % solver iterations
SaveParam.solver.specifics = false;         % solver specific information

if SaveParam.save
    BeSave(outdata,SaveParam)
end







 