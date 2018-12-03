function outdata = BuiSim(model, estim, ctrl, dist, refs, SimParam, PlotParam)

if nargin == 0  % model
   buildingType = 'Infrax';  
   ModelOrders.range = [100, 200, 600]; % add any reduced order model you wish to have
%    ModelOrders.choice = 200;            % insert model order or 'full' for full order SSM 
   ModelOrders.choice = 100;            % insert model order or 'full' for full order SSM 
   ModelOrders.off_free = 0;            %  augmented model
   reload = 0;
   model = BuiModel(buildingType, ModelOrders, reload);     %  construct the model
end
if nargin < 2   % estimator
   EstimParam.LOPP.use = 0;      %  Luenberger observer via pole placement 
   EstimParam.SKF.use = 0;    % stationary KF
   EstimParam.TVKF.use = 1;   % time varying KF
   EstimParam.MHE.use = 0;   % moving horizon estimation via yalmip
   EstimParam.MHE.Condensing = 1;   % moving horizon estimation via yalmip
   estim = BuiEstim(model, EstimParam);
end
if nargin < 3   % controller
   CtrlParam.precomputed = 1;
   CtrlParam.MPC.use = 0;
   CtrlParam.MPC.Condensing = 1;
   CtrlParam.RBC.use = 0;
   CtrlParam.PID.use = 0;
   ctrl = BuiCtrl(model, CtrlParam);
end
if nargin < 4  % disturbances
   dist = BuiDist(model.buildingType, model.reload);
end
if nargin < 5   % references
   refs = BuiRefs();      % TODO FIX to be general
end
if nargin < 6   % simulation parameters
    SimParam.run.start = 1;     % starting day
%     SimParam.run.end = 13;       % finishing day
    SimParam.run.end = 3;       % finishing day
    SimParam.verbose = 1;
    SimParam.flagSave = 0;
    SimParam.comfortTol = 1e-1;
    SimParam.profile = 0; 
end
if nargin < 7   % plotting  
    PlotParam.flagPlot = 0;     % plot 0 - no 1 - yes
    PlotParam.plotStates = 0;        % plot states
    PlotParam.plotDist = 1;        % plot disturbances
    PlotParam.plotEstim = 1;        % plot estimation
    PlotParam.plotCtrl = 1;        % plot control
    PlotParam.plotPrice = 1;        % plot price signal
end

% matlab profiler function for CPU evaluation
if SimParam.profile
    profile on 
end

%% Simulation setup   
fprintf('\n------------------ Simulation Setup -----------------\n');

fprintf('*** Building Type = %s\n' , model.buildingType);
fprintf('*** Prediction model order = %d, \n',   size(model.pred.Ad,1))
fprintf('*** Start day = %d , End day = %d \n', SimParam.run.start, SimParam.run.end);

% Simulation steps   
% starting and finishing  second:  24h = 86400 sec
SimStart_sec = (SimParam.run.start-1)*86400;
SimStop_sec = (SimParam.run.end)*86400;
% starting and finishing  step for simulation loop - MPC
SimStart = floor(SimStart_sec/model.plant.Ts)+1;
SimStop = ceil(SimStop_sec/model.plant.Ts);
% number of simulation steps for MPC
Nsim = length(SimStart:SimStop);

%% Initial values  

% if not(exist('ctrl.MLagent'))
%     ctrl.MLagent.use = 0;
% end

% preview setup
if ctrl.MPC.use
            N = ctrl.MPC.N;
            Nrp = ctrl.MPC.Nrp;
elseif ctrl.MLagent.use
            N = ctrl.MLagent.numDelays;
            Nrp = ctrl.MLagent.numDelays;
else
            N = 0;
            Nrp = 0;
end


X = zeros(model.plant.nx,Nsim+1);
D = dist.d(SimStart:SimStop+N,:)';

if  not(ctrl.use)  % precomputed inputs and outputs
    U = ctrl.precomputed.U(:,SimStart:SimStop);
%     Y = ctrl.precomputed.Y(:,SimStart:SimStop)+273.15;   % infrax - signal in deg C
    Y = ctrl.precomputed.Y(:,SimStart:SimStop);    
    uopt = 0*U(:,1); % initialize controls
    
else   % initialize matrices for closed loop control simulations
    Y = zeros(model.plant.ny,Nsim)+model.plant.Fd*1;
    U = zeros(model.plant.nu,Nsim);
    
    uopt = U(:,1); % initialize controls

    % ------ references ------
    R = refs.R(SimStart:SimStop,:)';
    % above and below threshold comfort zone
    wa = refs.wa(SimStart:SimStop+N,:)';
    wb = refs.wb(SimStart:SimStop+N,:)';
    % ------ PMV zone ------ 
    PMVub = refs.PMVub(SimStart:SimStop);
    PMVlb = refs.PMVlb(SimStart:SimStop);
    % ------ energy price ------
    Price = refs.Price(SimStart:SimStop+N,:)';

    if ctrl.RBC.use
        % supply water temperature
        TSup = refs.TSup(SimStart:SimStop,:)';
        heat = 0;
    end

end

if estim.use 
    %     estmator vector inits
    Xp = zeros(model.pred.nx,Nsim+1);
    Xe = zeros(model.pred.nx,Nsim);
    Ye = zeros(model.pred.ny,Nsim);
    Yp = zeros(model.pred.ny,Nsim);
%     Ye = 295.15*ones(model.pred.ny,Nsim);
%     Yp = 295.15*ones(model.pred.ny,Nsim);
%     Ye = 21.15*ones(model.pred.ny,Nsim);
%     Yp = 21.15*ones(model.pred.ny,Nsim);
    
    % current estim states
    xe = Xe(:, 1);  
    xp = Xp(:, 1); 
    ye = Ye(:, 1);  
    yp = Yp(:, 1); 
    
    if estim.TVKF.use
        EstimGain = cell(1,Nsim);
        ErrorCovar = cell(1,Nsim); 
    end

    if estim.MHE.use
        OBJ_MHE = zeros(1,Nsim);
        We = zeros(model.pred.nx,Nsim);
        Ve = zeros(model.pred.ny,Nsim);    
    end   
end


%     violation vectors
Viol = zeros(model.plant.ny,Nsim); 
AboveViol = zeros(model.plant.ny,Nsim);
BelowViol = zeros(model.plant.ny,Nsim); 
%     PMV vectors
PMV = zeros(model.plant.ny,Nsim); 
PMVViol =  zeros(model.plant.ny,Nsim); 
PMVAboveViol =  zeros(model.plant.ny,Nsim); 
PMVBelowViol =  zeros(model.plant.ny,Nsim); 
%  MPC objective function vectors
J = zeros(1,Nsim);     
J_v =  zeros(1,Nsim);   
J_u =  zeros(1,Nsim);   
OBJ =  zeros(1,Nsim);
   

%% ------ MAIN simulation loop ------
    % % ------ Verbose ------
fprintf('\n---------------- Simulation Running -----------------');
% % initiation clearing string 
reverseStr = '';

start_t = clock;

for k = 1:Nsim
    
%     current states, inputs, outputs and disturnances
    x0 = X(:,k);         % current sim states    - initialized to 0 at k = 1
    d0 = D(:,k);         % current disturbances  - from measurements    
        
    
%     TODO: implement 3 cases: 
% 1, plant simulation and control - all computed                 - DONE
% 2, measured (fixed) u and y                                    - DONE
% 3, measured (fixed) y - computed control and estimation        - TODO


    
%%  CONTROL  
% TODO standalone functions       
    if ctrl.use   
        if ctrl.RBC.use
            %  heat curve control 
            [uopt, heat] = BuiRBC(yn,R(:,k),heat,TSup(k),ctrl);
%             TODO - automatic tuning general case         
            
        elseif ctrl.PID.use
%             TODO - implementation
            
        elseif ctrl.MPC.use          
%             TODO: predictions as standalone function
            % preview of disturbance signals on the prediction horizon
            Dpreview = D(:, k:k+(ctrl.MPC.Ndp-1));   
            % preview of thresholds on the prediction horizon - Dynamic comfort zone
            wa_prev = wa(:, k:k+(ctrl.MPC.Nrp-1));
            wb_prev = wb(:, k:k+(ctrl.MPC.Nrp-1));
            % preview of the price signal
            Price_prev = Price(:, k:k+(ctrl.MPC.Nrp-1));
            
            
%             TODO:  adapt Dpreview wa_prev wb_prev
            if estim.use   % estimated states
                if model.plant.nd == 0  %  no disturbnances option
                     [opt_out, feasible, info1, info2] =  ctrl.MPC.optimizer{{xp, wa_prev, wb_prev, Price_prev}}; % optimizer with estimated states
                else
                     [opt_out, feasible, info1, info2] =  ctrl.MPC.optimizer{{xp, Dpreview, wa_prev, wb_prev, Price_prev}}; % optimizer with estimated states
                end
            else    % perfect state update
                if model.plant.nd == 0  %  no disturbnances option
                     [opt_out, feasible, info1, info2] =  ctrl.MPC.optimizer{{x0, wa_prev, wb_prev, Price_prev}}; % optimizer with estimated states
                else
                     [opt_out, feasible, info1, info2] =  ctrl.MPC.optimizer{{x0, Dpreview, wa_prev, wb_prev, Price_prev}}; % optimizer with measured states  
                end
                 
            end
            
%             %     feasibility check
%             if feasible == 3 %    3 Maximum iterations exceeded
% %                 uopt = does not change value
%                 obj = 0;
            if ~ismember(feasible, [0 3 4 5])
                k
                error('infeasible')      
            else
                uopt = opt_out{1};   % optimal control action
                obj =  opt_out{2};   %objective function value     
            end           
            
%             % Objective function value increments in each sim. step
%             J_v(k) = vb'*Qsb*vb + va'*Qsa*va;      % violations increments in objective
%             J_u(k) = uopt'*Qu*uopt;                % inputs increments in objective
%             J(k) =  J_v(k)+J_u(k);               % overalll increments in objective
%             %     yalmip objective function increments vector
%             OBJ(k) = obj;      


        
        elseif ctrl.MLagent.use   % machine learning controller
%             TODO: finish implementation of all ML ctrls
            if ctrl.MLagent.TDNN.use    % Time delay neural network           
                InputDelays = ctrl.MLagent.TDNN.net_apply.numInputDelays; % N-1
                
                if k > InputDelays
% features for TDNN             
%                     Dpreview = D(ctrl.MLagent.use_D, k:k+InputDelays);   
%                     wa_prev = wa(1, k:k+InputDelays);
%                     wb_prev = wb(1, k:k+InputDelays);
%                     y_past = Y(:,k-InputDelays:k-1);
%                      uopt =  TDNN_ctrl([yn; Dpreview(:,end); wb_prev(:,end)],...                              
%                      [ y_past; Dpreview(:,1:end-1); wb_prev(:,1:end-1)]);
%                      

%   current code - Working
                     uopt =  TDNN_ctrl([yn; D(ctrl.MLagent.use_D, k+InputDelays);wb(1,k+InputDelays)],...                              
                     [ Y(:,k-InputDelays:k-1); D(ctrl.MLagent.use_D, k:k+(InputDelays-1));wb(1,k:k+(InputDelays-1))]);

% previous code - WORKING
%                      uopt =  NN_TS_ctrl([yn; D(ctrl.MLagent.use_D, k+InputDelays);wb(1,k+InputDelays)],...                              
%                      [ Y(:,k-InputDelays:k-1); D(ctrl.MLagent.use_D, k:k+(InputDelays-1));wb(1,k:k+(InputDelays-1))]);

                    
                    
                else
                    uopt = zeros(model.pred.nu,1);
                end
                
            elseif ctrl.MLagent.RT.use
    
                
            end

     
        end
    end  
   
    
%% Controlled Sytem Dynamics
% TODO: wrap simulation and measurement in standalone function with standardized interface
    
% 1, Pre-computed controls
    if not(ctrl.use)
       uopt = U(:,k);       % current controls 
    end
         
% 2, EMULATOR - plant model
    if  SimParam.emulate
%    State and Output update
        xn = model.plant.Ad*x0 + model.plant.Bd*uopt+ model.plant.Ed*d0 +model.plant.Gd*1;
        yn = model.plant.Cd*x0 + model.plant.Dd*uopt + model.plant.Fd*1;

        % simulation model data vectors
        X(:,k+1) = xn;
        Y(:,k) = yn;
        U(:,k) = uopt;       
    end   
%     TODO: dymola co-simulation
    
    
% 3, Output measurements - real plant
    % TODO: implement real time measurement
    if  not(SimParam.emulate) 
       yn = Y(:,k);         % current outputs               
    end


%%   ESTIMATION  
% TODO standalone functions?
    if estim.use 
        if estim.SKF.use  % stationary KF
            
            % measurement update                              
            yp = model.pred.Cd*xp + model.pred.Dd*uopt + model.pred.Fd*1;          % output estimation
            ep = yn - yp;                                                       % estimation error
            xe = xp  + estim.SKF.L1*ep;                                       % estimated state
            
            % time update
            xp = model.pred.Ad*xe + model.pred.Bd*uopt + model.pred.Ed*d0 + model.pred.Gd*1;
            
            ye = model.pred.Cd*xe + model.pred.Dd*uopt + model.pred.Fd*1;     % output estimate with x[n|n]
             
        elseif estim.TVKF.use  % time varying KF           
            
              if k == 1
                  P = model.pred.Bd*estim.TVKF.Qe*model.pred.Bd';         % Initial error covariance   
              end
            
              % Measurement update
              L1 = P*model.pred.Cd'/(model.pred.Cd*P*model.pred.Cd'+estim.TVKF.Re); % observer gain
              yp = model.pred.Cd*xp + model.pred.Dd*uopt + model.pred.Fd*1;          % output estimation
              ep = yn - yp;                                                       % estimation error
              xe = xp + L1*ep;                                                    % x[n|n]
              P = (eye(model.pred.nx)-L1*model.pred.Cd)*P;                          % P[n|n]   estimation error covariance
              errcov = model.pred.Cd*P*model.pred.Cd';                              % output estimation error covariance
              
              % Time update
              xp = model.pred.Ad*xe + model.pred.Bd*uopt + model.pred.Ed*d0 + model.pred.Gd*1;        % x[n+1|n]
              P = model.pred.Ad*P*model.pred.Ad' + model.pred.Bd*estim.TVKF.Qe*model.pred.Bd';       % P[n+1|n]
            
              ye = model.pred.Cd*xe + model.pred.Dd*uopt + model.pred.Fd*1;     % output estimate with x[n|n]
              
              % time varying parameters data
              EstimGain{k} = L1;
              ErrorCovar{k} = errcov;
              
        elseif estim.MHE.use  % TODO: moving horizon estimation
   
            if k >= estim.MHE.N
             
                N = estim.MHE.N;
                                
                [opt_out, feasible, info1, info2] = estim.MHE.optimizer{{Y(:,k-N+1:k), U(:,k-N+1:k), D(:,k-N+1:k), Xp(:,k-N+1)}}; % optimizer with estimated states
               
                xe = opt_out{1};    % estimated state at x[n-N+1|n]
                ve =  opt_out{2};   % v decision variables at  x[n-N+1:n|n]
                obj_estim =  opt_out{3}; 
                
                if estim.MHE.Condensing                                       
                    we = zeros(model.pred.nx,estim.MHE.N);
                else
                    we =  opt_out{4};   % w decision variables at  x[n-N+1:n|n]
                end                          
                
            %    MHE post processing, integration of states at x[n-N+1|n] via state
            %    update and w to get x[n|n]           
                for j = 1:N
                    xe = model.pred.Ad*xe + model.pred.Bd*U(:,k-N+j) + model.pred.Ed*D(:,k-N+j) + model.pred.Gd*1 + we(:,j);
                end
             
                ye = model.pred.Cd*xe + model.pred.Dd*uopt + model.pred.Fd*1 + ve(:,N);     % output estimate with x[n|n]               
                
                OBJ_MHE(:,k) = obj_estim;
                We(:,k) = we(:,N);
                Ve(:,k) = ve(:,N);
            else
%              TODO:   put KF or growing horizon implementation for initial
%              estimate
                                                      
            end     
        end

    %     estimator data
        Xp(:,k+1) = xp;
        Xe(:,k) = xe;
        Ye(:,k) = ye;
        Yp(:,k) = yp;
    end
    
    
%% Control action post-processing
% TODO - heat flows to valve positions

% Q = m*cp*p*(T_sup - T_return)
% Q - 12 heat flows to zones computed by MPC
% m - prescribed nominal mass flow
% cp - thermal capcity of the water known
% p - all valve postition in pipe system to be computed by post-processing
% T_sup - supply water temperature to be measured
% T_return - return temperature to be measured
% % ALARM - CVUT used instead of T_return, T_concrete - why?

% 1, select only heat flows to zones as controls
% 2, compute necessary mass flows to individual zones based on given measurements
% 3, compute valve positions based on the mass balance of flows in whole pipe system  
% 4, consider variable mass flow rates
% 5, automate the process based on pipe topology

% problem: control input mismatch and measurements accuracy - minimized by
% extra state observer or lower level controller?

    
    
%% ---------- Comfort evaluation ----------

% move this to post-processing
% adapt computations to vector format to avoid loops

if  ctrl.use
    % VIOLATIONS of the thermal comfort zones of individual outputs
    va = zeros(model.plant.ny,1); vb = zeros(model.plant.ny,1); v = zeros(model.plant.ny,1);
    % VIOLATIONS of the PMV zone of individual outputs
    PMVva = zeros(model.plant.ny,1); PMVvb = zeros(model.plant.ny,1); PMVv = zeros(model.plant.ny,1);
    
    %     violation evaluation for j-th output
    for j = 1:model.plant.ny
        if yn(j) > wa(j,k) + SimParam.comfortTol    % above viol. condition
            v(j) = yn(j) - wa(j,k);        % above viol. magnitude
            va(j) = v(j);
        elseif yn(j)  <  wb(j,k) - SimParam.comfortTol  % below viol. condition
            v(j) = yn(j) - wb(j,k);           % below viol. magnitude
            vb(j) = v(j);
        end
        
        % PMV index for each zone, with Tr = 29 deg C
        PMV(j,k) = pmv_iso(yn(j)-273.15, 29);
        if  PMV(j,k) > PMVub(k)    % above PMV viol. condition
            PMVv(j) = PMV(j,k) - PMVub(k);        % above viol. magnitude
            PMVva(j) = PMVv(j);
        elseif  PMV(j,k)  <  PMVlb(k)  % below PMV viol. condition
            PMVv(j) = PMV(j,k) - PMVlb(k);        % above viol. magnitude
            PMVvb(j) = PMVv(j);
        end
    end
    %  comfort zone violations vectors
    Viol(:,k) = v;
    AboveViol(:,k) = va;
    BelowViol(:,k) = vb;
    % PMV violations vectors
    PMVViol(:,k) =  PMVv;
    PMVAboveViol(:,k) =  PMVva;
    PMVBelowViol(:,k) =  PMVvb;
end
            
%     REMAINING simulation time computation
    step_time = etime(clock, start_t);                  %  elapsed time of one sim. step
    av_step_time = step_time/k;                         % average_step_time
    rem_sim_time = av_step_time*(Nsim-k);           % remaining_sim_time
     
    msg = sprintf('\n*** estimated remaining simulation time = %.2f sec \n',rem_sim_time);    % statement   
    fprintf([reverseStr, msg]);                                                 % print statement
    reverseStr = repmat(sprintf('\b'), 1, length(msg));                         % clear line   
        
end

%% ---------- Simulation Output Data  ------------
% TODO: modify to be general

Ts = model.plant.Ts;

if  ctrl.use
    % -------- ENERGY COSTS ----------
    Uheat = U(U>0);
    Ucool = U(U<0);
    Qheat = 1;   % heat coefficient
    Qcool = 1;   % cool coefficient
    %  ENERGY COSTS for individudal inputs
        for j = 1:model.plant.ny
            outdata.info.HeatingCost(j) = sum(Qheat*Uheat(j:model.plant.ny:size(Uheat,1)))*Ts/1000/3600;       % heating cost [kW hours]
            outdata.info.CoolingCost(j) = sum(abs(Qcool*Ucool(j:model.plant.ny:size(Ucool,1))))*Ts/1000/3600;  % cooling cost [kW hours]
            outdata.info.TotalCost(j) = outdata.info.HeatingCost(j)+outdata.info.CoolingCost(j);             % total cost [kW hours]
        end
    % Overall ENERGY COST
    outdata.info.OverallHeatingCost = sum(outdata.info.HeatingCost);                           % heating cost [kW hours]
    outdata.info.OverallCoolingCost = sum(outdata.info.CoolingCost);                           % cooling cost [kW hours]
    outdata.info.OverallTotalCost = sum(outdata.info.OverallHeatingCost)+sum(outdata.info.OverallCoolingCost);     % total cost [kW hours]
end

% ------------ COMFORT -----------------
% COMFORT satisfaction for each output
if  ctrl.use
    for j = 1:model.plant.ny
        %  comfort satisfacion
        outdata.info.PositiveComfortRate(j) = 100*(1-nnz(AboveViol(j,:))/length(AboveViol(j,:)));   %  PCR
        outdata.info.NegativeComfortRate(j) = 100*(1-nnz(BelowViol(j,:))/length(BelowViol(j,:)));   %  NCR
        outdata.info.TotalComfortRate(j) = 100*(1-nnz(Viol(j,:))/length(Viol(j,:)));                %  TCR
        %   sum of comfort zone violations
        outdata.info.PositiveViolSum(j) =  norm(AboveViol(j,:), 1);  % sum of positive viol - PVS
        outdata.info.NegativeViolSum(j) =  norm(BelowViol(j,:), 1);  % sum of negative viol - NVS
        outdata.info.TotalViolSum(j) =  norm(Viol(j,:), 1);          % sum of total viol - TVS
        %   maximum of comfort zone violation
        outdata.info.PositiveViolMax(j) =  norm(AboveViol(j,:), Inf); % max of positive viol - PVM
        outdata.info.NegativeViolMax(j) =  norm(BelowViol(j,:), Inf); % max of negative viol - NVM
        outdata.info.TotalViolMax(j) =  norm(Viol(j,:), Inf);         % max of total viol - TVM
        %  kelvin hours
        outdata.info.KelvinHours(j) =  norm(Viol(j,:), 1)*900/3600;          % Kh
         %   sum of PMV violations
        outdata.info.PMVPositiveViolSum(j) =  norm(PMVAboveViol(j,:), 1);  % sum of positive viol - PVS
        outdata.info.PMVNegativeViolSum(j) =  norm(PMVBelowViol(j,:), 1);  % sum of negative viol - NVS
        outdata.info.PMVTotalViolSum(j) =  norm(PMVViol(j,:), 1);          % sum of total viol - TVS
        %   maximum of PMV violation
        outdata.info.PMVPositiveViolMax(j) =  norm(PMVAboveViol(j,:), Inf); % max of positive viol - PVM
        outdata.info.PMVNegativeViolMax(j) =  norm(PMVBelowViol(j,:), Inf); % max of negative viol - NVM
        outdata.info.PMVTotalViolMax(j) =  norm(PMVViol(j,:), Inf);         % max of total viol - TVM
    end
% Overall COMFORT zone satisfaction
outdata.info.Overall_PCR = sum(outdata.info.PositiveComfortRate)/length(outdata.info.PositiveComfortRate);  % positive
outdata.info.Overall_NCR = sum(outdata.info.NegativeComfortRate)/length(outdata.info.NegativeComfortRate);  % negative
outdata.info.Overall_TCR = sum(outdata.info.TotalComfortRate)/length(outdata.info.TotalComfortRate);        % total
% Overall sum of comfort zone violations
outdata.info.Overall_PVS = sum(outdata.info.PositiveViolSum);  % sum of positive viol
outdata.info.Overall_NVS = sum(outdata.info.NegativeViolSum);  % sum of negative viol
outdata.info.Overall_TVS = sum(outdata.info.TotalViolSum);     % sum of total viol
% Overall maximum of comfort zone violations
outdata.info.Overall_PVM =  max(outdata.info.PositiveViolMax);  % max of positive viol
outdata.info.Overall_NVM =  max(outdata.info.NegativeViolMax);  % max of negative viol
outdata.info.Overall_TVM =  max(outdata.info.TotalViolMax);     % max of total viol
% K/hour metric per zone
outdata.info.Overall_PV_metric = outdata.info.Overall_PVS/(Nsim*(Ts/3600))/model.plant.ny;  % K/h metric of positive viol
outdata.info.Overall_NV_metric = outdata.info.Overall_NVS/(Nsim*(Ts/3600))/model.plant.ny;  % K/h metric of negative viol
outdata.info.Overall_TV_metric = outdata.info.Overall_TVS/(Nsim*(Ts/3600))/model.plant.ny;  % K/h metric of total viol
% Overall kelvin hours
outdata.info.Overall_Kh = sum(outdata.info.KelvinHours);     % sum of kelvin hours
% Overall sum of PMV violations
outdata.info.Overall_PMV_PVS = sum(outdata.info.PMVPositiveViolSum);  % sum of positive viol
outdata.info.Overall_PMV_NVS = sum(outdata.info.PMVNegativeViolSum);  % sum of negative viol
outdata.info.Overall_PMV_TVS = sum(outdata.info.PMVTotalViolSum);     % sum of total viol
% Overall maximum of PMV violations
outdata.info.Overall_PMV_PVM =  max(outdata.info.PMVPositiveViolMax);  % max of positive viol
outdata.info.Overall_PMV_NVM =  max(outdata.info.PMVNegativeViolMax);  % max of negative viol
outdata.info.Overall_PMV_TVM =  max(outdata.info.PMVTotalViolMax);     % max of total viol
end

% ------------ DATA -----------------

% SIMILATION STRUCTURES
outdata.model = model;    %  model 
outdata.estim = estim;    % estimator
outdata.ctrl = ctrl;      %  controller
% outdata.dist = dist;      %  disturbances
% outdata.dist = refs;      %  references
outdata.SimParam = SimParam;        %  simulation parameters
outdata.SimParam.run.Nsim = Nsim;    % sim steps

% plant simulation data 
outdata.data.X = X;         %  state vector
outdata.data.Y = Y;         %  output vector
outdata.data.U = U;         %  input vector
outdata.data.D = D;         %  disturbance vector

% estimator data
if estim.use ==1
    outdata.data.Xe = Xe;         %  estimated state vector  [n|n]
    outdata.data.Ye = Ye;         %  estimated output vector [n|n]
    outdata.data.Xp = Xp;         %  previous estimaror state vector [n|n-1]
    outdata.data.Yp = Yp;         %  estimated output vector [n|n-1]
    
    if estim.TVKF.use
        outdata.data.EstimGain = EstimGain;
        outdata.data.ErrorCovar = ErrorCovar;
        
    elseif estim.MHE.use
        outdata.data.We = We;         
        outdata.data.Ve = Ve;    
    end

end

% % TODO: integrate this
if ctrl.use
    % referecne + comfort zone
    % outdata.data.R = refs_mpc;  %  references vector
    outdata.data.wa = wa;       %  above threshold
    outdata.data.wb = wb;       %  below threshold
%     % PMV zone
%     outdata.data.PMVub = PMVub; %  above threshold
%     outdata.data.PMVlb = PMVlb; %  below threshold
%     %  comfot zone violations
%     outdata.data.AV = AboveViol;
%     outdata.data.BV = BelowViol;
%     outdata.data.V = Viol;
%     % PMV index profiles
%     outdata.data.PMV = PMV;
%     %  PMV zone violations
%     outdata.data.PMV_AV = PMVAboveViol;
%     outdata.data.PMV_BV = PMVBelowViol;
%     outdata.data.PMV_V = PMVViol;
    
%     Price signal
    outdata.data.Price = Price(:,1:end-Nrp);       
    outdata.data.Cost = Price(:,1:end-Nrp).*U;   
    
%     if ctrl.MPC.use
%         % obj function
%         outdata.data.J = J;
%         outdata.data.J_v = J_v;
%         outdata.data.J_u = J_u;
%         outdata.data.OBJ = OBJ;
%     end
end

%  elapsed time of simulation
outdata.info.cmp = etime(clock, start_t);

%% Simulation results reports and plots

% print verbose condition
if SimParam.verbose
        % -------- PRINT --------------
        fprintf('\n------------------ Simulation Results ---------------\n');
        
    if ctrl.use
        %  energy cost
        fprintf('          Heating cost: %.2f kWh\n', outdata.info.OverallHeatingCost);
        fprintf('          Cooling cost: %.2f kWh\n', outdata.info.OverallCoolingCost);
        fprintf('            Total cost: %.2f kWh\n', outdata.info.OverallTotalCost);
        fprintf('               Comfort: %.2f Kh\n',  outdata.info.Overall_Kh);
        fprintf('        PMV violations: %.2f \n', outdata.info.Overall_PMV_TVS);
    end
    
        % compuation time
        fprintf('*** Simulation Time: %.1f secs\n', outdata.info.cmp);
    
end
% Simulation ends  
outdata.info.date =  datestr(now);
fprintf('*** Simulation finished at: %s \n',  outdata.info.date);

% SAVE data 
if SimParam.flagSave
    str = sprintf('../Data/outData%s_from%d_to%d.mat', model.buildingType, SimParam.run.start, SimParam.run.end);
    save(str,'outdata');
end

% PLOT the outut data 
if PlotParam.flagPlot
    BuiPlot(outdata,PlotParam)
end

% PROFILE CPU load
if SimParam.profile
    outdata.profile = profile('info');
    profile viewer
end


end