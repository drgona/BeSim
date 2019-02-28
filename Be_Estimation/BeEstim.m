function estimator = BeEstim(model, EstimParam)

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
   EstimParam.use = 1;       % use estimator if not take states from plant model
   EstimParam.LOPP.use = 0;  %  Luenberger observer via pole placement - Not implemented
   EstimParam.SKF.use = 0;    % stationary KF
   EstimParam.TVKF.use = 0;   % time varying KF
   EstimParam.MHE.use = 1;   % moving horizon estimation
   EstimParam.MHE.Condensing = 1;   % state condensing 
end

% estimator parameters
estimator.use = EstimParam.use;
estimator.SKF.use = EstimParam.SKF.use;
estimator.TVKF.use = EstimParam.TVKF.use;
estimator.MHE.use = EstimParam.MHE.use;
estimator.MHE.Condensing = EstimParam.MHE.Condensing;

%% ---------- State Estimation  ----------

if estimator.use
    fprintf('\n------------------ State Estimator ------------------\n');
    
 % % %  observability checks  % % % 
    OB = obsv(model.plant.Ad,model.plant.Cd);
        if rank(OB) == model.plant.nx
            fprintf('*** simulation model is observable \n')
        else    
            fprintf('*** simulation model is not observable \n')
        end

    OB = obsv(model.pred.Ad,model.pred.Cd);
        if rank(OB) == model.pred.nx
            fprintf('*** prediction model is observable \n')
        else    
            fprintf('*** prediction model is not observable \n')
        end
        
 % % % Estimator design  % % %
%     if estimator.LOPP.use   % pole placement -- REMOVE????
%         fprintf('*** Create Luenberger Observer via pole placement ... \n')
%          
%         dMat = zeros(size(pred_mod.c,1),size(pred_mod.B,2));
%         MPCObj = mpc(ss(pred_mod.a,pred_mod.b,pred_mod.c,dMat,model.plant.Ts));
%         [~,M,A1,Cm1] = getEstimator(MPCObj);
%         e = eig(pred_mod.a-pred_mod.a*M(1:model.pred.nx,:)*pred_mod.c);
%         new_poles = abs(e) .* 0.9;
% 
%         try % Try to place the pole
%             estimator.LOPP.L1 = place(pred_mod.a',pred_mod.c',new_poles)';
%             disp('*** offset free = no')
%         catch
%             Qe = 10*eye(size(model.pred.Ed_estim,2));
%             Re = 10*eye(model.pred.ny);   
%             Ge = model.pred.Ed_estim;
%             estimator.LOPP.L1 = dlqe(model.pred.Ad, Ge, model.pred.Cd, Qe, Re);
%             disp('*** No poles placement was possible for this model.')
%         end
% %        model.estim.M = model.pred.Ad \ model.estim.L1;
%         fprintf('*** Done.\n')
%         
   if  estimator.SKF.use   % stationary Kalman estimator design 
        fprintf('*** Create SKF estimator ... \n')    
        
        estimator.SKF.Qe = 10*eye(model.pred.nx);     % process noise covariance 
        estimator.SKF.Re = 1*eye(model.pred.ny);      % measurement noise covariance 
        estimator.SKF.Ge = 1*eye(model.pred.nx);
        estimator.SKF.L1 = dlqe(model.pred.Ad, estimator.SKF.Ge, model.pred.Cd, estimator.SKF.Qe, estimator.SKF.Re);
        
        fprintf('*** Done.\n')  
        
    elseif estimator.TVKF.use  % time varying Kalman estimator design 
      fprintf('*** Create TVKF estimator ... \n')    
      
%          estimator.TVKF.Qe = 1e6;        % process noise covariance - works good for Reno
%          estimator.TVKF.Qe = 1e15;        % process noise covariance - good output estimation error crazy states for Infrax
         estimator.TVKF.Qe = 1e10;        % process noise covariance - relatively ok for Infrax
         estimator.TVKF.Re = 1*eye(model.pred.ny); 
         
      fprintf('*** Done.\n')   
      
    elseif estimator.MHE.use  
      fprintf('*** Create MHE estimator ... \n')  
      % horizons  
      estimator.MHE.N = 3; %  estimation horizon     
      % weight diagonal matrices 
%       Qe = 1e10;                                                       % process noise covariance magnitude
      Qe = 1e6;                                                       % process noise covariance magnitude
      Re = 1e0;                                                      % measurement noise covariance magnitude
      P = Qe;
      estimator.MHE.Qe = Qe*eye(model.pred.nx);                            % process noise covariance 
      estimator.MHE.Re = Re*eye(model.pred.ny);                             % measurement noise covariance 
      estimator.MHE.P = model.pred.Bd*P*model.pred.Bd';                   % error covariance of the arrival cost  

         %  construct MHE optimizer 
      estimator.MHE.optimizer = BeMHEdesign(model, estimator.MHE);
     
      fprintf('*** Done.\n')
      
    end
    
    
%     ------------------------------------------------
    if model.pred.off_free           
            disp('*** offset free = yes')
    else      
            disp('*** offset free = no')
    end
    
    %     Kalman estimator design for the system
    %        x[n+1] = Ax[n] + Bu[n] + Gw[n]    {State equation}
    %        y[n]   = Cx[n] + Du[n] +  v[n]    {Measurements}

    if estimator.SKF.use
            % check the stability of the estimator matrix LS
        if estimator.SKF.use
           L1_stable = eig(model.pred.Ad-model.pred.Ad*estimator.SKF.L1*model.pred.Cd);
        elseif estimator.LOPP.use
           L1_stable = eig(model.pred.Ad-model.pred.Ad*estimator.LOPP.L1*model.pred.Cd); 
        end
        
        if sum(L1_stable > 1) > 0
            fprintf('*** estimator is unstable \n')
        else
            fprintf('*** estimator is stable \n')
        end
    end
    
end
% A linear discrete-time system described by the state equation x(k + 1) = A x(k) + B u(k) is
% asymptotically stable if and only if all eigenvalues have magnitude smaller than one.



end