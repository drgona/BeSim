function model = BeModel(buildingType, ModelParam)
%% Description


%% Initiation
% default offset free control setup  
if nargin == 0  
   buildingType = 'Reno'; 
end
if nargin < 2
    ModelParam.Orders.range = [7, 15, 30, 100];    % vector of model orders 
%     ModelParam.Orders.range = [4, 7, 10, 15, 20, 30, 40, 100];    % vector of model orders 
    ModelParam.Orders.choice = 'full';                            % model order selection for prediction
%    ModelParam.Orders.range = [100, 200, 600]; % add any reduced order model you wish to have
%    ModelParam.Orders.choice = 200;            % insert model order or 'full' for full order SSM 
%   alternative choice of the model order - adopt to be general, possibly
% TODO:  abandon this feature and adopt residential model
    ModelParam.Orders.ctrlModIndex = 9;
    ModelParam.Orders.plantModIndex = 9; 
    ModelParam.off_free = 0;    %  augmented model
    ModelParam.reload = 0;    % reload SSMs and regenerate ROMs flag
    ModelParam.analyze.SimSteps = 2*672; % Number of simulation steps (Ts = 900 s),  672 = one week
    ModelParam.analyze.openLoop.use = true;             %  open loop simulation   - TODO
    ModelParam.analyze.openLoop.start = 1;              % starting day of the analysis
    ModelParam.analyze.openLoop.end = 7;                % ending day of the analysis
    ModelParam.analyze.nStepAhead.use = true;           % n-step ahead predicion error  - TODO
    ModelParam.analyze.nStepAhead.steps = [1, 10, 40];  % x*Ts  
    ModelParam.analyze.HSV = false;                      %  hankel singular values of ROM
    ModelParam.analyze.frequency = false;                % frequency analysis - TODO

end

% building parameters
path = ['../buildings/', buildingType];
disturbanceType = ''; % can be '_lin' if used for linearization validation
model.buildingType = buildingType;
model.Orders.range = ModelParam.Orders.range;
model.Orders.choice =  ModelParam.Orders.choice;
model.reload = ModelParam.reload;
% model analysis
model.analyze.SimSteps = ModelParam.analyze.SimSteps;
model.analyze.openLoop.use = ModelParam.analyze.openLoop.use;
model.analyze.openLoop.start = ModelParam.analyze.openLoop.start;
model.analyze.openLoop.end = ModelParam.analyze.openLoop.end;
model.analyze.nStepAhead.use = ModelParam.analyze.nStepAhead.use;
model.analyze.nStepAhead.steps = ModelParam.analyze.nStepAhead.steps;
model.analyze.HSV = ModelParam.analyze.HSV;
model.analyze.frequency = ModelParam.analyze.frequency;

fprintf('\n------------------ Building Model -------------------\n');

	%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Load model 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if ModelParam.reload  
		% Model
		fprintf('*** Create ROM models ...\n')
        % Disturbances
        
%         TODO: get rid of this if else by unifying disturbances for house model 
        if strcmp(buildingType,'Reno') || strcmp(buildingType,'RenoLight') || strcmp(buildingType,'Old')  
            [t, v, x0] = disturbances_old(path, 0, 0);
        else
%             [t,  v, inputIndex, dictCtlInputs, dicValVar, dicOutputNameIndex, x0]= disturbances(path,disturbanceType,0, 0);
            [t,  v, ~, ~, ~, ~, x0]= disturbances(path,disturbanceType,0, 0);
        end
        Ts = t(2) - t(1);
        orders = ModelParam.Orders.range;
%        states are initalized to x0 = 293.15 K in original model, extended
%        model initalizes states to 0 which is equivalent with  293.15 K  via matrix extension  
		[sys_dExt, rom, redinfo] = fGenerateSysAndRom([path '/models/ssm.mat'], Ts, x0, orders);      
		save([path '/preComputed_matlab/mod.mat'], 'Ts', 'orders', 'sys_dExt', 'rom', 'v', 't', 'x0');
   		fprintf('*** Done.\n')
    else
        fprintf('*** Load ROM models ...\n')
		load([path '/preComputed_matlab/mod.mat']);
        fprintf('*** Done.\n')
%         TODO:  create mod.mat file also for 6-zone building - connect
%         models in one file
    end


%     load indexing of the models: inputs u_index, disturbances d_index,
%     measured outputs ym_index
load([path '/preComputed_matlab/indexing.mat']);


%% Linear SS Model
fprintf('*** Simulation and prediction model construction ...\n')
% ---------- Model orders selection ----------
NM = size(rom,1);       % number of investigated reduced order models
plantModIndex = NM+1;   % plant model index
% controller ROM index 
if ModelParam.Orders.choice == 'full' 
   ctrlModIndex = plantModIndex; 
else
   ctrlModIndex = find([orders, size(sys_dExt.a,1)] == ModelParam.Orders.choice);  
end

% ---------- Plant model  ----------      
    % % plant model structure
    % x_k+1 = Ad*x_k + Bd*u_k + Ed*d_0 + Gd*1 
    % yk = Cd*x_k + Dd*u_k + Fd*1

plant_mod = sys_dExt;       %  full SSM

model.plant.Ts = plant_mod.Ts;  % simulation sampling time

% construction of simulation model
model.plant.Ad = plant_mod.A; 
% % separation of the disturbances matrix E from control inputs matrix B
model.plant.Bd = plant_mod.B(:,u_index); 
model.plant.Ed = plant_mod.B(:,d_index); 
% reduction of the C and D matrix based on available measurements
model.plant.Cd = plant_mod.C(ym_index,:); 
model.plant.Dd = plant_mod.D(ym_index,u_index);
%  model initialization extension matrices
model.plant.Gd = plant_mod.b(:,end);             % extension of original matrix Bd - disturbances compensation matrix
model.plant.Fd = plant_mod.d(ym_index,end);      % extension of original matrix Dd - kelvins to celsius compensation matrix
% Overall sim. model dimensions
model.plant.nx = size(model.plant.Ad, 2);
model.plant.ny = size(model.plant.Cd, 1);
model.plant.nd = size(model.plant.Ed, 2);
model.plant.nu = size(model.plant.Bd, 2);


% ---------- Choice of the Controller model ----------
if ctrlModIndex < plantModIndex
    pred_mod = rom{ctrlModIndex};   % use reduced order model
elseif ctrlModIndex == plantModIndex
    pred_mod = sys_dExt;                % use original model
end

% ---------- Controller/Prediction model matrices and dimensions ----------
model.pred.Ts = pred_mod.Ts;  % simulation sampling time

    % construction of prediction model
    model.pred.Ad = pred_mod.A; 
    % % separation of the disturbances matrix E from control inputs matrix B
    model.pred.Bd = pred_mod.B(:,u_index); 
    model.pred.Ed = pred_mod.B(:,d_index); 
    % reduction of the C and D matrix based on available measurements
    model.pred.Cd = pred_mod.C(ym_index,:); 
    model.pred.Dd = pred_mod.D(ym_index,u_index);
    %  model initialization extension matrices
    model.pred.Gd = pred_mod.b(:,end);     % extension of original matrix Bd
    model.pred.Fd = pred_mod.d(ym_index,end);      % extension of original matrix Dd
    % Overall sim. model dimensions
    model.pred.nx = size(model.pred.Ad, 2);
    model.pred.ny = size(model.pred.Cd, 1);
    model.pred.nd = size(model.pred.Ed, 2);
    model.pred.nu = size(model.pred.Bd, 2);

    
%% control input constraints
    model.pred.umax = 1000*ones(model.pred.nu,1);
    model.pred.umin = -1000*ones(model.pred.nu,1);

%     TODO: generalize this based on loading constraint matrix
if  strcmp(buildingType,'Reno')
    model.pred.umax =[1680 685 154 1000 320 232]';
    model.pred.umin = zeros(model.pred.nu,1);
elseif strcmp(buildingType,'Old')
    model.pred.umax = [2940 960 300 1400 460 253]';
    model.pred.umin = zeros(model.pred.nu,1);  
elseif strcmp(buildingType,'RenoLight')
    model.pred.umax = [840 343 77 500 160 116]';
    model.pred.umin = zeros(model.pred.nu,1);   
end
    
%% post processing of individual models
if  strcmp(buildingType,'HollandschHuys')
    
    % supply vent temp as disturbance
    Tven_index = 17:28;
    model.plant.Ed = plant_mod.B(:,[d_index, Tven_index]);  %  adding supply temp for ventilation
    model.plant.nd = size(model.plant.Ed, 2);
    model.pred.Ed = pred_mod.B(:,[d_index, Tven_index]);  %  adding supply temp for ventilation
    model.pred.nd = size(model.pred.Ed, 2);
     
    % min-max heat flows circuits
    % Q = valve*m_nominal*1.159*dT;
    dT = 4;   % TODO: specify this - make this more accurate assumption
    model.m_nominal = [1283 1188 1393 1198 2471 3773 2619 2360 1885 1829 1350 1638 3952 1583 2786 2237 1433 1874 1792 1782];
    Qmax = model.m_nominal*1.159*dT;
    Qmin = -model.m_nominal*1.159*dT;

    model.pred.umax = Qmax';
    model.pred.umin = Qmin';

%     % total heat limits
%     model.pred.umax = 181000;
%     model.pred.umin = -90000;

    % model pre-processing
    % lumping interconnected inputs - see excel document
    oldB = model.plant.Bd;
    newB = [oldB(:,1:4), 0.5*oldB(:,5) + 0.5*oldB(:,12),0.5*oldB(:,6) + 0.5*oldB(:,13), ...
       0.5*oldB(:,7) + 0.5*oldB(:,15), oldB(:,8:11),  oldB(:,14), 0.5*oldB(:,16) + 0.5*oldB(:,28), ...
      0.5*oldB(:,17) + 0.5*oldB(:,29), oldB(:,18:19), 0.5*oldB(:,20) + 0.5*oldB(:,24), ...
      0.5*oldB(:,21) + 0.5*oldB(:,25), 0.5*oldB(:,22) + 0.5*oldB(:,26), 0.5*oldB(:,23) + 0.5*oldB(:,27)];
    % newB = [oldB(:,1:4), sum(oldB(:,[5,12]),2), sum(oldB(:,[6,13]),2), sum(oldB(:,[7,15]),2), ...
    %     oldB(:,8:11),  oldB(:,14), sum(oldB(:,[16,28]),2), sum(oldB(:,[17,29]),2), oldB(:,18:19), ...
    %     sum(oldB(:,[20,24]),2), sum(oldB(:,[21,25]),2), sum(oldB(:,[22,26]),2), sum(oldB(:,[23,27]),2)];
    model.plant.Bd = newB;
    model.plant.nu = size(model.plant.Bd,2);
    model.plant.Dd = model.plant.Dd(:,1:model.plant.nu);

    oldB = model.pred.Bd;
    % newB = [oldB(:,1:4), sum(oldB(:,[5,12]),2), sum(oldB(:,[6,13]),2), sum(oldB(:,[7,15]),2), ...
    %     oldB(:,8:11),  oldB(:,14), sum(oldB(:,[16,28]),2), sum(oldB(:,[17,29]),2), oldB(:,18:19), ...
    %     sum(oldB(:,[20,24]),2), sum(oldB(:,[21,25]),2), sum(oldB(:,[22,26]),2), sum(oldB(:,[23,27]),2)];
    newB = [oldB(:,1:4), 0.5*oldB(:,5) + 0.5*oldB(:,12),0.5*oldB(:,6) + 0.5*oldB(:,13), ...
       0.5*oldB(:,7) + 0.5*oldB(:,15), oldB(:,8:11),  oldB(:,14), 0.5*oldB(:,16) + 0.5*oldB(:,28), ...
      0.5*oldB(:,17) + 0.5*oldB(:,29), oldB(:,18:19), 0.5*oldB(:,20) + 0.5*oldB(:,24), ...
      0.5*oldB(:,21) + 0.5*oldB(:,25), 0.5*oldB(:,22) + 0.5*oldB(:,26), 0.5*oldB(:,23) + 0.5*oldB(:,27)];
    model.pred.Bd = newB;
    model.pred.nu = size(model.pred.Bd,2);
    model.pred.Dd = model.pred.Dd(:,1:model.pred.nu);
    
%     lump constraints
    model.pred.umax = model.pred.umax(1:model.pred.nu);
    model.pred.umin = model.pred.umin(1:model.pred.nu);  

end

if  ModelParam.off_free % use augmented prediction model for offset free control
    % number of output disturbances p_k =  number of outputs
    model.pred.np = size(model.pred.Cd, 1);
    % output disturbacne matrix- design conditions:  mag(Gp) < mag(Cd)
    model.pred.Gp = 0.1*eye(model.pred.ny,model.pred.np);
    % % augmented model with output disturbances p_k
    % xS_k+1 = AS*xS_k + BS*uk + ES*d_0 + GS*1 
    % yk = CS*xS_k + Fd*1 + Gp+p
    % xS_k = [x_k, p_k]
    model.pred.Ad = [model.pred.Ad, zeros(model.pred.nx,model.pred.np) ; zeros(model.pred.np,model.pred.nx), eye(model.pred.np)];
    model.pred.Bd = [model.pred.Bd; zeros(model.pred.np,model.pred.nu)];
    model.pred.Ed = [model.pred.Ed; zeros(model.pred.np,model.pred.nd)];
    model.pred.Gd = [model.pred.Gd; zeros(model.pred.np,1)];
    model.pred.Cd = [model.pred.Cd, model.pred.Gp];
    % Overall estim. model dimensions
    model.pred.nx = size(model.pred.Ad, 2);
    model.pred.ny = size(model.pred.Cd, 1);
    model.pred.nd = size(model.pred.Ed, 2);
    model.pred.nu = size(model.pred.Bd, 2); 
end
    %  offset free control indicator
    model.pred.off_free = ModelParam.off_free;   
    
%% Model analysis
%     analysis = BeModelAccuracy(model);

% TODO: implement functionality and plotting from function
% offLinePredPerf.m

% TODO: generate control actions U from simulation for every building
if (model.analyze.openLoop.use || model.analyze.nStepAhead.use || model.analyze.HSV ||  model.analyze.frequency)
    fprintf('*** Model analysis ...\n')
    
%     TODO: analysis error for Hol Huys model
    
        % % % % % % % %  open loop  analysis % % % % % % % %    
    if (model.analyze.openLoop.use || model.analyze.nStepAhead.use)   
        % Generate control actions by weekly oscilations between umin and umax
        sine_waves = time_transform(size(v,1),model.plant.Ts);
        U = zeros(size(v,1),model.plant.nu);  
        U(sine_waves(:,2)>=0,:) = sine_waves(sine_waves(:,2)>=0,2)*model.pred.umax'; 
        U(sine_waves(:,2)<0,:) = -sine_waves(sine_waves(:,2)<0,2)*model.pred.umin'; 
        % disturbance profiles
        D = v; 
        % lumped model inputs
        UExt = [U, D, ones(size(v,1),1)]';
        t_sim = (1:size(v,1))*sys_dExt.Ts;
        HSim = model.analyze.SimSteps;   
        
        %         Turn off lsim warnings
        id = 'Control:analysis:LsimStartTime';
        warning('off',id)
        
        % simulate the different models. 
        ny = size(sys_dExt.C, 1);
        yComp = zeros( length(orders)+1, HSim, ny );
            % Full SSM
        [yComp(end, :, :),t,xSSM] = lsim(sys_dExt,UExt(:,1:HSim),t_sim(1:HSim));
        leg = cell(length(orders)+1,1);
            % Reduced order models
        xROM = {};
        for i = 1:length(orders)
            [yComp(i, :, :),t,xROM{i}] = lsim(rom{i},UExt(:,1:HSim),t_sim(1:HSim));
            leg{i} = num2str(orders(i));
        end
        leg{end} = 'SSM';
        yComp(:, :, :) = yComp(:, :, :) -273.15;      
    end
    
     if model.analyze.openLoop.use       
        sampling = 14;
        markers = {'d','o','^','v','x','s','+','<','>','*','--','-*'};        
        % Plot simulation result for each zone
        fig = figure;
        time = t/3600/24; %in days  
        columns = ceil(ny/2);
        rows = 2;
        ha = tight_subplot(rows,columns,[.1 .03],[.15 .01],[.08 .01]);
        for i = 1:ny
            axes(ha(i));
            hold off
            for k = 1:length(orders)
                p = plot(time(1:sampling:end) , yComp(k, 1:sampling:end, i)); hold on;
                p.Marker = markers{k};
            end
            p = plot(time , yComp(end, :, i));
            p.LineWidth = 2;
            p.Color = 'k';
            p.LineStyle = '--';
%             ylim([15,25]);
%             xlim([0,7]);
            % Set title into plot
            titlePlot = title(['TZone ' num2str(i)]);
            pos = get(titlePlot,'Position');
            set(titlePlot,'Position',[pos(1) pos(2)-1.2 pos(3)])
            grid on
            if i==1 || i==4
                ylabel('T [^\circ C]')
            elseif i==5
                xlabel('Time [Day]')
                legend(leg, 'Orientation', 'horizontal')
            end
        end
%         % remove xLabel from first row
%         set(ha(1:3),'XTickLabel',''); 
%         % remove yLabel from second column 
%         set(ha(2:3),'YTickLabel','')
%         set(ha(5:6),'YTickLabel','')

        % Plot simulation result for zone 1 for all SSM and ROM
        % The black line is the reference, the dashed line is the full
        % order SSM
        % FIXME: legend symbole for SSM is not correct. Why?
        figure
        for i = 1:length(orders)
            p = plot(time(1:sampling:end) , yComp(i, 1:sampling:end, 1)); hold on;
            p.Marker = markers{i};
        end
        p = plot(time(1:sampling:end) , yComp(end, 1:sampling:end, 1));
        p.LineWidth = 2;
        p(end).Color = 'k';
        p(end).LineStyle = '--';
        legend(leg);
        ylabel('T [^\circ C]');
        xlabel('Day');
        title(['TZone ' num2str(1)]);
        grid on
%         plot(time,Rmin(:,1:HSim),'k','LineWidth',0.5);
        %plot(time,Rmax(:,1:opt.HSim),'k','LineWidth',0.5);

        % Plot zones average simulation result for all SSM and ROM
        % The black line is the reference, the dashed line is the full
        % order SSM
        % FIXME: legend symbole for SSM is not correct. Why?
        figure      
        for i = 1:length(orders)
            yAve = mean(yComp(i, 1:sampling:end,:),3);
            p = plot(time(1:sampling:end) ,yAve); hold on;
            p.Marker = markers{i};
        end
        p = plot(time(1:sampling:end) , yComp(end, 1:sampling:end, 1));
        p.LineWidth = 2;
        p(end).Color = 'k';
        p(end).LineStyle = '--';
        legend(leg);
        ylabel('T [^\circ C]');
        xlabel('Day');
        title(['Average(TZones)']);
        grid on
%         plot(time,Rmin(:,1:opt.HSim),'k','LineWidth',0.5);
     end
    
    
      if model.analyze.nStepAhead.use
        % generate and plot for n-steps ahead predictions 
        
%         initial parameters
        nStepAhead = model.analyze.nStepAhead.steps;
        H = length(nStepAhead);
        start = max(nStepAhead)+2;    
        yNStepAhead = zeros( length(orders) + 1, H, HSim, ny); % [order, time, zone]
        
        Mods = rom;
        Mods{length(rom)+1} = sys_dExt;
        
        for i = 1:(length(orders)+1)
                                
%       state initialization x0 of ROMs based on Kalman Filter
            sysMod = Mods{i};  
            nx = size(sysMod.A,1);
            ny = size(sysMod.C,1);
            Qe = 1e10;     % process noise covariance 
            Re = 1*eye(ny);      % measurement noise covariance    
            Xp = zeros(nx,HSim+1);
            Xe = zeros(nx,HSim);
            Ye = zeros(ny,HSim);
            Yp = zeros(ny,HSim);
            xp = Xp(:, 1); 
%             KF simulation
            for k = 1:HSim
                  if k == 1
                      P = sysMod.B*Qe*sysMod.B';         % Initial error covariance   
                  end
                  
                  uopt = UExt(:,k);
                  yn = squeeze(yComp(i, k, :))+273.15;
                  
                  % Measurement update
                  L1 = P*sysMod.C'/(sysMod.C*P*sysMod.C'+Re); % observer gain
                  yp = sysMod.C*xp + sysMod.D*uopt;          % output estimation
                  ep = yn - yp;                                                       % estimation error
                  xe = xp + L1*ep;                                                    % x[n|n]
                  P = (eye(nx)-L1*sysMod.C)*P;                          % P[n|n]   estimation error covariance
                  % Time update
                  xp = sysMod.A*xe + sysMod.B*uopt;        % x[n+1|n]
                  P = sysMod.A*P*sysMod.A' + sysMod.B*Qe*sysMod.B';       % P[n+1|n]          
                  ye = sysMod.C*xe + sysMod.D*uopt;     % output estimate with x[n|n]

                  Xp(:,k+1) = xp;
                  Xe(:,k) = xe;
                  Ye(:,k) = ye;
                  Yp(:,k) = yp;
            end
%             figure
%             plot(squeeze(yComp(i, :, :))+273.15-Ye')
%             figure
%             plot(Xe')
            
            
            for z=1:H   % --------- Different h-step ahead horizons
                h = nStepAhead(z);
                for t_end = start:HSim  % --------- Compute the h-step ahead prediction for all point from start to HSim
                    t_start = t_end - h;
            
                    if i == (length(orders)+1)
                        x_init = xSSM(t_start,:)';
                    else
            %              x_init =  rom{i}.C\squeeze(yComp(end,t_start,:));
                        x_init = Xe(:,t_start);
                    end
                    yTemp = lsim(Mods{i}, UExt(:,t_start:t_end),t_sim(t_start:t_end), x_init);
                    yNStepAhead(i, z, t_end, :) = yTemp(end,:)-273.15;
                    
                end      
            end      
        end
       
                 
% %         if opt.computeNSteps
%             for z=1:H   % --------- Different h-step ahead horizons
%                 h = nStepAhead(z);
%                 check = 0;
%                 for t_end = start:HSim  % --------- Compute the h-step ahead prediction for all point from start to HSim
%                     t_start = t_end - h;
%                     % Re-create ROM with correct initial values 
% %                     [sys_dExtXN, romXN] = fGenerateSysAndRom([path '/models/ssm.mat'], Ts, xSSM(t_start,:)' + x0, orders);
%                     for i = 1:length(orders) % --------- Compute prediction for each ROM
%                         % Simulate from t = t_end - h to t_end to get h-step ahead prediction      
%                         
% %                         TODO: project initial conditions of SSM onto ROM
% %                         solving system of lin equations
% %                         TODO: improve initialization with kalman filter  
%                         xROM_init =  rom{i}.C\squeeze(yComp(end,t_start,:));
%                         yTemp = lsim(rom{i}, UExt(:,t_start:t_end),t_sim(t_start:t_end), xROM_init);
% %                         yTemp = lsim(rom{i}, UExt(:,t_start:t_end),t_sim(t_start:t_end), xROM{i}(t_start,:)');
%                         % Save h-step ahead prediction
%                         
% %                         yTemp = lsim(romXN{i}, UExt(:,t_start:t_end),t_sim(t_start:t_end));                      
%                         yNStepAhead(i, z, t_end, :) = yTemp(end,:)-273.15;
%                     end
%                     yTemp = lsim(sys_dExt, UExt(:,t_start:t_end),t_sim(t_start:t_end), xSSM(t_start,:)');
% %                     yTemp = lsim( sys_dExtXN, UExt(:,t_start:t_end),t_sim(t_start:t_end));
%                     % Save h-step ahead prediction
%                     yNStepAhead(end, z, t_end, :) = yTemp(end,:)-273.15;
%                     check = max( check, abs(squeeze(yTemp(end,:)-273.15) - squeeze(yComp(end,t_end,:))') );
%                 end
%                 disp(['Check = ' num2str(check) ' for h = ' num2str(h) ' of building ' buildingType])
%             end    
            
            
        % Take temperature different between yNStepAhead and reference and
        % reshape yNSteAhead by appending all zones to 1 column.
        % Also compute 1-NRMSE
        yRef = squeeze(yComp(end,start:end,:));
        HSim_new = HSim - start + 1;
        deltaTNStepAhead = zeros( length(orders), H, HSim_new*ny);
        NRMSE_NStepAhead = zeros( length(orders), H, ny);
        for z=1:H
            for i=1:length(orders)
                yNStepAhead_iz = squeeze(yNStepAhead(i, z, start:end, :));
                deltaTNStepAhead(i, z, :) = reshape( yNStepAhead_iz - yRef, [], 1); %squeeze removes singleton dimensions
                ave = mean( yNStepAhead_iz );
                NRMSE_NStepAhead(i, z, :) = ones(1,ny) - sqrt( sum( (yNStepAhead_iz - yRef).^2 ) ) ./ ...
                    sqrt( sum( (yNStepAhead_iz - repmat(ave, HSim_new,1)).^2 ) );
            end
        end
        
        fig = figure;
        ylabel('[K]');
        leg = cell(length(orders),1);
        for i=1:length(orders)
            leg{i} = ['ROM ' num2str(orders(i))];
        end
        for z=1:H
            p(z) = subplot(1,H,z);
            deltaTNStepAhead_z = squeeze(deltaTNStepAhead(:,z,:));
            boxplot(deltaTNStepAhead_z', 'labels', leg, 'labelorientation','inline')
            title([num2str( nStepAhead(z) ) '-steps ahead'])
%             ylim([-3, 2])
        end
        % Add common label by adding invisible axes to subplot 1
        h=axes('position',[p(1).Position],'visible','off');
        ylabel('[K]','visible','on');
  
%         if opt.saveNSteps
%             % resize page size. This set the boundingBox of the eps file
%             figuresize(fig, 407, 279, 'points');
%             % Save figure as eps. Argument epsc is for color printing
%             saveas(fig,['PredPerf' buildingType{j} '.eps'],'epsc')
%         end
      end
  
    
    % % % % % % % %  HSV  analysis % % % % % % % %    
    % In state coordinates that equalize the input-to-state and state-to-output energy transfers
    % Hankel singular values  measure the contribution of each state to the input/output behavior.
    % Hankel singular values are to model order what singular values are  to matrix rank.
    if model.analyze.HSV           
%         TODO: HSV values lines paper form
        
         figure % HSV for full order model
        [hsvd_val_ssm, hsv_data_ssm] = hsvd(sys_dExt);
        bar(hsvd_val_ssm/max(hsvd_val_ssm),0.8)
          title(['HSV of full order model: ' buildingType],'FontSize',16)
%            title(['Original Model HSV (state contribution)'],'FontSize',12)          
%          legend hide
        xlabel('State []','FontSize',14)
         ylabel('State Energy []','FontSize',14)
         set(gca,'YScale', 'log','FontSize',14)
         grid on
         axis tight 
         
         figure % HSV for ROM models
         nCol = ceil( ( length(model.Orders.range )) /2  );
         for i=1:length(model.Orders.range )
            subplot(2, nCol ,i)
           [hsvd_val{i} hsv_data{i}] =  hsvd(rom{i});
%            hsvd(rom{i})
           HSV_rom{i} = hsvd_val{i}/max(hsvd_val{i});
           bar(HSV_rom{i},0.8)
%             title([buildingType(j) 'Model Order ' num2str(orders(i)) ]);
            title(['HSV of ROM: ' num2str(model.Orders.range(i)) ],'FontSize',12);
%             legend hide
            set(gca,'YScale', 'log','FontSize',9)
            grid on
            axis tight            
            if i ~= 1 && i ~= length(model.Orders.range)/2+1
                     ylabel('')
            else
               ylabel('State Energy','FontSize',12) 
            end
            if i > length(model.Orders.range )/2
                xlabel('State','FontSize',12)               
            else
                xlabel('')
            end                                         
         end
    end  

    % % % % % % % %  Frequency  analysis % % % % % % % %    
    if model.analyze.frequency
        % Option for bode plots
        opts = bodeoptions('cstprefs');
        opts.PhaseVisible = 'off';
        opts.MagLowerLimMode = 'manual';
        opts.MagLowerLim = -150;
        opts.Ylim = [-150, 0];
        opts.Xlim = [10^(-9), 10^(0) ];
        opts.FreqUnits = 'Hz';
        opts.Grid = 'on';
        opts.Title.String = '';
        
%         TODO: automate this to number of ROMs and all combinations of I/O
        figure
        h = bodeplot(sys_dExt,'k', opts, rom{1},'g--', opts, rom{2},'b:', opts);        
    end
    
    
% %     NOTES: initialization of ROM x0 via HSV projections of full SSM - unfinished 
    %         K = 10;
%         SimInputs = UExt(:,1:HSim)';
%         for i = 1:length(orders) 
%             SimOutputs{i} = squeeze(yComp(i,:,:))
%             SimData{i} = iddata(SimOutputs{i},SimInputs,Ts);
%             SimData{i} = [SimOutputs{i}, SimInputs];
% %             err = pe(rom{i},SimData{i},K);
% %             compare(SimData{i},rom{i},K);
% %             yp = predict(rom{i},SimData{i},K);
%         end           
        
% 
% % % [HSV_STAB,HSV_UNSTAB] = BeHankelsv(rom{2})
% [HSV_STAB,HSV_UNSTAB] = hankelsv(sys_dExt)
% 
% G = sys_dExt;
% 
% Ts = abs(get(G,'Ts'));
% 
% % Handle discrete case:
% 
% paaz = 2^-4*sqrt(eps);
% 
% if Ts
%     a = get(G,'a');
%     if (max(abs(eig(a)+1))>paaz)  % if no poles near z=-1
%         paaz = 0;
%     end
%     G = bilin(G,-1,'S_Tust',[Ts, 1-paaz]);
% end
% 
% hanksvstab = [];
% hanksvunstab = [];
% Gorig = G;
% 
% function r = localAssignRegion(p,jwtol)
% if real(p)<-jwtol
%    r = 1;
% elseif real(p)>jwtol
%    r = 2;
% else
%    r = 3;
% end

% % Split spectrum into stable/unstable/jw-axis
% jwtol = sqrt(eps) * norm(balance(G.A),1);
% [H,H0] = modsep(G,3,@(p) AssignRegion(p,jwtol));
% nx = order(H);
% m = nx(1);  no_unstable = nx(2);  no_jw_pole = nx(3); % no. stable, unstable, jw-axis poles
% mjws = m+no_jw_pole;   ra = mjws+no_unstable;
% [a,b,c,d] = ssdata(H,'cell');
% G = H0 + ss(a{1},b{1},c{1},d{1}) + ss(a{2},b{2},c{2},d{2});
% % G = H0 + subpar(H,{':',':',1}) + subpar(H,{':',':',2});
% [a,b,c,d] = ssdata(G);
% 
% 
% a = sys_dExt.a;
% b = sys_dExt.b;
% c = sys_dExt.c;
% 
% %            [a,b,c,d] = ssdata(sys_dExt);
%            p = gram(a,b);
%            q = gram(a',c');
%            [up,sp,vp] = svd(p);
%            lp = up*diag(sqrt(diag(sp)));
%            [uq,sq,vq] = svd(q);
%            lo = uq*diag(sqrt(diag(sq)));
%            [U,Sigma,V]  = svd(lo'*lp);
%            
% %            https://nl.mathworks.com/help/robust/ref/balancmr.html
%            SLbig = lo*U(:,1:30)*diag(sqrt(diag(Sigma(1:30,1:30))));
%            SRbig = lp*V(:,1:30)*diag(sqrt(diag(Sigma(1:30,1:30))));
%            
%            AA =SLbig'*a*SRbig;
%            
% 
% [hsvd_val_ssm, hsv_data_ssm] = hsvd(sys_dExt);

end

end