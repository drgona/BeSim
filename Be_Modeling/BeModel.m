function model = BeModel(buildingType, ModelParam)
%% Description


%% Initiation
% default offset free control setup  
if nargin == 0  
   buildingType = 'Infrax'; 
end
if nargin < 2
   ModelParam.Orders.range = [100, 200, 600]; % add any reduced order model you wish to have
   ModelParam.Orders.choice = 200;            % insert model order or 'full' for full order SSM 
%   alternative choice of the model order - adopt to be general, possibly
% TODO:  abandon this feature and adopt residential model
   ModelParam.Orders.ctrlModIndex = 9;
   ModelParam.Orders.plantModIndex = 9;
   
   ModelParam.Orders.off_free = 0;    %  augmented model
   ModelParam.reload = 0;    % reload SSMs and regenerate ROMs flag
end

% building parameters
path = ['../buildings/', buildingType];
disturbanceType = ''; % can be '_lin' if used for linearization validation
model.buildingType = buildingType;
model.Orders.range = ModelParam.Orders.range;
model.Orders.choice =  ModelParam.Orders.choice;
model.reload = ModelParam.reload;

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
            [t,  v, inputIndex, dictCtlInputs, dicValVar, dicOutputNameIndex, x0]= disturbances(path,disturbanceType,0, 0);
        end
        Ts = t(2) - t(1);
        orders = ModelParam.Orders.range;
%        states are initalized to x0 = 293.15 K in original model, extended
%        model initalizes states to 0 which is equivalent with  293.15 K  via matrix extension  
		[sys_dExt, rom] = fGenerateSysAndRom([path '/models/ssm.mat'], Ts, x0, orders);      
		save([path '/preComputed_matlab/mod.mat'], 'Ts', 'orders', 'sys_dExt', 'rom');
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

% u_index = 178:205;
% d_index = [1:177,206:287];
% ym_index = 1:19;


%% Linear SS Model
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

    
%% Control input constraints
    model.pred.umax = 10000*ones(model.pred.nu,1);
    model.pred.umin = -10000*ones(model.pred.nu,1);

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

if  ModelParam.Orders.off_free % use augmented prediction model for offset free control
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
    model.pred.off_free = ModelParam.Orders.off_free;   
    


end