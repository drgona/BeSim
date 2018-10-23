function model = BuiModel(buildingType, ModelOrders, reload)
%% Description


%% Initiation
% default offset free control setup  
if nargin == 0  
   buildingType = 'Infrax'; 
end
if nargin < 2
   ModelOrders.range = [100, 200, 600]; % add any reduced order model you wish to have
   ModelOrders.choice = 200;            % insert model order or 'full' for full order SSM 
%   alternative choice of the model order - adopt to be general, possibly
% TODO:  abandon this feature and adopt residential model
   ModelOrders.ctrlModIndex = 9;
   ModelOrders.plantModIndex = 9;
   
   ModelOrders.off_free = 0;    %  augmented model
end
if nargin < 3
   reload = 0;    % reload SSMs and regenerate ROMs flag
end

% building parameters
path = ['../buildings/', buildingType];
disturbanceType = ''; % can be '_lin' if used for linearization validation
model.buildingType = buildingType;
model.Orders.range = ModelOrders.range;
model.Orders.choice =  ModelOrders.choice;
model.reload = reload;

fprintf('\n------------------ Building Model -------------------\n');

	%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Load model 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if reload       
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
        orders = ModelOrders.range;
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
if ModelOrders.choice == 'full' 
   ctrlModIndex = plantModIndex; 
else
   ctrlModIndex = find([orders, size(sys_dExt.a,1)] == ModelOrders.choice);  
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

    
if  ModelOrders.off_free % use augmented prediction model for offset free control
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
    model.pred.off_free = ModelOrders.off_free;   
    
%% Control input constraints


    model.pred.umax = 10000*ones(model.pred.nu,1);
%     model.pred.umax(1:4) = 10;
%     model.pred.umax(25:28) = 0;
%     model.pred.umax([17,16,8,5,11,12]) = 0;
    model.pred.umin = -20000*ones(model.pred.nu,1);
%     model.pred.umax(1:4) = -5;
%     model.pred.uin(25:28) = 0;
%     model.pred.umin([17,16,8,5,11,12]) = 0;



end