function mhe = BeMHEdesign(model, MHEParam)

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
   MHEParam.N = 3; %  estimation horizon  
   % weights - covariances
   Qe = 1e6;                                     % process noise covariance magnitude
   Re = 1e0;                                     % measurement noise covariance magnitude
   MHEParam.Qe = Qe*eye(model.pred.nx);                            % process noise covariance 
   MHEParam.Re = Re*eye(model.pred.ny);                             % measurement noise covariance 
   MHEParam.P = model.pred.Bd*Qe*model.pred.Bd';                   % error covariance of the arrival cost  
   MHEParam.Condensing = 1;
end

 
%% MHE parameters
        
nx = model.pred.nx;
nu = model.pred.nu;
nd = model.pred.nd;
ny = model.pred.ny;
N =  MHEParam.N;

% variables
u = sdpvar(nu, N, 'full'); % ctrl actions 
d = sdpvar(nd, N, 'full'); % disturbances  
y = sdpvar(ny, N, 'full'); % outputs      
if not(MHEParam.Condensing)
    w = sdpvar(nx, N, 'full'); %  state update error slack
    x = sdpvar(nx, N+1, 'full'); % states of the building
else
    x = sdpvar(nx, 1, 'full'); % states of the building
end
v = sdpvar(ny, N, 'full'); %  output update error slack
s = sdpvar(nx, 1, 'full');           % arrival cost (x_{k-N+1})

%% MHE problem formulation
%  objective function+ constraints init
obj = 0;
con = [];
       
% state condensing matrices initialization
AB = zeros( nx , N*nu );
AE = zeros( nx , N*nd );
AG = zeros( nx, N*1);
AExpX0 = eye(nx) * x(:,1);
      
for k = 1:N   
 
        
        %   -------------  Constraints  -------------
        
       %     state + output update equations
        if MHEParam.Condensing 
            if k == 1
                AB(:, (N-k)*nu+1:(N-k+1)*nu ) = model.pred.Bd;        %  input matrix evolution
                AE(:, (N-k)*nd+1:(N-k+1)*nd ) =  model.pred.Ed;       %  disturbance matrix evolution
                AG(:, (N-k)*1+1:(N-k+1)*1 ) =  model.pred.Gd;         %  initial conditions matrix evolution
                con = con + [ y(:, k) == model.pred.Cd*x(:, k)  + model.pred.Dd*u(:, k) + model.pred.Fd*1  + v(:, k)];
            else
                AExpX0 = model.pred.Ad * AExpX0;
                con = con + [ y(:, k) == model.pred.Cd*( AExpX0 + AB(:, (N-k+1)*nu+1 : end ) * reshape( u(:,1:k-1) , nu * (k-1) , 1) + ...
                                                                 AE(:, (N-k+1)*nd+1 : end ) * reshape( d(:,1:k-1) , nd * (k-1) , 1) + ...
                                                                 AG(:, (N-k+1)*1+1 : end ) * ones(k-1,1) ) + ...
                                                                 model.pred.Dd*u(:, k)  + model.pred.Fd*1  + v(:, k)];

                AB(:, (N-k)*nu+1:(N-k+1)*nu ) = model.pred.Ad* AB(:, (N-k+1)*nu+1:(N-k+2)*nu );
                AE(:, (N-k)*nd+1:(N-k+1)*nd ) = model.pred.Ad* AE(:, (N-k+1)*nd+1:(N-k+2)*nd );
                AG(:, (N-k)*1+1:(N-k+1)*1 )   = model.pred.Ad* AG(:, (N-k+1)*1+1:(N-k+2)*1 );
            end    
        else                
            con = con + [ x(:, k+1) == model.pred.Ad*x(:, k) + model.pred.Bd*u(:, k) + model.pred.Ed*d(:, k) + model.pred.Gd*1 + w(:, k)];
            con = con + [ y(:, k) == model.pred.Cd*x(:, k)  + model.pred.Dd*u(:, k)  + model.pred.Fd*1 + v(:, k)];
        end
        
        % output estimation error penalization constraints        
        con = con + [  -5*ones(ny,1) <= v(:, k) <= 5*ones(ny,1)];
        
         %   state constraints   
        if not(MHEParam.Condensing)        
           con = con + [  -20*ones(nx,1) <= x(:, k) <= 30*ones(nx,1)];
           con = con + [  -30*ones(nx,1) <= w(:, k) <= 30*ones(nx,1)];   
        elseif k == 1
           con = con + [  -20*ones(nx,1) <= x(:, 1) <= 30*ones(nx,1)];
        end
        

        %   -------------  OBJECTIVE FUNCTION  -------------
         if MHEParam.Condensing 
             obj = obj +  v(:,k)'*MHEParam.Re*v(:,k);
         else
             obj = obj + w(:,k)'*MHEParam.Qe*w(:,k) + ...         %  state update error penalization
                        v(:,k)'*MHEParam.Re*v(:,k);               %  output update error penalization
         end
        
end
    
    obj = obj +  (x(:,1)-s)'*MHEParam.P*(x(:,1)-s);              %  arrival cost term to the first element in the horizon (x_{k-N+1})
    
     %% construction of object optimizer
     %   structure:  optimizer(constraints, objecttive, options, input_params, output_params)
 
    %  optimizer options
    % options = sdpsettings('verbose', 1, 'warning', 1, 'beeponproblem', 1, 'solver','cplex');
    options = sdpsettings('verbose', 1, 'solver','gurobi');

if MHEParam.Condensing 
   mhe = optimizer(con, obj, options,  {y, u, d, s}, {x(:,1); v; obj } );  
else
   mhe = optimizer(con, obj, options,  {y, u, d, s}, {x(:,1); v; obj; w } ); 
end


end