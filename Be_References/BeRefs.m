function references = BeRefs(model, dist, RefsParam)

% TODO:  FIX THIS 
% TODO: Variable TZoneSetLow not found in Data/preComputedInfrax.mat

if nargin < 3
   RefsParam.Price.variable = 0; %1 =  variable price profile, 0 = fixed to 1
end

if nargin == 0
   model.buildingType = 'Reno';  
   model.pred.ny = 6;
%    addpath('../common_files_all_sims/Data')
end

% TODO: generalize this temporary fix - use Te for infrax
if strcmp(model.buildingType,'Infrax') 
   model.pred.ny = 19; 
   model.buildingType = 'Reno';
end
  

bui_path = ['../buildings/', model.buildingType, '/disturbances/'];
fprintf('*** Load references ... \n')



if  strcmp(model.buildingType,'HollandschHuys')

    %% cofort zone based on standards ISO7730 and EN15251
    % for more details see Damians PhD page 168

    % Te(i)  - amient temperature at timestep i
    Te_index = 212;    % TODO: generalize ambient temperature indexing for all models in BuiSim
    Te = dist.d(:,Te_index) - 273.15;
    day_steps = 86400/model.pred.Ts;     % number of steps per day,  24h = 86400 sec
    steps = length(Te);                  % number of datapoints in dataset
    days = floor(steps/day_steps);  %  number of days in dataset

    % Tem(k) - one day average  amient temperature at day k
    for k = 0:days-1
       Tem(k+1) =  mean(Te(1+k*day_steps:(1+k)*day_steps));
    end
    % Tem = movmean(Te,day_steps);

    % TRM =  movmean(Tem,7);
    for k = 1:days
        % Trm(k) - moving mean ambient temperature at day k
        if k == 1
            Trm(k) = Tem(k);
        elseif k == 2
            Trm(k) = (Tem(k) + 0.8*Tem(k-1))/1.8;
        elseif k == 3
            Trm(3) = (Tem(k) + 0.8*Tem(k-1) + 0.6*Tem(k-2))/2.4;    
        elseif k == 4
            Trm(4) = (Tem(k) + 0.8*Tem(k-1) + 0.6*Tem(k-2) + 0.5*Tem(k-3))/2.9;       
        elseif k == 5
            Trm(5) = (Tem(k) + 0.8*Tem(k-1) + 0.6*Tem(k-2) + 0.5*Tem(k-3) + 0.4*Tem(k-4))/3.3;   
        elseif k == 6
            Trm(k) = (Tem(k) + 0.8*Tem(k-1) + 0.6*Tem(k-2) + 0.5*Tem(k-3) + 0.4*Tem(k-4) + 0.3*Tem(k-5))/3.6; 
        else
            Trm(k) = (Tem(k) + 0.8*Tem(k-1) + 0.6*Tem(k-2) + 0.5*Tem(k-3) + 0.4*Tem(k-4) + 0.3*Tem(k-5) + 0.2*Tem(k-6))/3.8;  
        end

        % upper and lower comfort bounds at individual days
        if Trm(k) < 10   % heating season comfort bounds
            TUp(k) = 24;
            TLow(k) = 20;
        elseif Trm(k) > 15  % cooling season comfort bounds
            TUp(k) = 26;
            TLow(k) = 23;
        else  % transient season comfort bounds
            TUp(k) = 24 + 2*(Trm(k)-10)/5;
            TLow(k) = 20+ 3*(Trm(k)-10)/5; 
        end
    end

    % repeat elements of vector to fit the sampling
    WB = repelem(TLow,day_steps)' + 273.15;
    WA = repelem(TUp,day_steps)' + 273.15;


    %%  night setbacks
    % thermal comfort start and end hour
    TCF_start_hour = 7;
    TCF_end_hour = 20;
    % hours to simsteps
    TCF_start_steps = TCF_start_hour*3600/model.plant.Ts;
    TCF_end_steps = TCF_end_hour*3600/model.plant.Ts;

    % therrmal comfort index of active timesteps
    TCF_index = [];
    for k = 0:days-1
       TCF_index = [TCF_index, k*day_steps+TCF_start_steps:k*day_steps+TCF_end_steps];
    end

    % night setback index
    NSB = ones(length(WA),1);
    NSB(TCF_index) = 0;
    % find([1:steps], TCF_index)


    % temperature comfort zone
    references.wb = WB*ones(1,model.pred.ny) - NSB*2;%  below threshold
    references.wa = WA*ones(1,model.pred.ny) + NSB*2; %  above threshold

    % PMV comfort zone
    references.PMVub(WB == min(WB)) = 1;
    references.PMVub(WA > min(WA)) = 0.5;
    references.PMVlb = -references.PMVub;


    %% variable price profiles
    if RefsParam.Price.variable
       references.Price = 1+sin(0.01*(1:length(WB)))'; % variable price profile
    %    references.Price = time_series;
    %    TODO:  load price profile interface
    else
       references.Price = ones(size(WB));  % standard fixed price 
    end

else
    %% comfort boundaries
    % function eval for full year comfort boundaries profiles in K
    [t_comf, TLow, TUp, TSup, TRefControl] = comfortTemperature(bui_path);

    % % visualisation of comfort constraints and reference
    % plot(t_comf,[TLow, TUp, TLow+(TUp-TLow)/2])
    % legend('TLow','TUp','ref')
    % figure
    % plot(t_comf,TSup)
    % legend('supply watter')

    references.R = (TLow+(TUp-TLow)/2)*ones(1,model.pred.ny);  % setpoint in K
    % references.ref = (TRefControl+2.5)*ones(1,model.pred.ny);  % setpoint in K
    references.wb = TLow*ones(1,model.pred.ny); %  below threshold
    references.wa = TUp*ones(1,model.pred.ny); %  above threshold
    % PMV bonds
    references.PMVub(TLow == min(TLow)) = 1;
    references.PMVub(TLow > min(TLow)) = 0.5;
    references.PMVlb = -references.PMVub;

    references.TSup = TSup;   % supply water temperature from HC

    %% variable price profiles
    if RefsParam.Price.variable
       references.Price = 1+sin(0.01*(1:length(TLow)))'; % variable price profile
    %    references.Price = time_series;
    %    TODO:  load price profile interface
    else
       references.Price = ones(size(TLow));  % standard fixed price 
    end
end

fprintf('*** Done.\n')
end