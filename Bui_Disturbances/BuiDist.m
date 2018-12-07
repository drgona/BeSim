function dusturb = BuiDist(model, DistParam)

if nargin == 0
   model.buildingType = 'Infrax';  
end
if nargin < 2
   DistParam.reload = 0;
end


path = ['../buildings/', model.buildingType];
disturbanceType = ''; % can be '_lin' if used for linearization validation

if DistParam.reload
        fprintf('*** Load disturbances ... \n')
		% Disturbances
        if strcmp(model.buildingType,'Reno') || strcmp(model.buildingType,'RenoLight') || strcmp(model.buildingType,'Old') 
            [t, v, x0] = disturbances_old(path, 0, 0);
            save([path '/preComputed_matlab/dis' disturbanceType '.mat'], 't', 'v', 'x0');
        else
            [t,  v, inputIndex, dictCtlInputs, dicValVar, dicOutputNameIndex, x0]= disturbances(path,disturbanceType,0, 0);
            save([path '/preComputed_matlab/dis' disturbanceType '.mat'], 't', 'v', 'inputIndex', 'dictCtlInputs', 'dicValVar','dicOutputNameIndex', 'x0');
        end	
		fprintf('*** Done.\n')
else
        fprintf('*** Load disturbances...\n')
		load([path '/preComputed_matlab/dis' disturbanceType '.mat']);
        fprintf('*** Done.\n')
%         TODO:  create mod.mat file also for 6-zone building - connect
%         models in one file
end

        
%% disturbances
dusturb.t = t;
dusturb.d = v;
% max and min disturbances
dusturb.dmin = min(dusturb.d, [], 1);
dusturb.dmax = max(dusturb.d, [], 1);

%% TODO: automatic evaluation of statistical properties of datasets
% probability distributions, mean, min, max, etc.


end