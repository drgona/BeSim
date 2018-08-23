function run_simulation(reload, nStep, saveFig, disturbanceType)
	if nargin < 2
        addpath(pwd);
		cd('../buildings/Infrax/');
        reload = 0;
        nStep = 400;
        saveFig = 0;
        disturbanceType = ''; % can be '_lin' if used for linearization validation
    elseif nargin < 3
        saveFig = 0;
    end;

	path = pwd;
    splitted_path = strsplit(path,'\'); buiName = splitted_path{end}; path_fig = 'Fig/';

	%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Load disturbances and model 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if reload
        fprintf('*** Load disturbances ... \n')
		% Disturbances
		[t, v, inputIndex, dictCtlInputs, dicValVar, dicOutputNameIndex, x0]= disturbances(path,disturbanceType,0, 0);
		save(['preComputed_matlab/dis' disturbanceType '.mat'], 't', 'v', 'inputIndex', 'dictCtlInputs', 'dicValVar','dicOutputNameIndex', 'x0');
		fprintf('*** Done.\n')
        
		% Model
		fprintf('*** Create ROM model ...\n')
        Ts = t(2) - t(1);
		orders = []; % add any reduced order model you wish to have
		[sys_dExt, rom] = fGenerateSysAndRom([path '/models/ssm.mat'], Ts, x0, orders);      
		save('preComputed_matlab/mod.mat', 'Ts', 'orders', 'sys_dExt', 'rom');
   		fprintf('*** Done.\n')
    else
        fprintf('*** Load disturbances and ROM model...\n')
		load(['preComputed_matlab/dis' disturbanceType '.mat']);
		load('preComputed_matlab/mod.mat');
        fprintf('*** Done.\n')
	end;


	%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Simulate
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('*** Simulate SSM model ...\n')
	% Parameters
	t_ind_start = 1;
	index_t = t_ind_start:1:t_ind_start+nStep-1;
	mod = sys_dExt;

	%%%% Input vector (only disturbances, zero inputs).
	nu = size(mod.B,2);
	% empty inputs
	u = zeros(nStep,nu);
	% Add input = 1 for initialization
	u(:,end) = ones(nStep,1);
	
	% Disturbances
	u(:,inputIndex) = v(1:nStep,:);

	% Simulate
	y_ssm = lsim(mod, u, t(index_t));
    fprintf('*** Done.\n')

	% ---- Plot operative temperatures
    fprintf('*** Plot results ...\n')
	% find output index and modelica validation data
	indy = []; %index of output for ssm
    nZones = sum(count(dicOutputNameIndex.keys,'outputBus.TZones['));
	legend_ssm={};
    for i_zone = 1:nZones
        key = ['TZones[' num2str(i_zone) ']'];
        indy = [indy; dicOutputNameIndex(['outputBus.' key])];
        legend_ssm{i_zone}=['ssm.' key];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    positionFig = [100, 100, 300, 200];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
	% plot time series
	fig1 = figure('Name',[buiName '-1'],'Position', positionFig);
	plot(t(index_t),y_ssm(:,indy)), title('SSM simulation')
	legend(legend_ssm{:})
    if saveFig
        set(fig1,'Units','Inches');
        pos = get(fig1,'Position');
        set(fig1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(fig1,[path_fig buiName '-ModSSMCompSim' disturbanceType],'-dpdf','-r0');
        %print(fig2,[path_fig buiName '-ModSSMCompSim' disturbanceType],'-depsc','-r0');
    end;
end
