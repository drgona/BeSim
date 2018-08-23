function [t, v, inputIndex, dictCtlInputs, dicValVar, dicOutputNameIndex, x0] = disturbances(pathBui,disturbanceType,flagPlot, debug)
    if nargin < 2 
        pathBui = '../buildings/Infrax';
        disturbanceType = ''; 
        flagPlot = true;
        debug = false;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INCLUDE VALIDATION
    include_validation = false;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    path_preComp = [pathBui '/disturbances/preComputed' disturbanceType '.mat'];
    path_unames = [pathBui '/models/uNames.txt'];
    path_ynames = [pathBui '/models/yNames.txt'];
    path_xnames = [pathBui '/models/xNames.txt'];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get disturbance variable names from uNames of model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Read variables names
    fileID = fopen(path_unames,'r');
    content = textscan(fileID,'%s','delimiter','\n');
    fclose(fileID);
    
    % find variables for disturbance
    dNamesInInputs = {'winBusIn','weaBus','prescribedBus'};
    dNamesInPreComputed = {'winBusOut','weaBusOut','prescribed'};
    inputIndex = [];
    for i=1:length(dNamesInInputs)
        % find string with keyword
        a = strfind(content{1},dNamesInInputs{i});
        % find corresponding index
        ind = find(~cellfun(@isempty,a));
        inputIndex = [inputIndex; ind];
    end
    % Replace input names by disturbance name
    varNames = content{1}(inputIndex);
    for i=1:length(inputIndex)
        for j=1:length(dNamesInInputs)
            varNames{i} = strrep(varNames{i},dNamesInInputs{j},dNamesInPreComputed{j});
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get initial value of all states
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Open xNames.txt to know number of states
    fileIDx = fopen(path_xnames,'r');
    contentx = textscan(fileIDx,'%s','delimiter','\n');
    fclose(fileIDx);
    nx = numel(contentx{1});
    % read initial state values
    x0 = zeros(nx,1);
    if include_validation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get validation variable names from uNames of model
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dicValVar = findMultipleOutputs(path_preComp,{'propsBusVal'},1);

        for i = 1:nx
            key = ['propsBusVal.TStaInit[' num2str(i) ']'];
            tmp = dicValVar(key);
            x0(i)=tmp(1);
        end
    else
        dicValVar = 0;
        x0 = ones(nx,1)*293.15;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save control inputs names and indexes into dictionary from uNames of model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dNamesControlInputs = {'inputBus'};
    keys = [];
    indices = [];
    for i=1:length(dNamesControlInputs)
        % find strings with keyword
        k = strfind(content{1},dNamesControlInputs{i});
        % find corresponding index
        ind = find(~cellfun(@isempty,k));
        indices = [indices; ind];
        keys = [keys; content{1}(ind)];
    end
    dictCtlInputs = containers.Map(keys, indices);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get disturbance variable values
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    t = findOutput(path_preComp, 'Time' )';
    v = zeros(length(t),length(varNames));
    for i=1:length(varNames)
        if debug
            temp1 = strfind(varNames(i),'Tenv');
            temp2 = strfind(varNames(i),'Te');
            temp3 = strfind(varNames(i),'TGroundDes');
            temp4 = strfind(varNames(i),'dummy');
            if ~isempty(temp1{1}) || ~isempty(temp2{1}) || ~isempty(temp3{1})
                v(:,i) = ones(length(t),1).*293.15;
            elseif ~isempty(temp4{1})
                v(:,i) = ones(length(t),1);
            else
                v(:,i) = zeros(length(t),1);
            end
        else
            v(:,i) = findOutput(path_preComp,varNames(i))';
        end
    end
    %Remove last point because not equidistant
    t = t(1:end-1);
    v = v(1:end-1,:); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get output names and indices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Read variables names
    fileID = fopen(path_ynames,'r');
    contenty = textscan(fileID,'%s','delimiter','\n');
    fclose(fileID);   
    dicOutputNameIndex = containers.Map(contenty{1}, 1:numel(contenty{1}));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot disturbances
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flagPlot
        indISolDir = find(~cellfun(@isempty,strfind(varNames,'iSolDir')));
        indISolDif = find(~cellfun(@isempty,strfind(varNames,'iSolDif')));
        figure()
        plot( t*ones(1,length(indISolDir)) , v(:,indISolDir) ), hold on
        plot(t,v(:,indISolDif)), hold off
        legend([varNames(indISolDir); varNames(indISolDif)])
        xlabel('t [s]')
        ylabel('Direct and diffuse [W/m2]')
        
        indTEnv = find(~cellfun(@isempty,strfind(varNames,'Tenv')));
        temp = find(~cellfun(@isempty,strfind(varNames,'Te')));
        temp2 = find(cellfun(@isempty,strfind(varNames(temp),'Tenv')));
        indTe = temp(temp2);
        figure()
        plot(t,v(:,indTEnv)), hold on
        plot(t,v(:,indTe),'k','linewidth',2), hold off
        legend([varNames(indTEnv); varNames(indTe)])
        xlabel('t [s]')
        ylabel('Environment and ambient temperature [K]')
        
        
        indQAbsAndISolWin = find(~cellfun(@isempty,strfind(varNames,'AbsQFlow')));
        figure()
        plot(t,v(:,indQAbsAndISolWin)), hold on
        legend(varNames(indQAbsAndISolWin))
        xlabel('t [s]')
        ylabel('Absorbed, direct and diffuse heat through window [W]')
        
        indQGai = find(~cellfun(@isempty,strfind(varNames,'prescribed.Qprescribed')));
        if length(indQGai)>1
            figure()
            plot(t,v(:,indQGai)), hold on
            legend(varNames(indQGai))
            xlabel('t [s]')
            ylabel('Radiative and convective gains [W]')
        end
        
        indTBou = find(~cellfun(@isempty,strfind(varNames,'prescribed.TBou')));
        if length(indTBou)>1
            figure()
            plot(t,v(:,indTBou)-273.15), hold on
            legend(varNames(indTBou))
            xlabel('t [s]')
            ylabel('Comfort temperatures [degC]')
        end
        
        indMVen = find(~cellfun(@isempty,strfind(varNames,'prescribed.m_flow_ven')));
        if length(indTBou)>1
            figure()
            plot(t,v(:,indMVen)), hold on
            legend(varNames(indMVen))
            xlabel('t [s]')
            ylabel('Mass flow ventilation [kg/s]')
        end
        
        if include_validation
            val = cell2mat(dicValVar.values');
            figure()
            plot(t,val-273.15), hold on
            legend(dicValVar.keys)
            xlabel('t [s]')
            ylabel('Validation Temperature [degC]')

            figure()
            plot(x0-273.15), title('Initial state values'), ylabel('[degC]')
            set(gca, 'XTick', 1:1:nx, 'XTickLabel', contentx{1},'XTickLabelRotation',45);
        end
    end
    
end