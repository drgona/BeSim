function outdataSampled = BeSample(outdata,SampleParam,PlotParam)

if nargin < 2 
    SampleParam.Random.use = 0; % 1, random sampling along the optimal trajectory
    SampleParam.StateInit = 1; % 2, offset state initialization per day  
end

if nargin < 3
    PlotParam.flagPlot = 1;     % plot 0 - no 1 - yes
    PlotParam.plotStates = 0;        % plot states
    PlotParam.plotDist = 0;        % plot disturbances
    PlotParam.plotEstim = 0;        % plot estimation
    PlotParam.plotCtrl = 1;        % plot control
    PlotParam.plotPrice = 0;        % plot price signal
end

if not(SampleParam.use) % DO NOT USE dataset extension
    outdataSampled = outdata;  
    
else    % USE dataset extension    
    % load disturbances and references
    DistParam.reload = 0;
    dist = BeDist(outdata.model, DistParam);        % construct a disturbances object  
    RefsParam.Price.variable = 0;       %1 =  variable price profile, 0 = fixed to 1
    refs = BeRefs(outdata.model, dist, RefsParam);     % construct a references object
%     initialization of a new sampling simulation
    if not(SampleParam.continue.use)
        % initialize output matrices
        outdataSampled.data.X = outdata.data.X(:,1:end-1);         %  state vector
        outdataSampled.data.Y = outdata.data.Y;         %  output vector
        outdataSampled.data.U = outdata.data.U;         %  input vector
        outdataSampled.data.D = outdata.data.D(:,1:end-outdata.ctrl.MPC.Ndp);         %  disturbance vector
        outdataSampled.data.wa = outdata.data.wa(:,1:end-outdata.ctrl.MPC.Nrp);       %  above threshold
        outdataSampled.data.wb = outdata.data.wb(:,1:end-outdata.ctrl.MPC.Nrp);       %  below threshold
        % outdata structures
        outdataSampled.model = outdata.model;
        outdataSampled.estim = outdata.estim;
        outdataSampled.ctrl = outdata.ctrl;
        outdataSampled.SimParam = outdata.SimParam;
    else
        outdataSampled = outdata;  
    end
    
    %%  MPC for a range of state initializations of the simulation dataset
    if SampleParam.StateInit
        %  plus minus offset initialization of states per day
        StateInitRange = 1:4;  
        % state initializing of the simulation dataset per day
        if not(SampleParam.continue.use)
            SimStart = outdata.SimParam.run.start;    
        else
            SimStart = SampleParam.continue.i;   
        end
        SimEnd = outdata.SimParam.run.end;
        SimDays = 1+SimEnd-SimStart;

        for i = 1:SimDays % iterate through days
            outdata.SimParam.run.start = SimStart+i-1;
            outdata.SimParam.run.end = SimStart+i-1;           
            for k = 1:2*length(StateInitRange)     %  iterate through state initialization offsets 
                fprintf('\n*** Iteration: day %d/%d, case %d/%d ***', outdata.SimParam.run.start, SimEnd,k,2*length(StateInitRange));
                if k<=length(StateInitRange)  % warm days
                     Offset = StateInitRange(k);          
                else  % cold days
                     Offset = -StateInitRange(k-length(StateInitRange));             
                end
                %   SIMULATE
                outdata.model.analyze.XinitOffset = Offset;
                outdataNew{i}{k} = BeSim(outdata.model, outdata.estim, outdata.ctrl, dist, refs, outdata.SimParam);   
                fprintf('*** State initialization offset: %.1f K\n', Offset);            
                %  plot
                if PlotParam.flagPlot
                     BePlot(outdataNew{i}{k},PlotParam)
                end             
                % concatenate data
                outdataSampled.data.X = [outdataSampled.data.X outdataNew{i}{k}.data.X(:,1:end-1)];  %  state vector
                outdataSampled.data.Y = [outdataSampled.data.Y outdataNew{i}{k}.data.Y];             %  output vector
                outdataSampled.data.U = [outdataSampled.data.U outdataNew{i}{k}.data.U];             %  input vector
                outdataSampled.data.D = [outdataSampled.data.D outdataNew{i}{k}.data.D(:,1:end-outdata.ctrl.MPC.Ndp)];             %  disturbance vector
                outdataSampled.data.wa = [outdataSampled.data.wa outdataNew{i}{k}.data.wa(:,1:end-outdata.ctrl.MPC.Nrp)];            %  above threshold
                outdataSampled.data.wb = [outdataSampled.data.wb outdataNew{i}{k}.data.wb(:,1:end-outdata.ctrl.MPC.Ndp)];            %  below threshold
            end
            % temporary save
%             save outdataSampled2.mat outdataSampled i k

            % TODO: comment if unused
            save('outdataSampled2.mat', 'outdataSampled', 'i', 'k', '-v7.3')
%             pause(120)  % cooling phase
        end
        % addinng forecast values for data structure compatibility
        outdataSampled.data.X = [outdataSampled.data.X outdataNew{i}{k}.data.X(:,end)];  
        outdataSampled.data.D = [outdataSampled.data.D outdataNew{i}{k}.data.D(:,end-outdata.ctrl.MPC.Ndp+1:end)];             %  disturbance vector
        outdataSampled.data.wa = [outdataSampled.data.wa outdataNew{i}{k}.data.wa(:,end-outdata.ctrl.MPC.Nrp+1:end)];            %  above threshold
        outdataSampled.data.wb = [outdataSampled.data.wb outdataNew{i}{k}.data.wb(:,end-outdata.ctrl.MPC.Nrp+1:end)]; 

%         TODO: create batches of training datasets due to the time
%         shifts - separate prepare data and reduction
        
    %%  MPC with a randomized state initialization of the optimal trajectories  
    elseif SampleParam.StateRandom
    % random sampling of states alongside optimal trajectory


    end
end
end