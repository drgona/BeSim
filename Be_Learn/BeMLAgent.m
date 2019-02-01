function MLagent = BeMLAgent(outdata, AgentParam)

% todo: generate agent untrained structure

% initialization
MLagent.RT.use = 0;
MLagent.regNN.use = 0;
MLagent.TDNN.use = 0;
MLagent.numDelays = outdata.ctrl.MPC.N;


fprintf('\n------------------ Constructing ML Agent -----------------\n');


if AgentParam.RT.use
    MLagent.RT.use = 1;
    MLagent.RT.QuadraticErrorTolerance = 1e-4;
    MLagent.numDelays = 2;  % shortened delayed features
    MLagent.type = 'Regression tree'; 

elseif AgentParam.regNN.use
     MLagent.regNN.use = 1;
     MLagent.regNN.hiddenLayerSize = [24, 12];  
     MLagent.regNN.trainFcn =  'trainscg';  
     MLagent.regNN.net = fitnet(MLagent.regNN.hiddenLayerSize,MLagent.regNN.trainFcn);   % initialize network
     % change  layer activation functions to ReLU  
%      MLagent.regNN.net.layers{1:end-1}.transferFcn =  'poslin';
     MLagent.type = 'Regression neural network'; 
     
elseif AgentParam.TDNN.use
    MLagent.TDNN.use = 1;
    MLagent.TDNN.hiddenLayerSize = [48,24,12];
    MLagent.TDNN.trainFcn =  'trainscg';  
    MLagent.TDNN.net = timedelaynet([1:MLagent.numDelays],MLagent.TDNN.hiddenLayerSize,MLagent.TDNN.trainFcn);  
    % change hidden layer activation functions to ReLU  
%     MLagent.TDNN.net.layers{1:end-1}.transferFcn =  'poslin';
    MLagent.type = 'Time delay neural network'; 
end

fprintf('*** Agent Type = %s\n', MLagent.type);
fprintf('*** Number of time delays = %d\n', MLagent.numDelays );
fprintf('*** Done.\n') 

end
