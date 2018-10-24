function MLagent = BuiLearn(MLagent,traindata)



if MLagent.TDNN.use
    fprintf('\n%s\n', repmat('-', 1, 40));
    fprintf('Time delayed neural network:\n');
    
    
%     end_time = 1000; %  NN works only with maximum 10000 samples 
    start_t = clock;
%     if size(traindata.features,1) > end_time 
%         net_ts = NN_TS(traindata.features(1:end_time,:),traindata.targets(1:end_time,:),MLagent.numDelays); 
%     else
%         net_ts = NN_TS(traindata.features,U,model.ml.numDelays);  
%     end

    % data for regression
    x = traindata.features';  % features
    t = traindata.targets';   % regression outputs
    % data for time series
    X = con2seq(x);
    T = con2seq(t);
% Learning parameters
    MLagent.TDNN.net.divideFcn = '';
    MLagent.TDNN.net.trainParam.showWindow = true;
    %             net.trainParam.showWindow = 0;
    %             net.trainParam.showCommandLine = 1;
    % % % % nr of iterations
    MLagent.TDNN.net.trainParam.epochs = 4000;
    % % % % nr of validation checks
    MLagent.TDNN.net.trainParam.max_fail = 100;
    [inputs,inputStates,layerStates,targets] = preparets(MLagent.TDNN.net,X,T);   % delay setup 
    % [net,tr] = train(net,inputs,targets);
    % [net,tr] = train(net,inputs,targets,inputStates,layerStates,'reduction',100);  % memory reduction option
    [MLagent.TDNN.net,tr] = train(MLagent.TDNN.net,inputs,targets,inputStates,layerStates);
    % CHECKPOINTS
    % http://www.mathworks.com/help/nnet/ug/checkpoint-saves-during-neural-network-training.html
    outputs = MLagent.TDNN.net(inputs,inputStates);
    errors = gsubtract(targets,outputs);
    performance = perform(MLagent.TDNN.net,targets,outputs)
    % Remove a delay from the network, to get the prediction one time step earlier.
    MLagent.TDNN.net = removedelay(MLagent.TDNN.net);
    % make net for simulations = closed-loop (parallel) configuration
    net_cl = closeloop(MLagent.TDNN.net);

    
    
    
    
    
    train_time = etime(clock, start_t);  
    fprintf('\nTime delayed neural network training time: %.2f sec:\n',train_time);

    
    
elseif MLagent.RT
    
    
end



end