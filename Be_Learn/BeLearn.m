function MLagent = BeLearn(MLagent,traindata)



% TODOO: if batch size for TDNN limited use consecutive batches of smaller
% size to update the network


fprintf('\n------------------ Training ML Agent -----------------\n');

if MLagent.TDNN.use
    fprintf('\n%s\n', repmat('-', 1, 40));
    fprintf('*** Time delayed neural network:\n');
    fprintf('*** Training in progress\n');
    
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
    % data for time series - vector of cells with rows as cell elements
    X = con2seq(x);          
    T = con2seq(t);
    % Learning parameters
    MLagent.TDNN.trainFcn =  'trainscg';  
    MLagent.TDNN.net.divideFcn = '';
    MLagent.TDNN.net.trainParam.showWindow = true;
    MLagent.TDNN.net.trainParam.epochs = 500; % nr of iterations
    MLagent.TDNN.net.trainParam.max_fail = 100;  % nr of validation checks
    %  prepare time series data via preparets 
    [inputs,inputStates,layerStates,targets] = preparets(MLagent.TDNN.net,X,T);   % delay setup 
%     https://www.mathworks.com/help/deeplearning/ref/preparets.html
    % [MLagent.TDNN.net,tr] = train(MLagent.TDNN.net,inputs,targets,inputStates,layerStates,'reduction',100);  % memory reduction option
    [MLagent.TDNN.net,tr] = train(MLagent.TDNN.net,inputs,targets,inputStates,layerStates);
    % http://www.mathworks.com/help/nnet/ug/checkpoint-saves-during-neural-network-training.html
    outputs = MLagent.TDNN.net(inputs,inputStates);
    MLagent.TDNN.errors = gsubtract(targets,outputs);  % training set errors
    MLagent.TDNN.performance = perform(MLagent.TDNN.net,targets,outputs);   % mse metric
    
%     TODO: make this non-recursive, do it only once before implementation
    % Remove a delay from the network, to get the prediction one time step earlier.
    MLagent.TDNN.net = removedelay(MLagent.TDNN.net);
%     % make net for simulations = closed-loop (parallel) configuration
    MLagent.TDNN.net_apply = closeloop(MLagent.TDNN.net);
    % training finisher
    MLagent.TDNN.train_time = etime(clock, start_t);  
    fprintf('*** Training time: %.2f sec:\n',MLagent.TDNN.train_time);

    genFunction(MLagent.TDNN.net_apply,'../Be_Learn/TDNN_ctrl.m','MatrixOnly','yes')  % evaluation function generation
%     genFunction(MLagent.TDNN.net_apply,'../Be_Learn/TDNN_ctrl.m','MatrixOnly','yes')  % evaluation function generation

    
elseif MLagent.RT
    
    
end

fprintf('*** Done.\n')

end