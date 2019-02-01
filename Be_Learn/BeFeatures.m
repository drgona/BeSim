function [traindata, MLagent] = BeFeatures(outdata, dist, MLagent, FeaturesParam)

% todo: finish


if nargin < 3
%   normalized  precisions on principal components and features  0 to 1
    % parameters for feature refuction function
    FeaturesParam.reduce.PCA.use = 1;
    FeaturesParam.reduce.PCA.component = 0.999;   % principal component weight threshold
    FeaturesParam.reduce.PCA.feature = 0.95;      % PCA features weight threshold
    FeaturesParam.reduce.D_model.use = 1;           %  
    FeaturesParam.reduce.D_model.feature = 0.99;    % model features weight threshold
    FeaturesParam.reduce.lincols.use = 1;
    FeaturesParam.reduce.flagPlot = 1;
end
    
fprintf('\n------------------ Constructing Features -----------------\n');
fprintf('*** Use PCA reduction = %d\n', FeaturesParam.reduce.PCA.use);
fprintf('*** Use disturbance model reduction = %d\n', FeaturesParam.reduce.D_model.use);
fprintf('*** Use linearly dependent columns reduction = %d\n', FeaturesParam.reduce.lincols.use);


%% prepare data
% ====== FEATURES FOR current time regressions ======
X = outdata.data.X(:,1:end-1)';  % X - feature (states) 
Y = outdata.data.Y';  % Y - feature (outputs)
U = outdata.data.U';  % U - targets (inputs)
D = outdata.data.D(:,1:end-outdata.ctrl.MPC.Ndp)';  % D - feature (disturbances)
wa = outdata.data.wa(1,1:end-outdata.ctrl.MPC.Nrp)';  % wa - feature (above threshold) 
wb = outdata.data.wb(1,1:end-outdata.ctrl.MPC.Nrp)';  % wa - feature (below threshold)

% ====== FEATURES FOR predictive regressions ======
% disturbances and thresholds are shifted N steps back to be suitable for time delayed NN and TS regressions
D_pred = outdata.data.D(:,1+(MLagent.numDelays-1):end-(outdata.ctrl.MPC.Ndp-MLagent.numDelays+1))';
wa_pred = outdata.data.wa(1,1+(MLagent.numDelays-1):end-(outdata.ctrl.MPC.Ndp-MLagent.numDelays+1))';
wb_pred = outdata.data.wb(1,1+(MLagent.numDelays-1):end-(outdata.ctrl.MPC.Ndp-MLagent.numDelays+1))';
% D_pred = outdata.data.D(:,1+MLagent.numDelays:end-(outdata.ctrl.MPC.Ndp-MLagent.numDelays))';
% wa_pred = outdata.data.wa(1,1+MLagent.numDelays:end-(outdata.ctrl.MPC.Ndp-MLagent.numDelays))';
% wb_pred = outdata.data.wb(1,1+MLagent.numDelays:end-(outdata.ctrl.MPC.Ndp-MLagent.numDelays))';


%% ====== FEATURES elimination ======
% elimination of disturbance variables
[D_use, D_discard, MLagent.use_D] = FeatureReduce(D,outdata,dist,FeaturesParam.reduce);
% [D_pred_use, D_pred_discard, MLagent.use_disturb] = FeatureSelect(D_pred,model,PrecisionParam,plotFlag,UseParam);

% TODO: elimination of state variables

%% Datasets construction
% rows = samples,  columns = feature variables


% %     TODO: FIX THIS delays   


%% Dataset selection
% TODO: select datasets based on MLagent
if MLagent.RT.use
    
    if FeaturesParam.time_transform
        sine_waves = time_transform(size(XX4,1),model.plant.Ts);
        XX6 = [XX4 sine_waves];
    else
%         XX4
    end
    
elseif MLagent.regNN.use
    %  XX4 suitable for regNN with predictive behavior 
    flipFlag = false;
    % create delayed features
    D_pred_delay = delayData(D_pred(:,MLagent.use_D),MLagent.numDelays,flipFlag);
    wb_pred_delay = delayData(wb_pred,MLagent.numDelays,flipFlag);
    % TODO: check pefrormance with delayed Y for  NN reg
    Y_delay = delayData(Y,MLagent.numDelays,flipFlag);
    %  suitable for TS regression trees
    traindata.features  = [Y(1:end-(MLagent.numDelays-1),:) D_pred_delay wb_pred_delay];
    % delayed response
    traindata.targets = U(MLagent.numDelays:end,:);
    
    
elseif MLagent.TDNN.use
    %  suitable for delayed TS NN with predictive behavior 
    % features with laged predictions for D and w  
    traindata.features = [Y D_pred(:,MLagent.use_D)  wb_pred]; 
    traindata.targets = U; 
end

fprintf('*** Done.\n') 

% TODO: generate index vector for feature selection from parametric space
% for MPC???


end
