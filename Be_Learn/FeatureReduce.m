function [Xreduced, Xdiscard, use_feature] = FeatureReduce(X,outdata,dist,ReduceParam)

% [Xreduced, Xdiscard, use_feature] = FeatureSelect(X,model,PrecisionParam,ReduceParam.flagPlot,UseParam)
% -----------INPUTS-------------------------------
%     X     - features matrix,  rows = samples,  columns = features
%     outdata - output data from simulation
% -----------PRECISIONS----------------------------
%     PrecisionParam.prec_pca           - PCA precision
%     PrecisionParam.prec_pcafeature    - influence of the features on PCA
%     PrecisionParam.prec_model         - influcene of the disturnances on the model dynamics
% -----------METHODS-------------------------------
%     1, ReduceParam.PCA.use       - Principal component analysis
%     2, ReduceParam.model.use   - disturbance effect on model dynamics via E matrix
%     3, ReduceParam.lincols.use   - linear dependency of the features
% -----------USE--------------------------------------
%     if input data X are disturbances than ReduceParam.model.use = true
%     otherwise, for more general data ReduceParam.model.use = false
% ---------------------------------------------

if nargin < 3
%   normalized  precisions on principal components and features  0 to 1
    % parameters for feature refuction function
    ReduceParam.PCA.use = 1;
    ReduceParam.PCA.component = 0.999;   % principal component weight threshold
    ReduceParam.PCA.feature = 0.95;      % PCA features weight threshold
    ReduceParam.D_model.use = 1;
    ReduceParam.D_model.feature = 0.99;   % model features weight threshold
    ReduceParam.lincols.use = 1;
    ReduceParam.flagPlot = 1;
end
    
nx = size(X,2); % number of features

 
%% Discard LINEARLY DEPENDENT COLUMNS
if ReduceParam.lincols.use
     [~,index_lincols]=licols(X);   
     use_feature_lincols = ismember(1:nx,index_lincols); % index of choosen features from lincols
else
     use_feature_lincols = ones(1,nx); % no features have been discarded
end
 
%% PRINCIPAL COMPONENT ANALYSIS
% http://stackoverflow.com/questions/19723993/which-variables-combine-to-form-most-of-the-variance-for-a-principal-component-i
% https://www.mathworks.com/matlabcentral/answers/49134-determining-variables-that-contribute-to-principal-components
 % http://www.mathworks.com/help/stats/pca.html
if ReduceParam.PCA.use
    [coeff,~,~,~,explained] = pca(X);  % Each column of coeff contains coefficients for one principal component.
    nu_pca = sum(explained > 100*(1-ReduceParam.PCA.component)); % number of relevant principal components
    var_weight = sum(abs(coeff(:, 1:nu_pca)')); % weight of the features
    var_weight = var_weight/max(var_weight);    % normalized weight of the features
    use_feature_PCA = var_weight >= 1-ReduceParam.PCA.feature; % indexdes of choosen features from PCA
else
    use_feature_PCA = ones(1,nx); % no features have been discarded
end

%% disturbance dynamics analysis  
%% TODO: rework disturbance dynamics analysis based on elementwise operations on A,E matrix and max(abs) of X,D or max(abs(diff)) of X,D
   % coefficients in E MATRIX show which disturbances are most influtential
   % max(model.dist.d,[],1)  = largest elements in disturbances or maybe use average values in the disturbances vector
if ReduceParam.D_model.use       
        influecneD = (sum(abs(outdata.model.plant.Ed),1).*max(dist.d,[],1));
        influecneD = influecneD/max(influecneD);  % normalized influecne of each disturbance in model dynamics  
        use_feature_D_model = influecneD > 1-ReduceParam.D_model.feature;   % index of choosen features from model dynamics  
else
        use_feature_D_model = ones(1,nx); % no features have been discarded
end

% if ReduceParam.X_model.use       
%         influecneX = (sum(abs(outdata.model.plant.Ad),1).*max(outdata.data.X',[],1));
%         influecneX = influecneX/max(influecneX);  % normalized influecne of each disturbance in model dynamics  
%         use_feature_X_model = influecneD > 1-ReduceParam.X_model.feature;   % index of choosen features from model dynamics  
% else
%         use_feature_X_model = ones(1,nx); % no features have been discarded
% end
        

% % if onlyPCA
%     use_feature = use_feature_PCA;
%     Xreduced = X(:,use_feature); % choosen features
%     Xdiscard = X(:,not(use_feature)); % discarded features
% %     TODO: incorporate linear dependency test also here!!!!
% % else 
%         strongD = influecneD > 1-PrecisionParam.prec_model;  
%         
%         % most significant disturbances based on PCA and E matrix
%         Dsignificant = find(use_feature_PCA + strongD >1);
%         
% %       removing lineary dependent rows
%         [Xreduced,use_feature_model]=licols(X(:,(use_feature_PCA + strongD >1)));      
% %         manual extraction for original model
% %         use_feature_model = setdiff(Dsignificant, [31 34 37]);   % 31 34 37 40 disturbances have the same profiles 
% 
% %         OUTPUTS from PCA + model analysis
% %         use_feature = Dsignificant;  % return significant features   without removing the linearly dependent columns   
%         use_feature = Dsignificant(use_feature_model);
%         Xdiscard = X(:,setdiff(1:size(X,2), use_feature)); % discarded features  

%% OUTPUTS 
use_feature = use_feature_D_model + use_feature_PCA + use_feature_lincols > 2; % index of choosen features
Xreduced = X(:,use_feature); % choosen features
Xdiscard = X(:,not(use_feature)); % discarded features

%% PLOTS     

% %  TODO: visualize thresholds
%  visualize A, E matrix coefficients

if ReduceParam.flagPlot

            figure         
            subplot(2,2,1)
            plot(Xreduced, 'linewidth', 2)
            title('Used features profiles')
            axis tight
            xlabel('samples []')
            set(gca,'fontsize',20)
            grid on
            
            subplot(2,2,2)
            plot(Xdiscard, 'linewidth', 2)
            title('Discarded features profiles')
            axis tight
            xlabel('samples []')
            set(gca,'fontsize',20)
            grid on
            
            %             subplot(2,2,2)
%             bar(explained)  
%             title('Importance of the principal components')
%             xlabel('Principal components')
%             ylabel('Importance [%]')
            % figure      % effect of the each variable on first nu_pca principal components 
            % bar(coeff(:, 1:nu_pca)) 
     
          if ReduceParam.PCA.use
            subplot(2,2,3)
            bar(var_weight)
            title('Importance of the features on principal components')
            xlabel('Features')
            ylabel('Relative importance [0-1]')
            set(gca,'fontsize',20)
            grid on
          end         

          if ReduceParam.D_model.use
            subplot(2,2,4)
            bar(influecneD)
            title('Importance of the disturbances in the model')
            xlabel('Disturbances (Features)')
            ylabel('Relative importance [0-1]')  
            set(gca,'fontsize',20)
            grid on
          end 
end
    
% % additional plots - presentation
% figure
% subplot(2,2,1)
% imagesc(coeff')
%             title('PCA coefficients')
%             xlabel('Features')
%             ylabel('Principal components') 
%             colorbar
%             caxis([-0.5, 1])
%     subplot(2,2,2) 
% % figure
%         barh(explained)    
%              title('Variance of the principal components')
%             ylabel('Principal components')
%             xlabel('Total variance percentage [%]')  
% %             set ( gca, 'ydir', 'reverse' )
%   subplot(2,2,3)
% % figure
% imagesc(coeff(:, 1:nu_pca)')              
%             title('Most significant PCA coefficients')
%             xlabel('Features')
%             ylabel('Principal components')   
%             colorbar
%             caxis([-0.5, 1])
% subplot(2,2,4)  
% % figure
% bar(var_weight)
%                title('Feature weights on most influential components')
%             xlabel('Features')
%             ylabel('Normalized weight')         

% 
% figure
% imagesc((outdata.model.plant.Ed))
% colorbar
% title('Disturbance matrix E')
% xlabel('Disturbances')
% ylabel('States')         

%% ADDITIONAL COMMENTS or directions

% % % %  ============ PCA comments ==============
% % % %  SOME EXPERIMENTS WITH PCA - if used adopt feature select function
% [coeff,score,~,tsquared,explained] = pca(XX3);
% tsquared_norm = tsquared/max(tsquared);
% % figure
% % bar(tsquared_norm)
% not_use = tsquared_norm < 0.01;
% % bar(tsquared_norm(not(not_use)))
% % bar(tsquared_norm((not_use)))
% XX3_reduced =  XX3(not(not_use),:);
% U_reduced = U(not(not_use),:);

% TODO: use tsquared - for analysis of the importance of the data samples
% for possible redictuction of dataset  = not beneficial


% % % %  ============ model comments ==============
        % % which states are influenced most by the disturbances?
        % bar(sum(abs(outdata.model.plant.Ed),2))
        % % what is the effect of states on states?
        % bar(sum(abs(outdata.model.plant.Ad),1))
        % todo: select the states which are affected  at most most by the disturbances
    

% % % %  ============ other comments ==============
% TODO: % % % % investigate other feature selection methods
% https://www.mathworks.com/help/stats/feature-selection.html
% https://www.mathworks.com/discovery/feature-selection.html
% https://www.mathworks.com/help/stats/feature-transformation.html#f72219

% % CHECK OUT THESE FUNCTIONS
% sequentialfs	Sequential feature selection
% relieff	Importance of attributes (predictors) using ReliefF algorithm
% stepwiselm	Create linear regression model using stepwise regression
% stepwiseglm	Create generalized linear regression model by stepwise regression
% fscnca	Feature selection using neighborhood component analysis for classification
% fsrnca	Feature selection using neighborhood component analysis for regression



end