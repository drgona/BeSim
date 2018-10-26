function [Xreduced, Xdiscard, use_feature] = FeatureSelect(X,model,dist,PrecisionParam,plotFlag,UseParam)

% [Xreduced, Xdiscard, use_feature] = FeatureSelect(X,model,PrecisionParam,plotFlag,UseParam)
% -----------INPUTS-------------------------------
%     X     - features matrix,  rows = samples,  columns = features
%     model - model structure generated via function b_model.m
% -----------PRECISIONS----------------------------
%     PrecisionParam.prec_pca           - PCA precision
%     PrecisionParam.prec_pcafeature    - influence of the features on PCA
%     PrecisionParam.prec_model         - influcene of the disturnances on the model dynamics
% -----------METHODS-------------------------------
%     1, UseParam.PCA       - Principal component analysis
%     2, UseParam.disturb   - disturbance effect on model dynamics via E matrix
%     3, UseParam.lincols   - linear dependency of the features
% -----------USE--------------------------------------
%     if input data X are disturbances than UseParam.disturb = true
%     otherwise, for more general data UseParam.disturb = false
% ---------------------------------------------

if nargin < 3
%   normalized  precisions on principal components and features  0 to 1
    PrecisionParam.PCAcomponent = 0.999;
    PrecisionParam.PCAfeature = 0.99;
    PrecisionParam.model = 0.999;
%     plotting flag
    plotFlag = false;
%     methods
    UseParam.PCA = true;
    UseParam.lincols = true;
    UseParam.disturb = true; 
end
    
if nargin < 2 
    UseParam.disturb = false;  % if there is no model info we can not perform disturbance analysis
end

nx = size(X,2); % number of features

 
%% USE EXCTRACTION OF LINEARLY INDEPENDENT COLUMNS
if UseParam.lincols
     [~,index_lincols]=licols(X);   
     use_feature_lincols = ismember(1:nx,index_lincols); % index of choosen features from lincols
else
     use_feature_lincols = ones(1,nx); % no features have been discarded
end
 
%% USE PRINCIPAL COMPONENT ANALYSIS
% http://stackoverflow.com/questions/19723993/which-variables-combine-to-form-most-of-the-variance-for-a-principal-component-i
% https://www.mathworks.com/matlabcentral/answers/49134-determining-variables-that-contribute-to-principal-components
 % http://www.mathworks.com/help/stats/pca.html
if UseParam.PCA
    [coeff,~,~,~,explained] = pca(X);  % Each column of coeff contains coefficients for one principal component.
    nu_pca = sum(explained > 100*(1-PrecisionParam.PCAcomponent)); % number of relevant principal components
    var_weight = sum(abs(coeff(:, 1:nu_pca)')); % weight of the features
    var_weight = var_weight/max(var_weight);    % normalized weight of the features
    use_feature_PCA = var_weight >= 1-PrecisionParam.PCAfeature; % indexdes of choosen features from PCA
else
    use_feature_PCA = ones(1,nx); % no features have been discarded
end

%% USE disturbance dinamics analysis in the model
   % coefficients in E MATRIX show which disturbances are most influtential
   % max(model.dist.d,[],1)  = largest elements in disturbances or maybe use average values in the disturbances vector
if UseParam.disturb       
        influecneD = (sum(abs(model.plant.Ed),1).*max(dist.d,[],1));
        influecneD = influecneD/max(influecneD);  % normalized influecne of each disturbance in model dynamics  
        use_feature_model = influecneD > 1-PrecisionParam.model;   % index of choosen features from model dynamics  
else
        use_feature_model = ones(1,nx); % no features have been discarded
end
        

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
use_feature = use_feature_model + use_feature_PCA + use_feature_lincols > 2; % index of choosen features
Xreduced = X(:,use_feature); % choosen features
Xdiscard = X(:,not(use_feature)); % discarded features

%% PLOTS     
if plotFlag

            figure         
            subplot(2,2,1)
            plot(Xreduced)
            title('Used features profiles')
                 
            subplot(2,2,2)
            plot(Xdiscard)
            title('Discarded features profiles')
            
            %             subplot(2,2,2)
%             bar(explained)  
%             title('Importance of the principal components')
%             xlabel('Principal components')
%             ylabel('Importance [%]')
            % figure      % effect of the each variable on first nu_pca principal components 
            % bar(coeff(:, 1:nu_pca)) 
     
          if UseParam.PCA
            subplot(2,2,3)
            bar(var_weight)
            title('Importance of the original features on principal components')
            xlabel('Features')
            ylabel('Relative importance [0-1]')
          end         

          if UseParam.disturb
            subplot(2,2,4)
            bar(influecneD)
            title('Importance of the disturbances in the model')
            xlabel('Disturbances (Features)')
            ylabel('Relative importance [0-1]')           
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
% imagesc((model.sim.Ed))
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
        % bar(sum(abs(model.sim.Ed),2))
        % % what is the effect of states on states?
        % bar(sum(abs(model.sim.Ad),1))
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