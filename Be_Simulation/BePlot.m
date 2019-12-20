function BePlot(outdata,PlotParam)

%% Description
% plot_data.m - general file for plotting the results from run files

if nargin < 2
    PlotParam.plotStates = 0;        % plot states
    PlotParam.plotDist = 0;        % plot disturbances
    PlotParam.plotEstim = 0;        % plot estimation
    PlotParam.plotCtrl = 0;        % plot control
    PlotParam.plotPrice = 0;        % plot price signal
end

% init
% Nsim = outdata.SimParam.run.Nsim; 
Nsim = length(outdata.data.Y);
Time = (1:Nsim)*outdata.model.plant.Ts/3600/24;  % days

font_use = 14;

% preview setup
if outdata.ctrl.MPC.use
    N = outdata.ctrl.MPC.N;
    Nrp = outdata.ctrl.MPC.Nrp;
elseif outdata.ctrl.MLagent.use
    N = outdata.ctrl.MLagent.numDelays;
    Nrp = N;
else
    N = 0;
end


  %% STATES   - TODO: categorize states based on zones
  
% initial sate conditions in K
x_init = outdata.model.plant.Fd(1);
  
if PlotParam.plotStates
    figure
    plot(Time, outdata.data.X(:,1:end-1)+x_init, 'linewidth', 2);
    title('States');
    axis tight
    grid on
    ylabel('Temp [K]')
    xlabel('time [days]')
end

if PlotParam.plotStates3D
    figure
    ribbon(Time', (outdata.data.X(:,1:end-1)+x_init)');
    title('State trajectories 3D');
    axis tight
    shading flat
    grid on
    zlabel('Temp [K]')
    ylabel('time [days]')
    xlabel('states [-]')
end

 %% DISTURBANCES - TODO: categorize disturbances based on magnitudes
if PlotParam.plotDist
    
    DistToPlot = outdata.data.D(:,1:end-N); 
    D = DistToPlot(any(DistToPlot,2),:);  % removing zero rows
    
%     TODO solve this recursively by discarding items from D!!!!
    if sum(max(abs(D),[],2) >1e4) ~= 0  % disturbances in high magnutudes 10 000 and above
        D1 = D(max(abs(D),[],2) >1e4,:);
        D = D(max(abs(D),[],2) <=1e4,:);
        figure
        plot(Time, D1, 'linewidth', 2);
        datetick
        axis tight
        title('Disturbances');
        grid on
        ylabel('[-]')
        xlabel('time [days]')
    end    
    if sum(max(abs(D),[],2) >1e2) ~= 0  % disturbances in magnutudes of hundreads and thousands 100 to 10 000
        D2 = D(max(abs(D),[],2) >1e2,:);
        D = D(max(abs(D),[],2) <=1e2,:);
        figure
        plot(Time, D2, 'linewidth', 2);
        datetick
        axis tight
        title('Disturbances');
        grid on
        ylabel('[-]')
        xlabel('time [days]')
    end
    if sum(max(abs(D),[],2) >1e0) ~= 0 % disturbances in magnutudes of hundreads and ones 1 to 100
        D3 = D(max(abs(D),[],2) >1e0,:);
        D = D(max(abs(D),[],2) <=1e0,:);
        figure
        plot(Time, D3, 'linewidth', 2);
        datetick
        axis tight
        title('Disturbances');
        grid on
        ylabel('[-]')
        xlabel('time [days]')
    end
    if sum(max(abs(D),[],2) <=1e0) ~= 0 % disturbances in magnutudes of fractions of 1
        D4 = D(max(abs(D),[],2) <=1e0,:);
        figure
        plot(Time, D4, 'linewidth', 2);
        datetick
        axis tight
        title('Disturbances');
        grid on
        ylabel('[-]')
        xlabel('time [days]')
    end
end

if PlotParam.plotDist3D
    DistToPlot = outdata.data.D(:,1:end-N); 
    D = DistToPlot(any(DistToPlot,2),:);  % removing zero rows
    
    figure
    ribbon(Time', D');
    title('Disturbance trajectories 3D');
    axis tight
    shading flat
    grid on
    zlabel('Disturbances [K,kW]')
    ylabel('time [days]')
    xlabel('disturbances [-]')
end


%%   ESTIMATOR

% if not(PlotParam.only_zone)
if outdata.estim.use 
 
    if PlotParam.plotEstim
            % outputs - simulated vs estimated
            figure
            subplot(2,1,1)
            title('Simulated vs estimated outputs');
            hold on
            if  strcmp(outdata.model.buildingType,'HollandschHuys')
                   plot(Time, outdata.data.Y-273.15,'-' ,'linewidth', 2);
                   plot(Time, outdata.data.Ye-273.15, '--', 'linewidth', 2);
            else
                   plot(Time, outdata.data.Y,'-' ,'linewidth', 2);
                   plot(Time, outdata.data.Ye, '--', 'linewidth', 2);
            end

            axis tight
            grid on
            ylabel('Temp. [^{\circ}C]')
    %         legend('simulated','estimated')
            xlabel('time [days]')   

            % output estimation error
            subplot(2,1,2)
            plot(Time,(outdata.data.Y-outdata.data.Ye),'linewidth', 2)
            title('Output estimation error')
            axis tight
            grid on
            ylabel('Temp. [^{\circ}C]')
            xlabel('time [days]')  

               % states - simulated vs estimated 
            figure
            subplot(2, 1, 1);
            plot(Time, outdata.data.X(:,1:end-1)+x_init-273.15, 'linewidth', 2);
            title('states: simulated');
            axis tight
            grid on
            xlabel('time [days]')
            ylabel('Temp. [^{\circ}C]')

            subplot(2, 1, 2);
            plot(Time, outdata.data.Xe(:,1:end)+x_init-273.15,'linewidth', 2);
            title('states:  estimated');
            axis tight
            grid on
            xlabel('time [days]')
            ylabel('Temp. [^{\circ}C]')

        if outdata.estim.MHE.use   % MHE variables 
            figure
            subplot(2,1,1)
            plot(Time, outdata.data.We,'linewidth', 2)
            title('We')
            grid on
            subplot(2,1,2)
            plot(Time, outdata.data.Ve,'linewidth', 2)
            title('Ve')
            grid on
    %         subplot(3,1,3)
    %         plot(Time(outdata.estim.MHE.N+1:end), outdata.data.Y(:,outdata.estim.MHE.N+1:end),'linewidth', 2)
    %         hold on
    %         plot(Time(outdata.estim.MHE.N+1:end), outdata.data.Yestim(:,outdata.estim.MHE.N+1:end),'--','linewidth', 2)
    %         title('Yestim')
    %         grid on
        end
    end
    
    if PlotParam.plotEstim3D
        figure
        ribbon(Time', outdata.data.Xe');
        title('State estimates trajectories 3D');
        axis tight
        shading flat
        grid on
        zlabel('State estimates [K]')
        ylabel('time [days]')
        xlabel('states [-]')
    end
    
        
end
% end
  
%%   CONTROL

if outdata.ctrl.use 
    
%     TODO implement objective function plotting
%      if outdata.ctrl.MPC.use  % mpc
%      % --- objective function ---
%         figure
%         stairs(Time, outdata.data.J, 'linewidth', 2);
%         hold on
%         stairs(Time, outdata.data.J_v, '--','linewidth', 2);
%         stairs(Time, outdata.data.J_u, '--','linewidth', 2);
%         stairs(Time, outdata.data.OBJ, '--','linewidth', 2);
%         legend('objective','obj. viol,','obj. inputs', 'yalmip obj.')
%         title('objective value increments');
%         axis tight
%         grid on       
%      end  
    
     
    if PlotParam.plotCtrl
    %      OUTPUTS
         figure
         subplot(2, 1, 1); 
         title('Indoor temperature','fontsize',font_use+2); 
         hold on
         plot(Time, outdata.data.Y-273.15, 'linewidth', 2);
         axis tight
         grid on
         ylabel('Temperatures [\circC]','fontsize',font_use)  
         set(gca,'fontsize',font_use)
         box on
    % %         slight rotation to prevent misplotting
    %      ax = gca;
    %      ax.XTickLabelRotation=1; 

         if outdata.ctrl.use
             Rmin = mean(outdata.data.wb(:,1:end-Nrp),1);
             Rmax = mean(outdata.data.wa(:,1:end-Nrp),1);
    %          R = outdata.data.R;
             stairs(Time, Rmin-273.15, 'k--', 'linewidth', 2);
             stairs(Time, Rmax-273.15, 'k--', 'linewidth', 2);
         end

    %      INPUTS
        subplot(2, 1, 2); 
        title('Heating','fontsize',font_use+2);   
        hold on
        h = stairs(Time, outdata.data.U');
        set(h, 'linewidth', 2, 'linestyle', '-');
        axis tight
        grid on
        ylabel('Heat flows [W]','fontsize',font_use)  
        xlabel('time [days]','fontsize',font_use)  
        set(gca,'fontsize',font_use)
        box on
    end
       
    if outdata.ctrl.MPC.use
        
        if strcmp(outdata.solver.MPC_options.solver,'+quadprog')
           rows = 3;
        else
           rows  = 2;
        end
        
        if PlotParam.plotPrimalDual
            %  MPC  Objective and Duals
            figure
            subplot(rows, 2, 1); 
            h = stairs(Time, outdata.solver.OBJ);
            set(h, 'linewidth', 2, 'linestyle', '-');
            axis tight
            grid on
            ylabel('Q value [-]')  
            xlabel('time [days]')  
            title('MPC objective value');   
            set(gca,'fontsize',font_use)
            box on 

            subplot(rows, 2, 2); 
            h = stairs(Time, outdata.solver.SolverTime);
            set(h, 'linewidth', 2, 'linestyle', '-');
            axis tight
            grid on
            ylabel('Time [s]')  
            xlabel('time [days]')  
            title('Solver Time');  
            set(gca,'fontsize',font_use)
            box on

            subplot(rows, 2, 3); 
            h = stairs(Time, outdata.solver.DUALS');
            set(h, 'linewidth', 2, 'linestyle', '-');
            axis tight
            grid on
            ylabel('Dual variables [-]')  
            xlabel('time [days]')  
            title('MPC dual variables');  
            set(gca,'fontsize',font_use)
            box on

            subplot(rows, 2, 4); 
            h = stairs(Time, outdata.solver.PRIMALS');
            set(h, 'linewidth', 2, 'linestyle', '-');
            axis tight
            grid on
            ylabel('Primal variables [-]')  
            xlabel('time [days]')  
            title('MPC primal variables');  
            set(gca,'fontsize',font_use)
            box on

            if strcmp(outdata.solver.MPC_options.solver,'+quadprog')
                subplot(3, 2, 5); 
                h = stairs(Time, outdata.solver.INEQLIN');
                set(h, 'linewidth', 2, 'linestyle', '-');
                axis tight
                grid on
                ylabel('INEQLIN [-]')  
                xlabel('time [days]')  
                title('quadprog INEQLIN');  
                set(gca,'fontsize',font_use)
                box on

                subplot(3, 2, 6); 
                h = stairs(Time, outdata.solver.EQLIN');
                set(h, 'linewidth', 2, 'linestyle', '-');
                axis tight
                grid on
                ylabel('EQLIN [-]')  
                xlabel('time [days]')  
                title('quadprog EQLIN');  
                set(gca,'fontsize',font_use)
                box on
            end              
        end

        if PlotParam.plotPrimalDual3D 
            figure
%           3D duals
            subplot(2, 1, 1); 
            ribbon(Time', outdata.solver.DUALS');
%             set(h, 'linewidth', 2, 'linestyle', '-');
%             shading flat
            axis tight
            grid on
            xlabel('Dual variables [-]')  
            ylabel('time [days]')  
            title('MPC dual variables 3D');  
            set(gca,'fontsize',font_use)
            shading(gca,'flat')
            box on

            subplot(2, 1, 2); 
            ribbon(Time', outdata.solver.PRIMALS');
            shading flat
%             set(h, 'linewidth', 2, 'linestyle', '-');
            axis tight
            grid on
            xlabel('Primal variables [-]')  
            ylabel('time [days]')  
            title('MPC primal variables 3D');  
            set(gca,'fontsize',font_use)
            box on               
        end  
        
        if PlotParam.plotDualActive
            tol = 0.01;
            DUAL_active = outdata.solver.DUALS>tol;
            DUAL_active_ineq = zeros(size(DUAL_active));
            DUAL_active_eq = zeros(size(DUAL_active));
            DUAL_active_eq(1:size(outdata.solver.EQLIN,1),:) =  (DUAL_active(1:size(outdata.solver.EQLIN,1),:) + 1) == 2;
            DUAL_active_ineq(size(outdata.solver.EQLIN,1)+1:end,:) =  (DUAL_active(size(outdata.solver.EQLIN,1)+1:end,:) +1) ==2;
            
%             DUAL_active_ineq = (DUAL_active + outdata.optimize_con_info.INDEX') ==2;
%             DUAL_active_eq = (DUAL_active + (outdata.optimize_con_info.INDEX ~=1)') ==2;
                
            figure
            spy(DUAL_active_ineq)
            hold on
            spy(DUAL_active_eq,'r')
            title('MPC dual variables activation');  
            ylabel('Dual variables [-]')  
            xlabel('time [steps]') 
            legend('ineq','eq')      
            
 
%       Labels_idx:            
%       'SSM_single' =  1 
%       'y_zone' = 2
%       'u_box' = 3        
%       'nonnegative_slacks' = 4  

            colours = {'r','b','g','m','c','k'};           
            Legend=cell(length(outdata.con_info.Label_types),1);
            figure
            hold on
            for k = 1:length(outdata.con_info.Label_types)
                temp = zeros(size(DUAL_active));
                temp(outdata.con_info.Labels_idx==k,:) = DUAL_active(outdata.con_info.Labels_idx==k,:);
                spy(temp,colours{k})
                Legend{k} = outdata.con_info.Label_types{k};
            end
            legend(Legend)   
            title('MPC dual variables constraints types');  
            ylabel('Dual variables [-]')  
            xlabel('time [steps]') 
            
            figure
            imagesc(outdata.solver.DUALS)
            title('MPC dual variables  heat map');  
            ylabel('Dual variables [-]')  
            xlabel('time [steps]') 
            colorbar
            
%             if duals reduced
            if (outdata.con_info.DiagnoseParam.Reduce.lincols.use || outdata.con_info.DiagnoseParam.Reduce.PCA.use)               
                Duals_reduced = zeros(size(DUAL_active));
                Duals_discarded = zeros(size(DUAL_active));
                Duals_reduced(outdata.con_info.use_duals,:) = outdata.solver.DUALS(outdata.con_info.use_duals,:); % choosen features (duals)
                Duals_discarded(not(outdata.con_info.use_duals),:) = outdata.solver.DUALS(not(outdata.con_info.use_duals),:); % discarded features (duals)                            
                figure
                spy(Duals_reduced>tol)
                hold on
                spy(Duals_discarded>tol,'r')
                title('MPC dual variables dim. reduction');  
                ylabel('Dual variables [-]')  
                xlabel('time [steps]') 
                legend('duals reduced','duals discarded')      
            
                
            end
            
            
        end
    end

    
    %           box on
% %         slight rotation to prevent misplotting
%           ax = gca;
%           ax.XTickLabelRotation=1;  

%% TODO IMPLEMENT - ctrl bounds
%      if outdata.ctrl.use
%          % control constraints - TODO unify for all controllers
%           Umax = outdata.model.pid.umax*ones(1,Nsim+1);
%           Umin = outdata.model.pid.umin*ones(1,Nsim+1);
%           plot(Time, Umax, 'k--', Time, Umin, 'k--', 'linewidth', 2 );
%           ylabel(['u [W]'],'fontsize',22)
%           ylim([-Umin(1,1)-Umax(1,1)/10 , Umax(1,1)+Umax(1,1)/10])    
%           xlabel('time [days]','fontsize',22)             
%      end

end
  
% energy price profile
if PlotParam.plotPrice && outdata.ctrl.use
    figure
%     price profile
    subplot(3, 1, 1); 
    plot(Time, outdata.data.Price, 'linewidth', 2);
    title('Energy Price Profile');
    axis tight
    grid on
    ylabel('Price [Euro/kW]')
    xlabel('time [days]')
%   Cost = price * energy consumed
    subplot(3, 1, 2); 
    stairs(Time, outdata.data.Cost', 'linewidth', 2);
    title('Energy Cost per Heat Flow');
    axis tight
    grid on
    ylabel('Cost [Euro]')
    xlabel('time [days]')
    subplot(3, 1, 3); 
    stairs(Time, sum(outdata.data.Cost), 'linewidth', 2);
    title('Total Energy Cost');
    axis tight
    grid on
    ylabel('Cost [Euro]')
    xlabel('time [days]')
end


end

