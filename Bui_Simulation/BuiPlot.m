function BuiPlot(outdata,PlotParam)

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
Nsim = outdata.SimParam.run.Nsim; 
Time = (1:Nsim)*outdata.model.plant.Ts/3600/24;  % days


% preview setup
if outdata.ctrl.MPC.use
    N = outdata.ctrl.MPC.N;
    Nrp = outdata.ctrl.MPC.Nrp;
elseif outdata.ctrl.MLagent.use
    N = outdata.ctrl.MLagent.numDelays;
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


%%   ESTIMATOR

% if not(PlotParam.only_zone)
if outdata.estim.use && PlotParam.plotEstim
 
        % outputs - simulated vs estimated
        figure
        subplot(2,1,1)
        title('Simulated vs estimated outputs');
        hold on
        plot(Time, outdata.data.Y,'-' ,'linewidth', 2);
        plot(Time, outdata.data.Ye, '--', 'linewidth', 2);
%         plot(Time, outdata.data.Y-273.15,'-' ,'linewidth', 2);
%         plot(Time, outdata.data.Ye-273.15, '--', 'linewidth', 2);
        axis tight
        grid on
        ylabel('Temp [^{\circ}C]')
%         legend('simulated','estimated')
        xlabel('time [days]')   
          
        % output estimation error
        subplot(2,1,2)
        plot(Time,(outdata.data.Y-outdata.data.Ye),'linewidth', 2)
        title('Output estimation error')
        axis tight
        grid on
        ylabel('Temp [^{\circ}C]')
        xlabel('time [days]')  
   
        % states - simulated vs estimated 
        figure
        subplot(2, 1, 1);
        plot(Time, outdata.data.X(:,1:end-1), 'linewidth', 2);
        title('states: simulated');
        axis tight
        grid on
        xlabel('time [days]')
        
        subplot(2, 1, 2);
        plot(Time, outdata.data.Xe(:,1:end),'linewidth', 2);
        title('states:  estimated');
        axis tight
        grid on
        xlabel('time [days]')
        
        
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
% end
  
%%   CONTROL

if outdata.ctrl.use && PlotParam.plotCtrl
    
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
    
     
%      OUTPUTS
     figure
     subplot(2, 1, 1); 
     title('Indoor temperature','fontsize',24); 
     
     plot(Time, outdata.data.Y-273.15, 'linewidth', 2);
     axis tight
     grid on
     ylabel('Zone temperatures [\circC]','fontsize',22)  
     set(gca,'fontsize',22)
     hold on
%      box on
% %         slight rotation to prevent misplotting
%      ax = gca;
%      ax.XTickLabelRotation=1; 


% TODO IMPLEMENT - comfort bounds
     if outdata.ctrl.use
         %      TODO: unify References for all controllers for plotting
%          Rmin = outdata.data.wb(1,1:end-N);
%          Rmax = outdata.data.wa(1,1:end-N);
         
         Rmin = mean(outdata.data.wb(:,1:end-Nrp),1);
         Rmax = mean(outdata.data.wa(:,1:end-Nrp),1);
         
%          R = outdata.data.R;
  
         stairs(Time, Rmin-273.15, 'k--', 'linewidth', 2);
         stairs(Time, Rmax-273.15, 'k--', 'linewidth', 2);
     end

%      INPUTS
    subplot(2, 1, 2); 
    title('Heating','fontsize',24);   
    hold on
    h = stairs(Time, outdata.data.U');
    set(h, 'linewidth', 2, 'linestyle', '-');
    axis tight
    grid on
    set(gca,'fontsize',22)
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
if PlotParam.plotPrice
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

