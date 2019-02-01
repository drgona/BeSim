% development and testing srcipt


run BeInit.m



return


% sparsity of A matrix
figure
 imagesc(model.plant.Ad)
 colorbar
 title('A')
 
figure
spy(model.plant.Ad>0.001)
title('A')

return
% paralel coordinates chart and box plots - reveal variance in data

figure
subplot(3,1,1)
plot(outdata.data.X)
subplot(3,1,2)
plot(outdata.data.X')
subplot(3,1,3)
boxplot(outdata.data.X')

figure
subplot(3,1,1)
plot(outdata.data.Y)
subplot(3,1,2)
plot(outdata.data.Y')
subplot(3,1,3)
boxplot(outdata.data.Y')


figure
subplot(3,1,1)
plot(outdata.data.U)
subplot(3,1,2)
plot(outdata.data.U')
subplot(3,1,3)
boxplot(outdata.data.U')


% todo: instead of 2 D plot use 3D plot based on A matrix - dynamical
% influence from different states

% TODO: investigate influence of A matrix on variance in box plots of x

return

% surface plots time series data

figure
subplot(2,2,1)
surf(outdata.data.Y')
shading flat
hold on
% surf(outdata.data.wa(:,1:end-10)')
% surf(outdata.data.wb(:,1:end-10)')
title('Y surface')
subplot(2,2,2)
ribbon(outdata.data.Y')
hold on 
% ribbon(outdata.data.wa(:,1:end-10)')
% ribbon(outdata.data.wb(:,1:end-10)')
shading flat
title('Y ribbon')
subplot(2,2,3)
imagesc(outdata.data.Y)
% contourf(outdata.data.Y',30,'LineColor','none')
title('Y heatmap')
subplot(2,2,4)
plot(outdata.data.Y')
title('Y trajectories')
hold on
plot(outdata.data.wa')
plot(outdata.data.wb')

% state profiles
% optimal state trajectories
figure
subplot(2,2,1)
surf(outdata.data.X')
shading flat
title('X surface')
subplot(2,2,2)
ribbon(outdata.data.X')
shading flat
title('X ribbon')
subplot(2,2,3)
imagesc(outdata.data.X)
% contourf(outdata.data.X',30,'LineColor','none')
title('X heatmap')
shading flat
subplot(2,2,4)
plot(outdata.data.X')
title('X trajectories')
% TODO: order states based on physical classification.
% e.g. floors, walls, roofs, zones

% TODO: CA plot elementwise to project contribution of x towards y
% one step backward propagation of y onto x
% or on step influence of x onto y

% do something similar
% one step influence of d onto y via CAE
% one step influence of u onto y via CAB

% control profiles
figure
subplot(2,2,1)
surf(outdata.data.U')
shading flat
title('U surface')
subplot(2,2,2)
ribbon(outdata.data.U')
shading flat
title('U ribbon')
subplot(2,2,3)
imagesc(outdata.data.U)
title('U heatmap')
% contourf(outdata.data.U',30,'LineColor','none')
subplot(2,2,4)
plot(outdata.data.U')
title('U trajectories')

% disturbance profiles
figure
subplot(2,2,1)
surf(outdata.data.D')
shading flat
title('D surface')
subplot(2,2,2)
ribbon(outdata.data.D')
shading flat
title('D ribbon')
subplot(2,2,3)
imagesc(outdata.data.D)
title('D heatmap')
% contourf(outdata.data.D',30,'LineColor','none')
subplot(2,2,4)
plot(outdata.data.D')
title('D trajectories')


% BU effect of controls u onto states x via B 
BU = (outdata.model.pred.Bd*outdata.data.U)';
figure
subplot(2,2,1)
surf(BU)
shading flat
title('BU surface')
subplot(2,2,2)
ribbon(BU)
shading flat
title('BU ribbon')
subplot(2,2,3)
imagesc(BU')
title('BU heat map')
subplot(2,2,4)
plot(BU)
title('BU trajectories')


% backward step influence of x_k on x_k-1
% x_k-1 = A^-1 x_k
% where I need to be in autonomous system to get to optimal trajectories in
% one step
xA = (outdata.model.pred.Ad'*outdata.data.X)';
figure
subplot(2,2,1)
surf(xA)
shading flat
title('xA surface')
subplot(2,2,2)
ribbon(xA)
shading flat
title('xA ribbon')
subplot(2,2,3)
imagesc(xA')
% contourf(outdata.data.X',30,'LineColor','none')
title('xA heatmap')
shading flat
subplot(2,2,4)
plot(xA)
title('xA trajectories')


% backward dX
% dXb =  x_k - x_k-1
% optimal one step state transition differences of autonomous system
dXb = (outdata.data.X-xA')';
figure
subplot(2,2,1)
surf(dXb)
shading flat
title('dXb surface')
subplot(2,2,2)
ribbon(dXb)
shading flat
title('dXb ribbon')
subplot(2,2,3)
imagesc(dXb')
% contourf(outdata.data.X',30,'LineColor','none')
title('dXb heatmap')
shading flat
subplot(2,2,4)
plot(dXb)
title('dXb trajectories')


% ED effect of disturbances d onto states x via E 
% E*d_k 
ED = (outdata.model.pred.Ed*outdata.data.D(:,1:end-ctrl.MPC.N))';
figure
subplot(2,2,1)
surf(ED)
shading flat
title('Ed_k  surface')
subplot(2,2,2)
ribbon(ED)
shading flat
title('Ed_k  ribbon')
subplot(2,2,3)
imagesc(ED')
title('Ed_k  heat map')
% contour((outdata.model.pred.Ed*outdata.data.D+outdata.model.plant.Gd)')
subplot(2,2,4)
plot(ED)
title('Ed_k  trajectories')


% EDG effect of disturbances d onto states x via E compensated with Gd
% E*d_k + G
EDG = (outdata.model.pred.Ed*outdata.data.D(:,1:end-ctrl.MPC.N)+outdata.model.plant.Gd)';
figure
subplot(2,2,1)
surf(EDG)
shading flat
title('Ed_k + G surface')
subplot(2,2,2)
ribbon(EDG)
shading flat
title('Ed_k + G ribbon')
subplot(2,2,3)
imagesc(EDG')
title('Ed_k + G heat map')
% contour((outdata.model.pred.Ed*outdata.data.D+outdata.model.plant.Gd)')
subplot(2,2,4)
plot(EDG)
title('Ed_k + G trajectories')

% TODO: off states e.g. in infrac - influence of EDG! How was G matrix
% computed?

% x_k+1 = Ax_k + Bu_k + Ed_k + G
X_opt =(outdata.model.pred.Ad*outdata.data.X(:,1:end-1))' + BU + EDG;
% -outdata.data.X(:,2:end)'
figure
subplot(2,2,1)
surf(X_opt)
shading flat
title('Ax_k + Bu_k + Ed_k + G surface')
subplot(2,2,2)
ribbon(X_opt)
shading flat
title('Ax_k + Bu_k + Ed_k + G ribbon')
subplot(2,2,3)
imagesc(X_opt')
title('Ax_k + Bu_k + Ed_k + G heat map')
% contour((outdata.model.pred.Ed*outdata.data.D+outdata.model.plant.Gd)')
subplot(2,2,4)
plot(X_opt)
title('Ax_k + Bu_k + Ed_k + G trajectories')


% backward step influence of EDG on x_k-1
% A'(-E*d_k - G)
% where I would need to be in autonomous system to compensate for upcomming
% disturbances toward zero state and no inputs
EDGA = (outdata.model.pred.Ad'*(-EDG)')';
figure
subplot(2,2,1)
surf(EDGA)
shading flat
title('A^T(-Ed_k - G) surface')
subplot(2,2,2)
ribbon(EDGA)
shading flat
title('A^T(-Ed_k - G) ribbon')
subplot(2,2,3)
imagesc(EDGA')
title('A^T(-Ed_k - G) heat map')
% contour((outdata.model.pred.Ed*outdata.data.D+outdata.model.plant.Gd)')
subplot(2,2,4)
plot(EDGA)
title('A^T(-Ed_k - G) trajectories')

% backward step influence of x_k and EDG on x_k-1 
% x_k = Ax_k-1 +  Edk  + G
% x_k - Edk - G = Ax_k-1
% x_k-1 = A'(x_k - Edk - G)
% where I need to be in previous step  x_k-1 of autonomous system 
% affected by disturnances  Edk + G
% to get to optimal trajectories x_k
xA_EDGA = xA(1:end-1,:) + EDGA;
figure
subplot(2,2,1)
surf(xA_EDGA)
shading flat
title('A^T(x_k - Ed_k - G) surface')
subplot(2,2,2)
ribbon(xA_EDGA)
shading flat
title('A^T(x_k - Ed_k - G) ribbon')
subplot(2,2,3)
imagesc(xA_EDGA')
title('A^T(x_k - Ed_k - G) heat map')
% contour((outdata.model.pred.Ed*outdata.data.D+outdata.model.plant.Gd)')
subplot(2,2,4)
plot(xA_EDGA)
title('A^T(x_k - Ed_k - G) trajectories')

% integrated effect of control actions via Ax dynamics
% dXb =  x_k - A'(x_k - Edk - G)
% optimal one step state transition differences of autonomous system
dX_EDGAb = (outdata.data.X(:,1:end-1)-xA_EDGA')';
figure
subplot(2,2,1)
surf(dX_EDGAb)
shading flat
title('x_k - A^T(x_k - Edk - G) surface')
subplot(2,2,2)
ribbon(dX_EDGAb)
shading flat
title('x_k - A^T(x_k - Edk - G) ribbon')
subplot(2,2,3)
imagesc(dX_EDGAb')
title('x_k - A^T(x_k - Edk - G) heat map')
% contour((outdata.model.pred.Ed*outdata.data.D+outdata.model.plant.Gd)')
subplot(2,2,4)
plot(dX_EDGAb)
title('x_k - A^T(x_k - Edk - G) trajectories')


% -BUA backward effect of controls u onto states x via B and A
% BUA = A'(- Buk)
BUA = (outdata.model.pred.Ad'*(-BU)')';
figure
subplot(2,2,1)
surf(BUA)
shading flat
title('A^T(- Buk) surface')
subplot(2,2,2)
ribbon(BUA)
shading flat
title('A^T(- Buk) ribbon')
subplot(2,2,3)
imagesc(BUA')
title('A^T(- Buk) heat map')
subplot(2,2,4)
plot(BUA)
title('A^T(- Buk) trajectories')


% EDGA + BUA
% EDGABUA = A'(-E*d_k - G) + A'(- Buk)
% influence of inputs and disturbances
% where would I need to be with my states to reach zero state trajectories
% w.r.t. given disturbances and inputs
EDGA_BUA = EDGA + BUA;
figure
subplot(2,2,1)
surf(EDGA_BUA)
shading flat
title('A^T(-Ed_k - G) + A^T(- Bu_k) surface')
subplot(2,2,2)
ribbon(EDGA_BUA)
shading flat
title('A^T(-Ed_k - G) + A^T(- Bu_k) ribbon')
subplot(2,2,3)
imagesc(EDGA_BUA')
title('A^T(-Ed_k - G) + A^T(- Bu_k) heat map')
subplot(2,2,4)
plot(EDGA_BUA)
title('A^T(-Ed_k - G) + A^T(- Bu_k) trajectories')


% TODO: issue inputs and disturbances shifted one step??
% x_k-1 = A^-1 x_k + A'(-E*d_k - G) + A'(- Buk)
% where would I need to be with my states to reach optimal trajectory
% w.r.t. given disturbances and inputs
% xA_EDGA_BUA = xA(1:end-1,:)+ EDGA_BUA;
xA_EDGA_BUA = (outdata.model.pred.Ad'*outdata.data.X(:,1:end-1))' + EDGA_BUA;
figure
subplot(2,2,1)
surf(xA_EDGA_BUA)
shading flat
title('A^Tx_{k+1} + A^T(-Ed_k - G) + A^T(- Bu_k) surface')
subplot(2,2,2)
ribbon(xA_EDGA_BUA)
shading flat
title('A^Tx_{k+1} + A^T(-Ed_k - G) + A^T(- Bu_k) ribbon')
subplot(2,2,3)
imagesc(xA_EDGA_BUA')
title('A^Tx_{k+1} + A^T(-Ed_k - G) + A^T(- Bu_k) heat map')
subplot(2,2,4)
plot(xA_EDGA_BUA)
title('A^Tx_{k+1} + A^T(-Ed_k - G) + A^T(- Bu_k) trajectories')


% plot(outdata.data.X(:,1:end-1)')

% reconstructed optimal trajectories???
% X_opt = (outdata.model.pred.Ad'*xA_EDGA_BUA')'; 
% figure
% subplot(2,2,1)
% surf(X_opt)
% shading flat
% title('X_opt surface')
% subplot(2,2,2)
% ribbon(X_opt)
% shading flat
% title('X_opt ribbon')
% subplot(2,2,3)
% imagesc(X_opt')
% title('X_opt heat map')
% subplot(2,2,4)
% plot(X_opt)
% title('X_opt trajectories')


% CX effect of states x onto outputs via C
figure
subplot(2,2,1)
surf((outdata.model.plant.Cd*outdata.data.X)')
shading flat
title('C*X surface')
subplot(2,2,2)
ribbon((outdata.model.plant.Cd*outdata.data.X)')
shading flat
title('C*X  ribbon')
subplot(2,2,3)
imagesc(outdata.model.plant.Cd*outdata.data.X)
title('C*X heat map')
% contour((outdata.model.pred.Ed*outdata.data.D+outdata.model.plant.Gd)')
subplot(2,2,4)
plot((outdata.model.plant.Cd*outdata.data.X)')
title('C*X trajectories')


% references projected onto state profiles 
% what we want to be
% influence of r_k on x_k
% r = Cx_k + F
% x_k = C'*(r - F)
RC = (outdata.data.wb-model.plant.Fd)'*outdata.model.pred.Cd;
% RC = outdata.model.pred.Cd'*(outdata.data.wb-model.plant.Fd);
figure
subplot(2,2,1)
surf(RC)
shading flat
title('RC surface')
subplot(2,2,2)
ribbon(RC)
shading flat
title('RC ribbon')
subplot(2,2,3)
imagesc(RC')
% contourf(outdata.data.X',30,'LineColor','none')
title('RC heatmap')
shading flat
subplot(2,2,4)
plot(RC)
title('RC trajectories')

% backward r influence on previous state: r_k -> x_k -> x_k-1
% where would we need to be in autonomous system to get to RC
% r = Cx_k + F
% x_k = C'*(r - F)   - (RC)
% x_k = Ax_k-1 
% x_k-1 = A'*x_k     - (RCA)
RCA = (outdata.model.pred.Ad'*RC')';
% RCA = RC*outdata.model.pred.Ad;
figure
subplot(2,2,1)
surf(RCA)
shading flat
title('RCA surface')
subplot(2,2,2)
ribbon(RCA)
shading flat
title('RCA ribbon')
subplot(2,2,3)
imagesc(RCA')
% contourf(outdata.data.X',30,'LineColor','none')
title('RCA heatmap')
shading flat
subplot(2,2,4)
plot(RCA)
title('RCA trajectories')

% verified one step transition of references
RC1 = (outdata.model.pred.Ad*RCA')';
figure
subplot(2,2,1)
surf(RC1)
shading flat
title('RC1 surface')
subplot(2,2,2)
ribbon(RC)
shading flat
title('RC1 ribbon')
subplot(2,2,3)
imagesc(RC1')
% contourf(outdata.data.X',30,'LineColor','none')
title('RC1 heatmap')
shading flat
subplot(2,2,4)
plot(RC1)
title('RC1 trajectories')

% difference in RC - RCA
% visualized transfer from walls onto zones
dRC =  RC - RCA;
figure
subplot(2,2,1)
surf(dRC)
shading flat
title('dRC surface')
subplot(2,2,2)
ribbon(dRC)
shading flat
title('dRC ribbon')
subplot(2,2,3)
imagesc(dRC')
% contourf(outdata.data.X',30,'LineColor','none')
title('dRC heatmap')
shading flat
subplot(2,2,4)
plot(dRC)
title('dRC trajectories')


% AX - forward in time state dynamics propagation
AX = (outdata.model.plant.Ad*outdata.data.X)';
figure
subplot(2,2,1)
surf(AX)
shading flat
title('AX surface')
subplot(2,2,2)
ribbon(AX)
shading flat
title('AX ribbon')
subplot(2,2,3)
imagesc(AX')
% contourf(outdata.data.X',30,'LineColor','none')
title('AX heatmap')
shading flat
subplot(2,2,4)
plot(AX)
title('AX trajectories')


% dAX = X_k+1 - AX_k
dAX = (outdata.data.X(:,2:end)-outdata.model.plant.Ad*outdata.data.X(:,1:end-1))';
figure
subplot(2,2,1)
surf(dAX)
shading flat
title('X-AX surface')
subplot(2,2,2)
ribbon(dAX)
shading flat
title('X-AX ribbon')
subplot(2,2,3)
imagesc(dAX')
% contourf(outdata.data.X',30,'LineColor','none')
title('X-AX heatmap')
shading flat
subplot(2,2,4)
plot(dAX)
title('X-AX trajectories')

% BU + ED + G
BUEDG =  BU+EDG;
% BUEDG = (outdata.model.pred.Bd*outdata.data.U+outdata.model.pred.Ed*outdata.data.D(:,1:end-ctrl.MPC.N)+outdata.model.plant.Gd)';
figure
subplot(2,2,1)
surf(BUEDG)
shading flat
title('BU+E*D+G surface')
subplot(2,2,2)
ribbon(BUEDG)
shading flat
title('BU+E*D+G ribbon')
subplot(2,2,3)
imagesc(BUEDG')
title('BU+E*D+G heat map')
% contour((outdata.model.pred.Ed*outdata.data.D+outdata.model.plant.Gd)')
subplot(2,2,4)
plot(BUEDG)
title('BU+E*D+G trajectories')


%  dAX  = BU + ED + G
%  BU = dAX - ED - G
BU_calc = dAX-EDG;
figure
subplot(2,2,1)
surf(BU_calc)
shading flat
title('BU = dAX - ED - G surface')
subplot(2,2,2)
ribbon(BU_calc)
shading flat
title('BU = dAX - ED - G ribbon')
subplot(2,2,3)
imagesc(BU_calc')
title(' BU = dAX - ED - Gheat map')
% contour((outdata.model.pred.Ed*outdata.data.D+outdata.model.plant.Gd)')
subplot(2,2,4)
plot(BU_calc)
title(' BU = dAX - ED - G trajectories')

% % Solving system of linear equations to obtain u
% B*u = BU
% u = B\BU
% u = outdata.model.plant.Bd\BU(1,:)'

% reconstructing u
u = zeros(model.plant.nu,length(BU_calc));
for jj = 1:length(BU_calc)
    u(:,jj) =outdata.model.plant.Bd\BU_calc(jj,:)';
end

% match with MPC u and reconstructed u
figure
plot(u','--','LineWidth',2)
hold on
plot(outdata.data.U')

% difference reconsturcted and computed controls
sum(sum(outdata.data.U-u))


return

%% 4D plots transition matrix 
% https://matlabnewbie.blogspot.com/2013/07/recently-i-read-post-from-dr.html

figure
i = 2
imagesc(outdata.model.pred.Ad.*outdata.data.X(:,i-1)'-outdata.model.pred.Ad.*outdata.data.X(:,i)')
colorbar

figure
imagesc(outdata.model.pred.Ad.*outdata.data.X(:,i)')
colorbar


AX_time = zeros([size(outdata.model.pred.Ad), 100]);
for i = 2:length(outdata.data.X)
    AX_time(:,:,i) = imagesc(outdata.model.pred.Ad.*outdata.data.X(:,i)');
    
end

figure


% todo: compute influence of states on another states - contribution of A
% matrix elements normalized with state values


return

% FIX ISSUE HERE - use estimator to calculate u from BU???
%  or simple LP with bounded u

% u = zeros()
u = zeros(6,length(BU_calc))
for jj = 1:length(BU_calc)
    u(:,jj) = BU_calc(jj,:)'\outdata.model.plant.Bd;
end

figure
plot(u')

% u = sdpvar(6,length(BU_calc));
% Objective = sum(u)
% for jj = 1:length(BU_calc)
%     Constraints = [BU_calc(jj,:)' == outdata.model.plant.Bd*u(:,jj)];
%     Constraints = [Constraints, zeros(6,1) <= u(:,jj) <= 1000*ones(6,1)];
% end
% optimize(Constraints,Objective)


% TODO: solve issue with this constraints satisfaction problem
u = sdpvar(6,1);
UUU = zeros(6,200)

ops = sdpsettings('verbose',0);
for jj = 201:400
    Constraints = [];
    Objective = [];
    
    Objective = sum(u);
    Constraints = [BU_calc(jj,:)' == outdata.model.plant.Bd*u(:,1)];
    Constraints = [Constraints, zeros(6,1) <= u(:,1) <= 1000*ones(6,1)];
    optimize(Constraints,Objective,ops);
    
    jj
    UUU(:,jj) = value(u);
end

plot(UUU')

% min R surface - Y surface
% Y surface = C* (X surface + D surface + U surface)
% system identification:  mapping od U surface onto X surface under known D
% surface and x0 init surface edge

% CHALLENGE:  define optimal surface dAX
% known = x0, reference trajectories for selected x
% idea:  train on past data of optimal control


% TODO: change d and compute u based on surface algebra BU = dAX - ED - G



%% eigenvalues
return


[A_eig_vec, A_eig_val] = eig(outdata.model.pred.Ad);

figure
subplot(3,1,1)
imagesc(outdata.model.pred.Ad)
    colorbar
subplot(3,1,2)
imagesc(A_eig_vec)
    colorbar
subplot(3,1,3)
imagesc(A_eig_val)
    colorbar
    
    
% todo: plot eigenvectors in 3D    
    
return
    %% dynamic plots
    figure(1)
    subplot(2,1,1)
    imagesc(outdata.model.pred.Ad)
    colorbar
%     subplot(3,1,2)
%     imagesc(model.pred.Bd)
%     colorbar
    subplot(2,1,2)
% %     plot(X(:,i)')
%     subplot(2,2,4)
%     plot(uk')   

 figure(2)
 subplot(2,1,1)
 imagesc(outdata.model.pred.Bd)
 subplot(2,1,2)
 
  figure(3)
 subplot(2,1,1)
 imagesc(outdata.model.pred.Ed)
 subplot(2,1,2)
 

for i = 2:length(outdata.data.X)

%%     dynamic plots   
figure(1)
    pause(0.0001)
    subplot(2,1,1)
 % actual heat map of Ax
%     imagesc(outdata.model.pred.Ad.*outdata.data.X(:,i)')
%       caxis([-25 0])
 % actual heat map of A-Ax
%     imagesc(outdata.model.pred.Ad-outdata.model.pred.Ad.*outdata.data.X(:,i)')
%     colorbar
%     caxis([0 15])
    
     % actual heat map of Ax_k-Ax_k+1
    imagesc(outdata.model.pred.Ad.*outdata.data.X(:,i-1)'-outdata.model.pred.Ad.*outdata.data.X(:,i)')
    colorbar
    caxis([-0.1 0.1])
%     caxis([0 15])
%     subplot(3,1,2)
%     imagesc(model.pred.Bd.*uk(:,i)')
%     colorbar
%     caxis([0 12])
 pause(0.0001)
    subplot(2,1,2)
     plot(outdata.data.X(:,1:i)')
%     subplot(2,2,4)
%     plot(uk')


figure(2)
    pause(0.0001)
    subplot(2,1,1)
      % actual heat map of Bu_k-Bu_k+1
%     imagesc(outdata.model.pred.Bd.*outdata.data.U(:,i-1)'-outdata.model.pred.Bd.*outdata.data.U(:,i)')
%     colorbar
%     caxis([-0.1  0.1])
%     
    imagesc(outdata.model.pred.Bd.*outdata.data.U(:,i)')
    colorbar
    caxis([0  1])
    
    pause(0.0001)
    subplot(2,1,2)
    plot(outdata.data.U(:,1:i)')
    
    
     figure(3)
     pause(0.0001)
     subplot(2,1,1)
%      imagesc((outdata.model.pred.Ed.*outdata.data.D(:,i)'))
%      colorbar
%       caxis([-10  1])
    imagesc(outdata.model.pred.Ed.*outdata.data.D(:,i-1)'-outdata.model.pred.Ed.*outdata.data.D(:,i)')
     colorbar
      caxis([-0.05  0.05])
     
     pause(0.0001)
     subplot(2,1,2)
    plot(outdata.data.D(:,1:i)')
    
end

return
% state, disturbances and iput trajectories
figure
    subplot(5,1,1)
 plot((outdata.model.pred.Bd*outdata.data.U)')
     subplot(5,1,2)
 plot((outdata.model.pred.Ed*outdata.data.D(:,1:end-ctrl.MPC.N)+outdata.model.plant.Gd)')
subplot(5,1,3)
 plot((outdata.model.pred.Ad*outdata.data.X)')
subplot(5,1,4)
 plot((outdata.model.pred.Ed*outdata.data.D(:,1:end-ctrl.MPC.N)+outdata.model.pred.Bd*outdata.data.U)')
subplot(5,1,5)
 plot((outdata.model.pred.Ad*outdata.data.X(:,1:end-1)+outdata.model.pred.Ed*outdata.data.D(:,1:end-ctrl.MPC.N)+outdata.model.pred.Bd*outdata.data.U +outdata.model.plant.Gd)')
 
