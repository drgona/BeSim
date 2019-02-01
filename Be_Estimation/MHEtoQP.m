function [H_val,G_val]=MHEtoQP(model, estim, y, u, d, s_k)
% C=[0 1 0 0];
% x_{k+1}=A_kx_k+B_ku_k+b_k
%y_k=C_kx_k+D_i

% linear MHE time invariant model

% model dimensions
nx=model.pred.nx;
nu=model.pred.nu;
ny=model.pred.ny;
% Model matrices
A_k = model.pred.Ad;
B_k = [model.pred.Bd model.pred.Ed];
C_k = model.pred.Cd;
% D_k = model.pred.Dd;
b_k = model.pred.Gd;
d_k = model.pred.Fd;
% MHE parameters
N = estim.MHEParam.N;
Q_k = estim.MHEParam.Qe;
R_k = estim.MHEParam.Re;
P_k = estim.MHEParam.P;

% data
Y_k = y;
U_k = [u d];

% hessian and gradient initialization
H_k=zeros(nx*N);
G_k=zeros(nx*N,1);

lx=[1:nx];
H_k(lx,lx)=P_k+C_k'*R_k*C_k+A_k'*Q_k*A_k;
for r=1:N-2
     H_k(r*nx+lx,r*nx+lx)=Q_k+A_k'*Q_k*A_k+C_k'*R_k*C_k;
end
H_k((N-1)*nx+lx,(N-1)*nx+lx)=Q_k+C_k'*R_k*C_k;

for r=1:N-1
     H_k((r-1)*nx+lx,r*nx+lx)=-A_k'*Q_k;
     H_k(r*nx+lx,(r-1)*nx+lx)=-Q_k*A_k;
end

G_k(lx)=-2*(P_k*s_k+C_k'*R_k*(Y_k(1,:)'-d_k)-A_k'*Q_k*(B_k*U_k(1,:)'+b_k));
for r=1:N-2
     G_k(r*nx+lx)=2*A_k'*Q_k*(B_k*U_k(r+1,:)'+b_k)-2*C_k'*R_k*(Y_k(r+1,:)'-D_k')-2*Q_k*(B_k*U_k(r,:)'+b_k);
end

% H_val=H_k(:);
% G_val=G_k(:);

H_val=reshape(2*H_k',nx*N,nx*N);
G_val=reshape(G_k',nx*N,1);




