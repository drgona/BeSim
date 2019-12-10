clear all;
close all;
define_constants;
casedetails = 'case6ww';
mpc = loadcase(casedetails);
Y = makeYbus(mpc);
D = makeIncidence(mpc.bus,mpc.branch);
nbus = size(Y,1);
%results = runopf(mpc);
B = real(Y);
slack_ref = 1;

nbus = size(B,1);
pcap = zeros(nbus,1);
pcap(1) = 2;
pcap(5) = 1;
rcap = zeros(nbus,1);
rcap(1) = 1;
rcap(5) = 0;

T = 24;
d = (2/nbus)*ones(nbus,T);

x_init = rand(nbus,1);

%Problem definition
theta = sdpvar(nbus,T);
p = sdpvar(nbus,T);
dflex = sdpvar(nbus,T);
r = sdpvar(nbus,T);
r_ext = sdpvar(1,T);
x = sdpvar(nbus,T+1);
cons =[];
cons = [theta(slack_ref,1:T) == 1];
cons = [cons B*theta == p - (d + dflex)];
cons = [cons -0.1 <= dflex <= 0.1];
cons = [cons sum(r) + r_ext >= 1];
cons = [cons r_ext >= 0];
cons = [cons p+r <= repmat(pcap,1,T)];
cons = [cons p+r >= 0];
cons = [cons r >= 0];
cons = [cons r <= repmat(rcap,1,T)];
cons = [cons p >= 0];
 cons = [cons x >= 0; x <= 1];
 cons = [cons x(:,1) == x_init];
 cons = [cons x(:,2:end) == x(:,1:end-1) + dflex];

nbranch = size(mpc.branch,1);
for i = 1:nbranch
    f = mpc.branch(i,F_BUS);
    t = mpc.branch(i,T_BUS);
    cons = [cons B(f,t)*(theta(f,:) - theta(t,:)) >= -0.4 B(f,t)*(theta(f,:) - theta(t,:)) <= 0.4];    
end
%cons = [cons ones(1,nbus)*(p-d) == 0];
%cons = [cons H*(d-p) <= 0.5];
%cons = [cons H*(d-p) >= -0.5];
c = 0.1*1:nbus;
lambda_ext = mean(c);

obj = 0;
for i = 1:T
    obj = obj + p(:,i)'*diag(c)*p(:,i) + 1*(x(:,i) - 0.5*ones(nbus,1))'*(x(:,i) - 0.5*ones(nbus,1));
    %obj = obj + p(:,i)' * diag(c) * p(:,i);
end
obj = obj + 10*(x(:,T+1) - 0.5*ones(nbus,1))'*(x(:,T+1) - 0.5*ones(nbus,1)) + sum(lambda_ext*r_ext);
ops=sdpsettings('solver','mosek','verbose',1);
sol=optimize(cons,obj,ops);

gamma = dual(cons(2));
lambda = dual(cons(4)); 
u = dual(cons(1));

reserve_p = lambda;
energy_price = gamma;
