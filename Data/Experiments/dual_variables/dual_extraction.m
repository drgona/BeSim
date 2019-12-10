

X = sdpvar(2,1);

F = [sum(X) == 1];
F = F+ [ [-2;-2] <=  X <= [2;2] ];
obj = X'*X

% optimize(F,obj);
% double(X)
% dual(F(1))

options = sdpsettings('verbose', 1, 'solver','quadprog', 'savesolveroutput', 1, 'saveduals', 1 );
QP_optim = optimizer(F, obj, options, [], X );


return

sdpvar a x
Constraints = [a+1 <= x];
Objective = x^2;
P = optimizer(Constraints,Objective,[],a,x)
P(1)
z = (-5:0.1:5);
plot(z,P(z))

return

    % optimizer for dynamic comfort zone
%     if MPCParam.dualize
% %         OutputParams = cell(length(con)+2,1);
% %         OutputParams{1} = u(:,1);
% %         OutputParams{2} = obj;
% %     
% %         for i = 3:size(con)
% %            OutputParams{i} = dual(con(i))      
% %         end      
%         mpc.con = con;
%         mpc.obj =  obj;
%         mpc.options = options;
%         mpc.params = { x(:, 1), d_prev, wa_prev, wb_prev, price };
%         mpc.out = {u(:,1); obj







