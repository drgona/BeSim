function [sys_dExt, rom] = fGenerateSysAndRom(path_ssm, Ts, x0_value, orders)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate the extended discrete state space model and a set of
    % reduce order models of order = orders.
    % Hankel singular decomposition is used for as MOR technic.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    load(path_ssm);
    nx = size(A,1);

    %% Extend state space to include initial conditions
    % x+ = Ax + [B x0] [u 1]'
    % y = Cx + [D Cx0] [u 1]'
    if length(x0_value) == 1
        x0 = x0_value .* ones(nx,1);
    else
        x0 = x0_value;
    end
    BExt = [B A*x0];
    DExt = [D C*x0];
    
    sys_dExt = c2d(ss(A,BExt,C,DExt),Ts);
    
    rom = cell(length(orders),1);
    for i = 1:length(orders)
        rom{i} = reduce(sys_dExt, orders(i));
    end

end
