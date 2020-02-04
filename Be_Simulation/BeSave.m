function BeSave(outdata,SaveParam)
%  creates table T and saves it as a .csv file
%  rows = simulation steps
%  columns = selected variables 

if SaveParam.save && outdata.ctrl.MPC.use
    T = table;
    if SaveParam.data.states
        T.states = outdata.data.X(:,1:end-1)';
    end
    if SaveParam.data.outputs
        T.outputs = outdata.data.Y';    
    end
    if SaveParam.data.inputs
        T.inputs = outdata.data.U';    
    end
    if SaveParam.data.disturbances
        T.disturbances = outdata.data.D(:,1:end-outdata.ctrl.MPC.N)';   
    end
    if SaveParam.data.references
        T.wa = outdata.data.wa(:,1:end-outdata.ctrl.MPC.N)';   
        T.wb = outdata.data.wb(:,1:end-outdata.ctrl.MPC.N)';  
    end
    if SaveParam.solver.objective
        T.objective = outdata.solver.OBJ';   
    end
    if SaveParam.solver.duals
        T.duals = outdata.solver.DUALS';   
    end
    if SaveParam.solver.primals
        T.primals = outdata.solver.PRIMALS';   
    end
    if SaveParam.solver.PCA_duals
        T.PrincipalDual = outdata.con_info.PrincipalDual;         
    end
    if SaveParam.solver.SolverTime
        T.SolverTime = outdata.solver.SolverTime';   
    end
    if SaveParam.solver.iters
        T.SolverIters = outdata.solver.ITERS';   
    end
    if SaveParam.solver.specifics
        T.SolverINEQLIN = outdata.solver.INEQLIN';   
        T.SolverEQLIN = outdata.solver.EQLIN';  
    end
  
%     save csv
    writetable(T,[SaveParam.path '/Dataset_Order_',outdata.model.Orders.choice,...
        '_Days_',int2str(outdata.SimParam.run.start),'_',int2str(outdata.SimParam.run.end),'_',outdata.model.buildingType,'.csv'],'Delimiter',',')
    
    if SaveParam.data.ActiveSets
        T_AS = table;
        T_AS.ActiveSets = outdata.con_info.ActiveSets';
        writetable(T_AS,[SaveParam.path '/ActiveSets_Order_',outdata.model.Orders.choice,...
        '_Days_',int2str(outdata.SimParam.run.start),'_',int2str(outdata.SimParam.run.end),'_',outdata.model.buildingType,'.csv'],'Delimiter',',')
    end
    
end                    


end