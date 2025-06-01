function options=mpcOpt()
    % Initialize the options structure for the MPC class
    options=struct('solver','quadprog','verbose',0,'ref',true,'Nc',[],'K',[],'uBounds',[0,1],'xBounds',[],'name',[],'P',[],'Ts',0.1);
end