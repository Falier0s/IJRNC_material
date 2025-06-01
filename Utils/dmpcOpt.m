function options=dmpcOpt()
    % Initialize the options structure for the DMPC class
    options=mpcOpt();
    options.rho = 10;
    options.tol = 1e-3;
    options.maxIter = 200;
    options.updtMode = 'standard';
end