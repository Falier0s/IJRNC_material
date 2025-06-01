
%    +#%%%%%%%%%%%%%%%.      =#%%%#+        %%%%+              :%%%%:  +#%%%%%%%%%%%%%%%. %%%%%%%%%%%%%%%%#=   +#%%%%%%%%%%%%%%%*:  :*%%%%%%%%%%%%%%%%%: 
%   *@@@@@@@@@@@@@@@@@.     #@@@@@@@#       @@@@*              -@@@@- *@@@@@@@@@@@@@@@@@. @@@@@@@@@@@@@@@@@@- *@@@@@@@@@@@@@@@@@@@  @@@@@@@@@@@@@@@@@@@: 
%   *@@@@-------------     #@@@@#@@@@#      @@@@*              -@@@@- *@@@@-------------  -------------+@@@@= *@@@@----=@@@@@@@@@@  @@@@+--------------  
%   *@@@@%%%%%%%%%%%%%    #@@@@= =@@@@#     @@@@*              -@@@@- *@@@@%%%%%%%%%%%%%  %%%%%%%%%%%%%%@@@@= *@@@@   +@@@@*.*@@@@  @@@@%%%%%%%%%%%%%#-  
%   *@@@@@@@@@@@@@@@@@   *@@@@*...+@@@@*    @@@@*              -@@@@- *@@@@@@@@@@@@@@@@@  @@@@@@@@@@@@@@@@@#. *@@@@ -%@@@%:  *@@@@  +@@@@@@@@@@@@@@@@@@: 
%   *@@@@.............  *@@@@@@@@@@@@@@@+   @@@@*              -@@@@- *@@@@.............  @@@@#....=@@@@@-    *@@@@*@@@@=    *@@@@    ............=@@@@: 
%   *@@@@              *@@@@%%%%%%%%%@@@@-  @@@@@@@@@@@@@@@@@* -@@@@- *@@@@@@@@@@@@@@@@@. @@@@*     .#@@@@+   *@@@@@@@@@@@@@@@@@@@  @@@@@@@@@@@@@@@@@@@: 
%   *@@@@             +@@@@+         +@@@@- +@@@@@@@@@@@@@@@@* -@@@@- :%@@@@@@@@@@@@@@@@. @@@@*       -@@@@%. :%@@@@@@@@@@@@@@@@@+  @@@@@@@@@@@@@@@@@@*  
%   .....             .....           .....   ................  ....    ................  .....         .....   ................    .................  
%==========================================================================
% FILE: yDMPC.m
%==========================================================================
% AUTHOR: Fabio Faliero
% 
% DESCRIPTION: Distributed Linear Quadratic MPC using YALMIP
% This function is intended to build a simple MPC controller using YALMIP
% for both regulation and tracking problems of a distributed system using ADMM. 
% 
% COPYRIGHT (c) 2024 Politecnico di Torino
% All rights reserved.




classdef yDMPC < yMPC
    properties
        neighbors       % Cell containing the neighbors objects
        Ni              % Number of neighbors
        An              % Interaction matrices
        Bn
        rho             % ADMM parameter
        xN              % neighbor states
        xaN
        xrN             
        uN              % neighbor inputs
        uaN
        urN
        x0N             % Initial neighbor states
        zx               % Consensus variables
        zu
        lx               % Dual variables
        lu
        x_bar
        u_bar
        x_barP
        u_barP
        xN_bar
        uN_bar            
        lambdax
        lambdau          
        tol             % Tolerance for convergence
        max_iter        % Maximum number of iterations
        Kc_star
        ADMM_star
    end

    properties (Access=public)
        primalRes       % Primal residual
        dualRes         % Dual residual
        gammaX           % Local copy of the states decision variables
        gammaU           % Local copy of the inputs decision variables
        iter            % Iteration counter
        totalTime
        mu              
        tau             
        eps     
        rRho            % reset ADMM parameter        
        Rho             
        updtMode
        nN               % Dimension of auxiliar state variables
        mN               % Dimension of auxiliar decision variables
        status   
        sConvergence    % Self Convergence flag   
        nConvergence    % Self Convergence flag  
        gConvergence    % Self Convergence flag
        K_cost
        ADMM_cost 
        Xn                  % Neighbours state Domain
        Zn
        Un
        Fn                  % Neighbours constraint matrix
        idx
        resetted
        Kn
        Akn
    end
    
    methods
        function obj = yDMPC(sys, Q, R, N, options)
            % Chiamata al costruttore della superclasse
            obj@yMPC(sys, Q, R, N, options);


            if isfield(options, 'rho')&&~isempty(options.rho)
                obj.Rho = options.rho;
                obj.rRho = options.rho;
            else
                % For penalizing heavely consensus instead of using low rho
                % it has been increased to be comparable with the MPC
                % weights
                obj.Rho = 10;
                obj.rRho = 10;
            end
            if isfield(options, 'updtMode')&&~isempty(options.updtMode)
                obj.updtMode = options.updtMode;
            else
                obj.updtMode = 'standard';
            end
            
            if isfield(options, 'tol')&&~isempty(options.tol)
                obj.tol = options.tol;
            else
                obj.tol = 1e-3;
            end

            if isfield(options, 'maxIter')&&~isempty(options.maxIter)
                obj.max_iter = options.maxIter;
            else
                obj.max_iter = 200;
            end

            obj.mu = 1e3;
            obj.tau = 2;
            obj.eps = 1e-4;

            % Convergence flag must be normally set to true at the beginning 
            % to ensure the initialization of the variables in the method solve
            obj.status = 0;
            obj.sConvergence = true;
            obj.K_cost = 0;
            obj.ADMM_cost = 0;
            obj.resetted = false;
        end
        
        function obj = decompose(obj)

            % In here we want to find the minimum parameter decomposition of the steady states tuple (x_s,u_s) 
            % such that those states can be associated to a given setpoint y_r. i.e. E*[x_s;u_s]=F*y_r


            E = [obj.A-eye(obj.n) obj.B obj.An;obj.C zeros(obj.p,obj.m+obj.nN)];
            F = [zeros(obj.n,obj.p);eye(obj.p)];

            [U,S,V]=svd(E);

            % Now we want to define the matricies M and N such that:
            % [x_s;u_s]=M*\nu and y_r=N*\nu  where \nu is the new set of minimal required parameters
            Vp=V(:,rank(S)+1:end); 
            Up=U(:,rank(S)+1:end);
            if size(U,2)==(obj.p+obj.n)
                G=eye(obj.p);
            elseif size(U,2)<(obj.p+obj.n)
                G=(F'*Up);%_p; % add the orthonormal complement of the matrix to be defined
                G=G(:,rank(S)+1:end);
            end

            if size(U,2)<(obj.m+obj.n)
                obj.Mn=[V*pinv(S)*U'*F*G Vp]; % Requires the orthonormal complement Vp to be defined
                obj.Nn=[G zeros(obj.p,obj.m+obj.n+obj.nN-size(U,2))];
            elseif size(U,2)==(obj.m+obj.n)
                obj.Mn=V*pinv(S)*U'*F*G;
                obj.Nn=G;
            end

            obj.nu=sdpvar(repmat(size(obj.Mn,2),1,2),ones(1,2));
%             obj.svd_flag=true;
        end

        function obj = initialize(obj) %#ok<*MANU>
            % get feasible space info from neighbours
            obj.Xn=Polyhedron();
            obj.Un=Polyhedron();
            for i=1:obj.Ni
                obj.Xn=obj.Xn*obj.neighbors{i}.K.X;
                obj.Un=obj.Un*obj.neighbors{i}.K.U;
            end
            
            obj.Zn=obj.Z*obj.Xn;
            obj.Zs= (obj.sigma*obj.Zn);
            obj=obj.addHConstraints();
            if obj.svd_flag
                obj=obj.decompose();
            end
            if ~obj.svd_flag
                obj.Mn=eye(obj.n+obj.m+obj.nN);
            else
                obj.Mx=[eye(obj.n),zeros(obj.n,obj.m+obj.nN)]*obj.Mn;
            end

            if all(obj.K==0,"all")
                obj=obj.stabGain(obj.A,obj.B,[],mode);
                obj.Kn=zeros(obj.m,obj.nN);
            end
             %obj.lam=.99;
            obj.Ak = obj.A+obj.B*obj.K;
            obj.Akn= obj.An+obj.B*obj.Kn;
            L = [-obj.K eye(obj.m)  -obj.Kn]*obj.Mn;
            Aa=[obj.Ak obj.B*L;zeros(size(obj.Mn,2),obj.n) eye(size(obj.Mn,2))];
            Ba=[obj.Akn;zeros(size(obj.Mn,2),obj.nN)];
             %if nargin>1 && ~isempty(K)
                % obj.Ak=[obj.A+obj.B*obj.K obj.An+obj.B*obj.Kn];           
            % OCP construction and solver definition

            obj=obj.contractSet(Aa,Ba, ...
                                Polyhedron([obj.Z.A*[eye(obj.n), zeros(obj.n,size(obj.Mn,2));obj.K, L];[zeros(size(obj.Zs.A,1),obj.n) obj.Zs.A*obj.Mn]],[obj.Z.b;obj.Zs.b]),...
                                obj.Xn,1);

            obj.partSol=false;
            obj=obj.costBuild();
            obj=obj.sysBuild();

            if strcmp(obj.solver,'gurobi')
                Options = sdpsettings('solver',obj.solver,'verbose',obj.verbose,'savesolveroutput',1,'warmstart',1,'gurobi.FeasRelax',1);
                Options.gurobi.FeasibilityTol=0.01;
                Options.gurobi.FeasRelaxBigM=1e4;
            else
                Options = sdpsettings('solver',obj.solver,'verbose',obj.verbose,'savesolveroutput',1);
            end
            % Creazione del controller Yalmip
            obj=obj.setController(Options);
            for j=1:obj.Ni  % For each Neighbour
                    for k=1:length(obj.neighbors{j}.neighbors) % Check in its entire neighbour list
                        if obj.neighbors{j}.neighbors{k}.K==obj % find current agent handle
                            obj.idx(j)=k;                       % Retrieve the relative index
                        end
                    end
            end
        end

        function obj = getNeighbors(obj, neighbours, An, Bn)
            % Obtains the cell of neighbors pointers and the cell of interaction matrices
            
            obj.neighbors = neighbours;
            
            if obj.Ts
                obj.An = obj.Ts*An; 
                obj.Bn = obj.Ts*Bn;
            else
                obj.An = An;
                obj.Bn = Bn;
            end 

            obj.Ni = length(neighbours);
            obj.nN=size(obj.An,2);
            obj.mN=size(obj.Bn,2);
            obj.nConvergence = true(1,obj.Ni);%%
            obj.gConvergence = true(1,obj.Ni);%%

            if ~isempty(An)    
                obj.xN = sdpvar(repmat(obj.nN,1,obj.N+1),ones(1,obj.N+1));
                if obj.track
                    %obj.xaN=sdpvar(obj.nN,1);
                    if iscell(obj.xr)
                        obj.xaN=sdpvar(repmat(obj.nN,1,obj.N+1),ones(1,obj.N+1));
                        obj.xrN=sdpvar(repmat(obj.nN,1,obj.N+1),ones(1,obj.N+1));
                    else
                        obj.xaN=sdpvar(obj.nN,1);
                        obj.xrN=sdpvar(obj.nN,1);
                    end
                else
                    obj.xaN=zeros(obj.nN,1);
                    obj.xrN=zeros(obj.nN,1);
                end
                % Global state vector
                
            else
                obj.xN=cell(1,obj.N+1);
                obj.xN=cellfun(@(x) zeros(obj.nN,1), obj.xN, 'UniformOutput', false);
                obj.xaN=zeros(obj.nN,1);
                obj.xrN=zeros(obj.nN,1);
            end
        
            obj.xN_bar = sdpvar(obj.n+obj.nN,length(obj.xN)+(length(obj.xa)*iscell(obj.xa)+1*isa(obj.xa,'sdpvar')));%+(length(obj.xr)*iscell(obj.xr)+1*isa(obj.xr,'sdpvar'))); %+N is for xa xaN, +N is for xr xrN 
            %obj.xN_bar = sdpvar(obj.n+obj.nN,length(obj.xN)+(1*isa(obj.xa,'sdpvar'))+(length(obj.xr)*iscell(obj.xr)+1*isa(obj.xr,'sdpvar'))); 
            % Dual state variable
            obj.lambdax = sdpvar(obj.n+obj.nN,length(obj.xN)+(length(obj.xa)*iscell(obj.xa)+1*isa(obj.xa,'sdpvar')));%+(length(obj.xr)*iscell(obj.xr)+1*isa(obj.xr,'sdpvar')));
            %obj.lambdax = sdpvar(obj.n+obj.nN,length(obj.xN)+(1*isa(obj.xa,'sdpvar'))+(length(obj.xr)*iscell(obj.xr)+1*isa(obj.xr,'sdpvar')));


            if ~isempty(Bn)
                obj.uN = sdpvar(repmat(obj.mN,1,obj.N),ones(1,obj.N));
                if obj.track
                    %obj.uaN=sdpvar(obj.mN,1);
                    if iscell(obj.ua)
                        obj.uaN=sdpvar(repmat(obj.mN,1,obj.N),ones(1,obj.N));
                        obj.urN=sdpvar(repmat(obj.mN,1,obj.N),ones(1,obj.N));
                    else
                        obj.uaN=sdpvar(obj.mN,1);
                        obj.urN=sdpvar(obj.mN,1);
                    end
                else
                    obj.uaN=zeros(obj.mN,1);
                    obj.urN=zeros(obj.mN,1);
                end

            else
                obj.uN=cell(1,obj.N+1);
                obj.uN=cellfun(@(x) zeros(obj.mN,1), obj.uN, 'UniformOutput', false);
                obj.uaN=zeros(obj.mN,1);
                obj.urN=zeros(obj.mN,1);
            end
                % Global input vector
                obj.uN_bar = sdpvar(obj.m+obj.mN,length(obj.uN)+(length(obj.ua)*iscell(obj.ua)+1*isa(obj.ua,'sdpvar')));%+(length(obj.ur)*iscell(obj.ur)+1*isa(obj.ur,'sdpvar')));
                %obj.uN_bar = sdpvar(obj.m+obj.mN,length(obj.uN)+(1*isa(obj.ua,'sdpvar'))+(length(obj.ur)*iscell(obj.ur)+1*isa(obj.ur,'sdpvar')));
                
                % Dual input variable
                obj.lambdau = sdpvar(obj.m+obj.mN,length(obj.uN)+(length(obj.ua)*iscell(obj.ua)+1*isa(obj.ua,'sdpvar')));%+(length(obj.ur)*iscell(obj.ur)+1*isa(obj.ur,'sdpvar')));
                %obj.lambdau = sdpvar(obj.m+obj.mN,length(obj.uN)+(1*isa(obj.ua,'sdpvar'))+(length(obj.ur)*iscell(obj.ur)+1*isa(obj.ur,'sdpvar')));

                obj.rho = sdpvar;
                obj.tol = sqrt((obj.n+obj.m+obj.nN+obj.mN)*(obj.N+1+(length(obj.xa)*iscell(obj.xa)+1*isa(obj.xa,'sdpvar'))));%+(length(obj.xr)*iscell(obj.xr)+1*isa(obj.xr,'sdpvar'))))*obj.tol;
                %obj.tol = sqrt((obj.n+obj.m+obj.nN+obj.mN)*(obj.N+1+(1*isa(obj.xa,'sdpvar'))+(length(obj.xr)*iscell(obj.xr)+1*isa(obj.xr,'sdpvar'))))*obj.tol;

        end

        function obj= contractSet(obj,Aa,Ba,X0,H,lam)
            % A is the state matrix
            % B is the input matrix
            % K is the control gain matrix
            % Z is the Polyhedron of state and input tuple
            v=size(Aa,3);

            % Initialize the admissible Set 

            t=0;
            Contractive=false;
            na=size(Aa,2);
            Aa=permute(Aa,[1 3 2])./lam;
            
            if ~isempty(Ba)
                %Non-autonomous -> neighbours influence
                ma=size(Ba,2);
                Ba=permute(Ba,[1 3 2])./lam;
            else
                ma = 0;
                Ba = zeros(na,v,ma);
            end
            %figure
            %hold on
            %X0.projection(1:2).plot();
             %c1=X0.chebyCenter().x;
            while ~Contractive
                t=t+1;
                Xa_1=Polyhedron([X0.A,                                         zeros(size(X0.A,1),ma); ...
                                 reshape(pagemtimes(X0.A,Aa),[],na), reshape(pagemtimes(X0.A,Ba),[],ma); ...
                                 zeros(size(H.A,1),na),                         H.A],...
                                 [X0.b;repmat(X0.b,v,1);H.b]);
                if ma~=0
                    Xa_1p=Xa_1.projection(1:na,'fourier','OSQP').minHRep();
                    %[A_proj, b_proj] = obj.project(Xa_1.A, Xa_1.b,na);
                    %Xa_1p = Polyhedron(A_proj, b_proj).minHRep();  % Final polyhedron
                else
                    Xa_1p=Xa_1.minHRep();
                end
                fprintf("Iteration: %d, Constraints: %d\n", t, size(Xa_1p.A, 1));    

                 %c2=Xa_1p.chebyCenter().x;
                %Xa_1.projection(1:2).plot();
                 %if norm(c1-c2,Inf)<1e-6
                    %if Xa_1p==X0
%                 if norm([Xa_1p.chebyCenter().x-X0.chebyCenter().x;Xa_1p.chebyCenter().r-X0.chebyCenter().r],'inf')<1e-4
                    if  Xa_1p.eq(X0)
                        Contractive=true;
                    end
%                 end

                X0=Xa_1p;
                 %c1=c2;
            end
            obj.Xa=X0.minHRep();
        end
        
        function obj = addHConstraints(obj, varargin)
            obj=addHConstraints@yMPC(obj);
            if nargin < 2 || isempty(varargin)
                % Define neighbours constraints in the form of F*x+G*u+F_n*x_n<=1 
                obj.Fn=obj.Zn.A(:,obj.n+obj.m+1:end);                   obj.Fn=obj.Fn(abs(obj.Zn.b)~=1e12,:)./obj.Zn.b(abs(obj.Zn.b)~=1e12);
                % append zeros to F and G to match the length of Fn knowing that Fn has the same number of row of F and G + the number of rows of the neighbours constraints i.e. size(Xn.A,1)
                 %obj.F=[obj.F;zeros(size(obj.Xn.A,1),size(obj.F,2))];     obj.G=[obj.G;zeros(size(obj.Xn.A,1),size(obj.G,2))];
            end
        end    

        function obj = optimize(obj, x0,ref)
            
            % % If no reference is provided, set the final state and input to zero
            % if nargin<2 || isempty(ref)
            %     obj.x{end}=zeros(obj.n,1);
            %     obj.u{end}=zeros(obj.m,1);
            %     obj.xN{end}=zeros(obj.nN,1);
            %     obj.uN{end}=zeros(obj.mN,1);
            % end % For Regulation Problems
            
            obj = obj.reset(x0);
            % Solve the optimization problem
            
            In = {obj.x0};
            if ~isempty(obj.An)
                In{end+1} = obj.x0N;                          % Position 2
                In{end+1} = obj.zx;                         % Position 3
                In{end+1} = obj.lx;                        % Position 4
            end
            if ~isempty(obj.Bn)
                In{end+1} = obj.zu;                         % Position 2/5
                In{end+1} = obj.lu;                        % Position 3/6
            end

            In{end+1} = ref;
            In{end+1} = obj.Rho;
            
            if ~all([obj.sConvergence,obj.nConvergence],'all')%,obj.gConvergence
             % Solve the optimization problem
                [Sol,~,~,~,~,Data] = obj.controller(In);
            
                obj.u_star = Sol{5}(:,1);
                if ~isempty(obj.An)
                    obj.gammaX = [Sol{4};   %x xa  xr
                                  Sol{6}];  %xN xaN xrN
                    if ~isempty(obj.Bn)
                        obj.gammaU = [Sol{5};   %u_sol
                                      Sol{7}];  %uN_sol
                    else
                        obj.gammaU = [Sol{5};
                                      zeros(obj.Ni,size(Sol{5},2))];
                    end
    
                elseif ~isempty(obj.Bn)
                    if ~isempty(obj.Bn)
                        obj.gammaU = [Sol{5};   %u_sol
                                      Sol{6}];  %uN_sol
                    end
                else
                    obj.gammaX = [Sol{4};
                                  zeros(obj.Ni,size(Sol{4},2))];
                    obj.gammaU = [Sol{5};
                                  zeros(obj.Ni,size(Sol{5},2))];
                end
    
                    
                obj.J_star = Sol{1};
                obj.Kc_star = Sol{2};
                obj.ADMM_star = Sol{3};
                obj.time = Data.solvertime;
            else
                obj.time=0;
            end
        end

        % Receive the the self variables from the neighbours and update the self consensus variables

        function obj = updateConsensus(obj,xJ,uJ) 

            % compute the self shared decision variable average. 
            % As [obj.gammaX(1:obj.n,:);xJ] represents the vercat of the self and the neighbors copies of the self decision variables
            % The average is computed by multiplying the horizontal concatenation of identity matrices of the size of the self decision variables 
            % by the vercat of the self and the neighbors copies of the self decision variables

            obj.x_bar=1/abs(obj.Ni+1)*repmat(eye(obj.n),1,obj.Ni+1)*[obj.gammaX(1:obj.n,:);xJ];
            obj.u_bar=1/abs(obj.Ni+1)*repmat(eye(obj.m),1,obj.Ni+1)*[obj.gammaU(1:obj.m,:);uJ]; 
            obj.time=toc(obj.timer)+obj.time;
        end

        % Receive the consensus variables from the neighbours and update the dual variables

        function obj = updateDual(obj,xN_bar,uN_bar)

            % Define the global copy of the decision variables concatenating the self and the neighbors average decision variables

            obj.zx=([obj.x_bar;xN_bar]);
            obj.zu=([obj.u_bar;uN_bar]);
            
            
            % Update the Dual Variables summing to the old one the weighted error b/w the local copy (solution of the optimization problem) and the global copy (average of the neighbors copies)

            obj.lx=obj.lx+obj.Rho*(obj.gammaX-obj.zx);                          
            obj.lu=obj.lu+obj.Rho*(obj.gammaU-obj.zu);                          % update the dual variables
            obj.primalRes = (~isempty(obj.An).*(vec(obj.gammaX-obj.zx)'*vec(obj.gammaX-obj.zx)))+(~isempty(obj.Bn).*(vec(obj.gammaU-obj.zu)'*vec(obj.gammaU-obj.zu)));
            obj.time=toc(obj.timer)+obj.time;
            
        end
        
        function obj = checkConvergence(obj)
            % Compute the primal residual
            
            primEps=obj.tol+obj.eps*max([sqrt((~isempty(obj.An).*(vec(obj.gammaX)'*vec(obj.gammaX)))+(~isempty(obj.Bn).*(vec(obj.gammaU)'*vec(obj.gammaU)))),...
                                         sqrt((~isempty(obj.An).*(vec(obj.zx)'*vec(obj.zx)))+(~isempty(obj.Bn).*(vec(obj.zu)'*vec(obj.zu))))]);
            dualEps=obj.tol+obj.eps*sqrt((~isempty(obj.An).*(vec(obj.lx)'*vec(obj.lx)))+(~isempty(obj.Bn).*(vec(obj.lu)'*vec(obj.lu))));
            % Residuals norms 
            %obj.primalRes = sqrt((~isempty(obj.An).*(vec(obj.gammaX-obj.zx)'*vec(obj.gammaX-obj.zx)))+(~isempty(obj.Bn).*(vec(obj.gammaU-obj.zu)'*vec(obj.gammaU-obj.zu)))); 
            primalRes=obj.primalRes;

            for i=1:obj.Ni
                primalRes=primalRes+obj.neighbors{i}.K.primalRes;
            end

            primalRes = sqrt(primalRes);
            obj.dualRes = sqrt(obj.Ni+1).*obj.Rho*sqrt((~isempty(obj.An).*(vec(obj.zx-obj.x_barP)'*vec(obj.zx-obj.x_barP)))+(~isempty(obj.Bn).*(vec(obj.zu-obj.u_barP)'*vec(obj.zu-obj.u_barP)))); 
            obj.x_barP=obj.zx;
            obj.u_barP=obj.zu;
            %Below the residual is computed using as the 
            %obj.dualRes = sqrt(obj.Ni+1).*obj.Rho*sqrt((~isempty(obj.An).*(vec(obj.x_bar-obj.x_barP)'*vec(obj.x_bar-obj.x_barP)))+(~isempty(obj.Bn).*(vec(obj.u_bar-obj.u_barP)'*vec(obj.u_bar-obj.u_barP)))); 
            %obj.x_barP=obj.x_bar;
            %obj.u_barP=obj.u_bar;
            dRes=obj.dualRes;
            obj.sConvergence = ((primalRes<=primEps && obj.dualRes<=dualEps)||(obj.iter>=obj.max_iter));

            for j=1:obj.Ni  % For each Neighbour
                obj.neighbors{j}.K.nConvergence(obj.idx(j))=obj.sConvergence;
                obj.neighbors{j}.K.gConvergence(obj.idx(j))=all([obj.sConvergence,obj.nConvergence],'all');
            end

            if ~all([obj.sConvergence,obj.nConvergence],'all')%,obj.gConvergence
                obj = obj.continueSolve();
            else
                obj.iter=obj.iter+1;
            end
                      
        end

        function obj = continueSolve(obj)
             %obj.sConvergence = false;

            obj.iter=obj.iter+1;
            % And update the ADMM parameter
            obj.Rho = obj.updateRho(obj.updtMode);
        end

        function obj = solve(obj, x0, ref)

            switch obj.status 
                case 0
                    obj = obj.optimize(x0,ref);
                case 1
                    obj.timer = tic;
                    [xJ, uJ]= obj.syncronization();
                    obj = obj.updateConsensus(xJ,uJ);
                case 2
                    obj.timer = tic;
                    [zetaX, zetaU] = obj.syncronization();
                    obj = obj.updateDual(zetaX,zetaU);
                case 3
                    obj.timer = tic;
                    obj = obj.checkConvergence();
                    obj.verboseMode();
            end

            if obj.status<3
                obj.status = obj.status+1;
            else
                obj.status = 0;
            end
        end 

        function [message1, message2] = syncronization(obj)
        
            switch obj.status 
                case 1
                    % Retrieve the self variable from the neighbours local copies
                    %xJ=zeros(obj.n,obj.N+1+1*isa(obj.xa,'sdpvar')+length(obj.xr)*iscell(obj.xr)+1*isa(obj.xr,'sdpvar'));
                    %uJ=zeros(obj.m,obj.N+1*isa(obj.ua,'sdpvar'));
                    xJ=zeros(obj.n,obj.N+1+length(obj.xa)*iscell(obj.xa)+1*isa(obj.xa,'sdpvar'));%+length(obj.xr)*iscell(obj.xr)+1*isa(obj.xr,'sdpvar'));
                    uJ=zeros(obj.m,obj.N+length(obj.ua)*iscell(obj.ua)+1*isa(obj.ua,'sdpvar'));
                    for j=1:obj.Ni
                        xJ(1+(j-1)*obj.n:j*obj.n,:) = obj.neighbors{j}.K.gammaX((obj.idx(j)*obj.neighbors{j}.K.n)+1:(obj.idx(j)+1)*obj.neighbors{j}.K.n,:); % retrieve the gamma of the neighbour
                        uJ(1+(j-1)*obj.m:j*obj.m,:) = obj.neighbors{j}.K.gammaU((obj.idx(j)*obj.neighbors{j}.K.m)+1:(obj.idx(j)+1)*obj.neighbors{j}.K.m,:); % retrieve the gamma of the neighbour
                    end
                    message1 = xJ;
                    message2 = uJ;

                case 2
                    xN_b=zeros(obj.nN,obj.N+1+length(obj.xa)*iscell(obj.xa)+1*isa(obj.xa,'sdpvar'));%+length(obj.xr)*iscell(obj.xr)+1*isa(obj.xr,'sdpvar'));
                    uN_b=zeros(obj.mN,obj.N+length(obj.ua)*iscell(obj.ua)+1*isa(obj.ua,'sdpvar'));
                    %xN_b=zeros(obj.nN,obj.N+1+1*isa(obj.xa,'sdpvar')+length(obj.xr)*iscell(obj.xr)+1*isa(obj.xr,'sdpvar'));
                    %uN_b=zeros(obj.mN,obj.N+1*isa(obj.ua,'sdpvar'));
                    for j=1:obj.Ni
                        if ~isempty(obj.An)
                            xN_b(1+(j-1)*obj.neighbors{j}.K.n:j*obj.neighbors{j}.K.n,:) = obj.neighbors{j}.K.x_bar;
                        end
                        if ~isempty(obj.Bn)
                           uN_b(1+(j-1)*obj.neighbors{j}.K.m:j*obj.neighbors{j}.K.m,:) = obj.neighbors{j}.K.u_bar;
                        end
                    end
                    message1 = xN_b;
                    message2 = uN_b;
            end

        end
       
    end

    methods (Access=public)
    
        function obj = costBuild(obj)

            % Copy the cost function from the superclass
            obj = costBuild@yMPC(obj);
            % Save the contribute for control purposes

            obj.K_cost=0.5.*obj.J;
         
            % Add the ADMM cost
            for k=1:obj.N+1
                if ~isempty(obj.An)
                    obj.ADMM_cost = obj.ADMM_cost + obj.lambdax(:,k)'*([obj.x{k};obj.xN{k}]-obj.xN_bar(:,k)) ... % Dual term 
                                                  + obj.rho/2*([obj.x{k};obj.xN{k}]-obj.xN_bar(:,k))'*([obj.x{k};obj.xN{k}]-obj.xN_bar(:,k)); %Regularization
                end
                if ~isempty(obj.Bn)
                    obj.ADMM_cost = obj.ADMM_cost+ obj.lambdau(:,k)'*([obj.u{k};obj.uN{k}]-obj.uN_bar(:,k)) + obj.rho/2*([obj.u{k};obj.uN{k}]-obj.uN_bar(:,k))'*([obj.u{k};obj.uN{k}]-obj.uN_bar(:,k));
                end
            end

            if iscell(obj.xr)
                des_sp = [obj.xr{:};obj.xrN{:}];
            else
                des_sp = [obj.xr;obj.xrN];
            end

            if iscell(obj.xa)
                art_sp = [obj.xa{:};obj.xaN{:}];
            else
                art_sp = [obj.xa;obj.xaN];
            end
            obj.ADMM_cost = obj.ADMM_cost + vec(obj.lambdax(:,end-(length(obj.xa)*iscell(obj.xa)+1*isa(obj.xa,'sdpvar')-1):end))'*vec(art_sp-obj.xN_bar(:,end-(length(obj.xa)*iscell(obj.xa)+1*isa(obj.xa,'sdpvar')-1):end)) + obj.rho/2*vec(art_sp-obj.xN_bar(:,end-(length(obj.xa)*iscell(obj.xa)+1*isa(obj.xa,'sdpvar')-1):end))'*vec(art_sp-obj.xN_bar(:,end-(length(obj.xa)*iscell(obj.xa)+1*isa(obj.xa,'sdpvar')-1):end));
            %obj.ADMM_cost = obj.ADMM_cost + vec(obj.lambdax(:,end-(length(obj.xr)*iscell(obj.xr)+1*isa(obj.xr,'sdpvar')-1):end))'*vec(des_sp-obj.xN_bar(:,end-(length(obj.xr)*iscell(obj.xr)+1*isa(obj.xr,'sdpvar')-1):end)) + obj.rho/2*vec(des_sp-obj.xN_bar(:,end-(length(obj.xr)*iscell(obj.xr)+1*isa(obj.xr,'sdpvar')-1):end))'*vec(des_sp-obj.xN_bar(:,end-(length(obj.xr)*iscell(obj.xr)+1*isa(obj.xr,'sdpvar')-1):end));
            %obj.ADMM_cost = obj.ADMM_cost + vec(obj.lambdax(:,end-(length(obj.xa)*iscell(obj.xa)+1*isa(obj.xa,'sdpvar')+length(obj.xr)*iscell(obj.xr)+1*isa(obj.xr,'sdpvar')-1):end-(length(obj.xr)*iscell(obj.xr)+1*isa(obj.xr,'sdpvar'))))'*vec(art_sp-obj.xN_bar(:,end-(length(obj.xa)*iscell(obj.xa)+1*isa(obj.xa,'sdpvar')+length(obj.xr)*iscell(obj.xr)+1*isa(obj.xr,'sdpvar')-1):end-(length(obj.xr)*iscell(obj.xr)+1*isa(obj.xr,'sdpvar')))) + obj.rho/2*vec(art_sp-obj.xN_bar(:,end-(length(obj.xa)*iscell(obj.xa)+1*isa(obj.xa,'sdpvar')+length(obj.xr)*iscell(obj.xr)+1*isa(obj.xr,'sdpvar')-1):end-(length(obj.xr)*iscell(obj.xr)+1*isa(obj.xr,'sdpvar'))))'*vec(art_sp-obj.xN_bar(:,end-(length(obj.xa)*iscell(obj.xa)+1*isa(obj.xa,'sdpvar')+length(obj.xr)*iscell(obj.xr)+1*isa(obj.xr,'sdpvar')-1):end-(length(obj.xr)*iscell(obj.xr)+1*isa(obj.xr,'sdpvar'))));
        
            obj.J=obj.K_cost+obj.ADMM_cost;
        end
        
        function obj = sysBuild(obj)
            
            % Constraints initialization
           
            obj.constraints = [];
            if obj.relax
                obj.s=sdpvar(size([obj.F; obj.Fc * obj.C], 1), 1);
                obj.J=obj.J+0e0*sum(obj.s);
                obj.constraints = [obj.s>=0];
            else
                obj.s=zeros(size([obj.F; obj.Fc * obj.C], 1), 1);
                obj.constraints=[];
            end
            for k=1:obj.N
                obj.constraints = [obj.constraints, ...
                                  ...[[obj.F;obj.Fc*obj.C],[obj.G;zeros(size(obj.Fc,1),size(obj.G,2))],[obj.Fn;zeros(size(obj.Fc,1),size(obj.Fn,2))]]*[obj.x{k};obj.u{k};obj.xN{k}]<=ones(size([obj.F;obj.Fc*obj.C],1),1), ... % inequality constraints Fx+Gu<=1
                                  [[obj.F;obj.Fc*obj.C],[obj.G;zeros(size(obj.Fc,1),size(obj.G,2))]]*[obj.x{k};obj.u{k}]<=ones(size([obj.F;obj.Fc*obj.C],1),1)+obj.s, ... % inequality constraints Fx+Gu<=1
                                  ...[obj.F;obj.Fc]*obj.x{k}+...
                                  ...[obj.G;zeros(size(obj.Fc,1),size(obj.G,2))]*obj.u{k}+...
                                  ...[obj.Fn;zeros(size(obj.Fc,1),size(obj.Fn,2))]*obj.xN{k}<=ones(size([obj.F;obj.Fc],1),1), ... % inequality constraints Fx+Gu<=1 
                                   obj.x{k+1} == obj.A*obj.x{k} + obj.B*obj.u{k} + obj.An*obj.xN{k}];  % System dynamics
                                  
            end

            if obj.track
                if iscell(obj.xa)
                    for k=1:obj.N
                        obj.constraints = [ obj.constraints,...
                                            obj.xa{k+1}==[obj.A,obj.B,obj.An]*[obj.xa{k};obj.ua{k};obj.xaN{k}],... % Artificial Trajectories
                                            ...[obj.A-eye(obj.n),obj.B]*[obj.xa{k};obj.ua{k}] + obj.An*obj.xaN{k}==zeros(obj.n,1),... 
                                            ...[obj.A obj.B obj.An;obj.C zeros(obj.p,obj.m) zeros(obj.p,obj.nN)]*[obj.xr{k};obj.ur{k};obj.xrN{k}] == [obj.xr{k+1};obj.r(:,k)]
                                            ];
                    end
                    obj.constraints = [ obj.constraints, ismember([obj.x{obj.N+1};obj.xa{end};obj.ua{end};obj.xaN{end}],obj.Xa),...
                                                            [obj.A-eye(obj.n),obj.B,obj.An]*[obj.xa{end};obj.ua{end};obj.xaN{end}] ==zeros(obj.n,1),...
                                                            ...[obj.A-eye(obj.n) obj.B obj.An;obj.C zeros(obj.p,obj.m) zeros(obj.p,obj.nN)]*[obj.xr{end};obj.ur{end};obj.xrN{end}] == [zeros(obj.n,1);obj.r(:,end)]];
                                                            ...obj.C*obj.xr{end}==obj.r(:,end)
                                                            ];
                elseif obj.svd_flag
                    obj.constraints = [ obj.constraints, ismember([obj.x{obj.N+1};obj.nu{1}],obj.Xa),...
                                                            [obj.x{end};obj.u{end};obj.xN{end}]==obj.Mn*obj.nu{1},...
                                                            obj.r==obj.Nn*obj.nu{2},...
                                                            ...[obj.xr;obj.ur;obj.xrN]==obj.Mn*obj.nu{2}
                                                            ];
                else
                    obj.constraints = [ obj.constraints, ismember([obj.x{obj.N+1};obj.x{end};obj.u{end};obj.xN{end}],obj.Xa),...
                                                            [obj.A-eye(obj.n),obj.B]*[obj.x{end};obj.u{end}] + obj.An*obj.xN{end}==zeros(obj.n,1),...
                                                            ...[obj.xr;obj.ur;obj.xrN] ==zeros(obj.n+obj.m+obj.nN,1)];
                                                            ...[obj.A-eye(obj.n) obj.B obj.An;obj.C zeros(obj.p,obj.m) zeros(obj.p,obj.nN)]*[obj.xr;obj.ur;obj.xrN] == [zeros(obj.n,1);obj.r]
                                                            ];
                end

            end
        end
        
        function [Inputs,Outputs] = setArguments(obj)
        
            % Takes as input the update of consensus and dual variables and the initial state of self states and neighbors states
            Inputs = {obj.x{1}};                                    % Position 1
            Outputs = {obj.J,obj.K_cost,obj.ADMM_cost};             % Positions 1  2  3
            if iscell(obj.xa)&&iscell(obj.xr)
                Outputs{end+1} = [obj.x{:}, obj.xa{:}];%, obj.xr{:}];                        % Position 4
                Outputs{end+1} = [obj.u{:}, obj.ua{:}];                        % Position 5
            elseif isa(obj.xr,'sdpvar')&&isa(obj.xa,'sdpvar')
                Outputs{end+1} = [obj.x{:}, obj.xa];%, obj.xr];                        % Position 4
                Outputs{end+1} = [obj.u{:}, obj.ua];                        % Position 5
            end
             %Outputs{end+1} = [obj.u{:}, obj.ua];                        % Position 5
            if ~isempty(obj.An)
                Inputs{end+1} = obj.xN{1};                          % Position 2
                Inputs{end+1} = obj.xN_bar;                         % Position 3
                Inputs{end+1} = obj.lambdax;                        % Position 4
                if iscell(obj.xaN)&&iscell(obj.xrN)
                    Outputs{end+1} = [obj.xN{:}, obj.xaN{:}];%, obj.xrN{:}];                        % Position 6
                elseif isa(obj.xrN,'sdpvar')&&isa(obj.xaN,'sdpvar')
                    Outputs{end+1} = [obj.xN{:}, obj.xaN];%, obj.xrN];                        % Position 6
                end
            end
            if ~isempty(obj.Bn)
                Inputs{end+1} = obj.uN_bar;                         % Position 2/5
                Inputs{end+1} = obj.lambdau;                        % Position 3/6
                Outputs{end+1} = [obj.uN{:}, obj.uaN];                       % Position 6/7
            end
            
            Inputs{end+1} = obj.r;                                  % Position end-1
            Inputs{end+1} = obj.rho;                                % Position end
        
        end

        function obj = setController(obj, Options)
            [Inputs,Outputs] = obj.setArguments;
            obj.controller = optimizer(obj.constraints, obj.J, Options, Inputs, Outputs);        
        end

        function rho = updateRho(obj,mode)
            switch mode
                   
                case 'bounded'
                    maxRho = 1e3;  % Valori massimi e minimi per rho
                    minRho = 1e-3;
                    
                    if obj.primalRes > obj.mu * obj.dualRes
                        rho = min(obj.tau * obj.Rho, maxRho);  % Limita l'incremento
                    elseif obj.dualRes > obj.mu * obj.primalRes
                        rho = max(obj.Rho / obj.tau, minRho);  % Limita la diminuzione
                    else
                        rho = obj.Rho;
                    end
                case 'gradual'
                    ratio = obj.primalRes / obj.dualRes;
                    if ratio > obj.mu
                        rho = obj.Rho * min(obj.tau * ratio, 1.5);  % Aumenta gradualmente fino a un massimo di 1.5
                    elseif 1/ratio > obj.mu
                        rho = obj.Rho / min(obj.tau / ratio, 1.5);  % Diminuisci gradualmente fino a un massimo di 1/1.5
                    else
                        rho = obj.Rho;
                    end
                case 'frequency'
                    updateFrequency = 5;  % Aggiorna rho ogni 10 iterazioni
                    residualThreshold = 1e-3;  % Soglia di variazione dei residui 
                    if mod(obj.iter, updateFrequency) == 0 || ...
                            abs(obj.primalRes - obj.dualRes) > residualThreshold
                        if obj.primalRes > obj.mu * obj.dualRes
                            rho = obj.tau * obj.Rho;
                        elseif obj.dualRes > obj.mu * obj.primalRes
                            rho = obj.Rho / obj.tau;
                        else
                            rho = obj.Rho;
                        end
                    else
                        rho = obj.Rho;  % Mantieni rho invariato
                    end
                case 'fixed'
                    rho = obj.Rho;
                otherwise
                    if obj.primalRes>obj.mu*obj.dualRes
                        rho=obj.tau*obj.Rho;
                    elseif obj.dualRes>obj.mu*obj.primalRes
                        rho=obj.Rho/obj.tau;
                    else
                        rho=obj.Rho;
                    end
            end



        end

        function obj = reset(obj,x0)
            % If the whole network has converged or the maximum number of iterations has been reached
            flag=all(obj.gConvergence,'all');
            for j=1:obj.Ni
                flag=flag&&all(obj.neighbors{:,j}.K.gConvergence,'all');
            end
              if ~obj.resetted %|| flag
             %if all([obj.gConvergence,obj.nConvergence,obj.sConvergence])
                obj.resetted=true;
                obj.x0=x0;                                      % Update the initial state
                obj.x0N=[];                                     % Reset the neighbor states     
                for j=1:obj.Ni
                    obj.x0N=[obj.x0N;obj.neighbors{j}.x0];      % Update the neighbor states
                end

                if ~isempty(obj.x_barP)
                    obj.x_barP(:,1:obj.N)=obj.x_barP(:,2:obj.N+1);             % Warmstart the consensus variables about states
                    last=obj.N+1;
                    if iscell(obj.xa)
                        obj.x_barP(:,(last+1):(obj.N+last))=obj.x_barP(:,(last+2):(last+obj.N+1)); % Warmstart the consensus variables about artificial setpoints
%                         last=last+obj.N;
                    end
%                     last=last+1;
%                     if iscell(obj.xr)
%                         obj.x_barP(:,(last+1):(obj.N+last))=obj.x_barP(:,(last+2):(last+obj.N+1)); % Warmstart the consensus variables about artificial setpoints
%                     end
                else
                    obj.x_barP=zeros(obj.n+obj.nN,obj.N+1+length(obj.xa)*iscell(obj.xa)+1*isa(obj.xa,'sdpvar'));%+length(obj.xr)*iscell(obj.xr)+1*isa(obj.xr,'sdpvar'));                    % Reset the consensus variables
                end
                if ~isempty(obj.u_barP)
                    obj.u_barP(:,1:obj.N-1)=obj.u_barP(:,2:obj.N);             % Warmstart the consensus variables
                else
                    obj.u_barP=zeros(obj.m+obj.mN,obj.N+length(obj.ua)*iscell(obj.ua)+1*isa(obj.ua,'sdpvar'));                    % Reset the consensus variables
                end
                if ~isempty(obj.lx)
                    obj.lx(:,1:obj.N)=obj.lx(:,2:obj.N+1);          % Warmstart the dual variables
                    last=obj.N+1;
                    if iscell(obj.xa)
                        obj.lx(:,(last+1):(obj.N+last))=obj.lx(:,(last+2):(last+obj.N+1)); % Warmstart the consensus variables about artificial setpoints
%                         last=last+obj.N;
                    end
                    last=last+1;
%                     if iscell(obj.xr)
%                         obj.lx(:,(last+1):(obj.N+last))=obj.lx(:,(last+2):(last+obj.N+1)); % Warmstart the consensus variables about artificial setpoints
%                     end
                else
                    obj.lx=zeros(obj.n+obj.nN,obj.N+1+length(obj.xa)*iscell(obj.xa)+1*isa(obj.xa,'sdpvar'));%+length(obj.xr)*iscell(obj.xr)+1*isa(obj.xr,'sdpvar'));             % Reset the dual variables
                end
                if ~isempty(obj.lu)
                    obj.lu(:,1:obj.N-1)=obj.lu(:,2:obj.N);          % Warmstart the dual variables
                else
                    obj.lu=zeros(obj.m+obj.mN,obj.N+length(obj.ua)*iscell(obj.ua)+1*isa(obj.ua,'sdpvar'));             % Reset the dual variables
                end
                if ~isempty(obj.zx)
                    obj.zx(:,1:obj.N)=obj.zx(:,2:obj.N+1);          % Warmstart the global variables
                    last=obj.N+1;
                    if iscell(obj.xa)
                        obj.zx(:,(last+1):(obj.N+last))=obj.zx(:,(last+2):(last+obj.N+1)); % Warmstart the consensus variables about artificial setpoints
%                         last=last+obj.N;
                    end
%                     last=last+1;
%                     if iscell(obj.xr)
%                         obj.zx(:,(last+1):(obj.N+last))=obj.zx(:,(last+2):(last+obj.N+1)); % Warmstart the consensus variables about artificial setpoints
%                     end
                else
                    obj.zx=zeros(obj.n+obj.nN,obj.N+1+length(obj.xa)*iscell(obj.xa)+1*isa(obj.xa,'sdpvar'));%+length(obj.xr)*iscell(obj.xr)+1*isa(obj.xr,'sdpvar'));             % Reset the global variables 
                end
                if ~isempty(obj.zu)
                    obj.zu(:,1:obj.N-1)=obj.zu(:,2:obj.N);          % Warmstart the global variables
                else
                    obj.zu=zeros(obj.m+obj.mN,obj.N+length(obj.ua)*iscell(obj.ua)+1*isa(obj.ua,'sdpvar'));             % Reset the global variables
                end
                obj.iter=0;                                         % Reset the iteration counter
                obj.sConvergence=false;                             % Reset the convergence flags
                obj.nConvergence = false(1,obj.Ni);
                obj.gConvergence = false(1,obj.Ni);
                obj.Rho = obj.rRho;
                obj.u_star=[];                                      % Reset the ADMM parameter  
                obj.totalTime = 0;
            end 

            obj.time=0;
        end

        function obj = updateTime(obj)
            obj.time=toc(obj.timer)+obj.time;
            obj.totalTime = obj.totalTime + obj.time;
        end
        
        function saveToFile(obj)
            fname = obj.name;
            fid = fopen([fname '.csv'], 'a'); % Open the file
            
            if fid == -1    % If the file does not exist create it
                error('Cannot open file: %s', fname);
            end
            
            fseek(fid, 0, 'eof');
            file_size = ftell(fid);
            if file_size == 0
                fprintf(fid,'Iteration, IterTime(s), ElapsTime(s), (Sub)OpInput, (Sub)OpCost, kCost, consCost, PrimalRes, DualRes, Conv, rho\n'); % Write the header
            end
            
            fprintf(fid,'%d, %.5f, %.5f, %s, %.3f, %.3f, %.3f, %.4f, %.4f, %d, %.2f\n',obj.iter,obj.time,obj.totalTime,mat2str(obj.u_star),obj.J_star,double(obj.K_cost),double(obj.ADMM_cost),obj.primalRes,obj.dualRes,all(obj.gConvergence,'all'),obj.Rho); % Write the data
            fclose(fid); % Close the file
        end
        
        function printVerbose(obj)
            s=sprintf(' [%s] |\t Iter: %d |\t Solver: %s  |\t Iter Time: %f |\t Elapsed Time: %f |\t (Sub) Optimal Input: %s |\t (Sub) Optimal Cost: %f |\t Primal Res: %f |\t Dual Res: %f |\t Converged: %d\n', obj.name,obj.iter, obj.solver, obj.time, obj.totalTime, mat2str(obj.u_star),obj.J_star,obj.primalRes,obj.dualRes,all(obj.gConvergence,'all'));
            disp(s) %#ok<DSPSP>
        end
    
    end
end
