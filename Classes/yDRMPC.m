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
% FILE: yDRMPC.m
%==========================================================================
% AUTHOR: Fabio Faliero
% 
% DESCRIPTION: Distributed Linear Quadratic MPC using YALMIP
% This function is intended to build a simple MPC controller using YALMIP
% for both regulation and tracking problems of a distributed system using ADMM. 
% 
% COPYRIGHT (c) 2024 Politecnico di Torino
% All rights reserved.

classdef yDRMPC < yDMPC
    properties
        W                   % Disturbance set
        c                   % Real input
        th                  % uncertain state parameters estimate
        thC                 % uncertain output parameters estimate
        Th                  % uncertain state parameters set
        ThC                 % uncertain output parameters set 
        q                   % number of state parameters
        qC                  % number of output parameters
        alpha               % tube scaling
        V                   % tube Shape
        Av                  % state matrix evaluated at the vertices
        Bv                  % input matrix evaluated at the vertices
        Cv                  % output matrix evaluated at the vertices
        Ak_v                % closed loop system evaluated at the vertices of the uncertain set
        Akn_v
        w_bar               % disturbance maximum value
        Hbar                
        HCbar   
        th_v                % vertices of the uncertain state parameter set (transposed)
        v                   % number of vertices of the uncertain state parameter set
        thC_v               % vertices of the uncertain output parameter set (transposed)
        vC                  % number of vertices of the uncertain output parameter set
        na
        mask  

    end

    properties (Access=public)
    
        Anv                 % Neighbours state matrix evaluated at verteces
        

    end


    
    methods
        function obj = yDRMPC(sys, Q, R, N, options)
            % Chiamata al costruttore della superclasse
            obj@yDMPC(sys, Q, R, N, options);

            % Handling uncertain system , noice set and 
            
            % An must be taken using the same mask of the A matrix such that depends only on common parameters
            C_mask = squeeze(any(sys.C ~= 0, [1, 2]));
            obj.mask= squeeze(any([any(sys.A ~= 0, [1, 2]),any(sys.B ~= 0, [1, 2])],2));

            obj.A=obj.A(:,:,obj.mask);
            obj.B=obj.B(:,:,obj.mask);
            obj.C=obj.C(:,:,C_mask);
            
            
            obj.q=sum(obj.mask)-1;
            obj.qC=sum(C_mask)-1;

            if (isfield(options,'W')||~isempty(options.W))&&isa('options.W','Polyhedron')
                % if the disturbance set is provided integrate it in the controller
                obj.W = options.W;
            else
                % else set the disturbance set as a small ball around the origin
                obj.W = Polyhedron('lb',-ones(obj.n,1).*1e-1,'ub',ones(obj.n,1).*1e-1);
            end

            obj.c=sdpvar(repmat(obj.m,1,obj.N),ones(1,obj.N));

        end

        function obj = initUnc(obj,Th,ThC)
        
            if ~isempty(Th)
                   obj.Th = Th;        % State Parameter Set
                   [obj.th_v, obj.v] = obj.setParamSet(Th);    
            else
               obj.Th = Polyhedron('V',eye(obj.q));
               if obj.q>0
                   [obj.th_v, obj.v] = obj.setParamSet(Th);
               else
                   obj.th_v = 1;
                   obj.v = 1;
               end
            end
            obj.th=sdpvar(obj.q+1,1);
            obj.Av=permute(pagemtimes(permute(obj.A,[1 3 2]),obj.th_v'),[1 3 2]);
            obj.Bv=permute(pagemtimes(permute(obj.B,[1 3 2]),obj.th_v'),[1 3 2]);
            obj.Anv=permute(pagemtimes(permute(obj.An,[1 3 2]),obj.th_v'),[1 3 2]);
            obj.Xn=Polyhedron();
            obj.Un=Polyhedron();
            for i=1:obj.Ni
                obj.Xn=obj.Xn*obj.neighbors{i}.K.X;
                obj.Un=obj.Un*obj.neighbors{i}.K.U;
            end
            obj.Zn=obj.Z*obj.Xn;
            obj.Zs= (obj.sigma*obj.Zn);
            obj=obj.addHConstraints();
            if nargin>2 && ~isempty(ThC)
                obj.ThC = ThC;       % Output Parameter Set
                [obj.thC_v, obj.vC] = obj.setParamSet(ThC);
            else
                obj.ThC = Polyhedron('V',eye(obj.qC));
                if obj.qC>0
                    [obj.thC_v, obj.vC] = obj.setParamSet(ThC);
                else
                    obj.thC_v = 1;
                    obj.vC = 1;
                end
            end
            obj.thC=sdpvar(obj.qC+1,1);
            obj.Cv=obj.paramEval(obj.C,obj.thC_v); % Output matrix evaluated at the vertices
        end

        function obj = initialize(obj,mode)
        
            % In the initialization we do all the OFFLINE computations
            % and we build the OCP problem to be solved online using YALMIP

            % Offline computations:
            % 1. Compute the tube gain K_tube and the terminal weight P    **OFFLINE:STEP1**
            % 2. Compute the tube shape V as a contractive set
            % 3. Compute the bounding worst case disturbace (i.e. w_bar),
            %    tube evolution condition (i.e. matrix H) 
            %    and the constraint tightening (i.e. matrix Hc)
            if all(obj.K==0,"all")
                obj=obj.stabGain(obj.Av,obj.Bv,obj.W.V,mode); 
                obj.Kn=zeros(obj.m,obj.nN);
            end
            obj.Ak = obj.A+pagemtimes(obj.B,obj.K);
            obj.Akn= obj.An+pagemtimes(obj.B,obj.Kn);

            % obj.Ak=[obj.A+pagemtimes(obj.B,obj.K),obj.An+pagemtimes(obj.B,obj.Kn)];
    
            obj.Ak_v=permute(pagemtimes(permute(obj.Ak,[1 3 2]),squeeze(obj.th_v)'),[1 3 2]);                                    % Define closed loop system evaluated at the vertices of the uncertain set
            obj.Akn_v=permute(pagemtimes(permute(obj.Akn,[1 3 2]),squeeze(obj.th_v)'),[1 3 2]);   
            
            % 2. Compute the tube shape V as a contractive set
            obj=obj.contractSet(obj.Ak_v,obj.Akn_v,Polyhedron(obj.Z.A*[eye(obj.n);obj.K],obj.Z.b),obj.Xn,obj.lam);                                              % 2. Compute the tube shape V as a contractive set
            obj.V=obj.Xa;
            obj.na=length(obj.V.A);
            obj.alpha=sdpvar(repmat(obj.na,1,obj.N+1),ones(1,obj.N+1));
            obj=obj.tubeInclusion();                                                    % 3. Compute the bounding worst case disturbace (i.e. w_bar), tube evolution condition (i.e. matrix Hbar) 
            obj=obj.tightConstraints();                                                 % and the constraint tightening (i.e. matrix HCbar)
            % OCP construction and solver definition
            for j=1:obj.Ni  % For each Neighbour
                for k=1:length(obj.neighbors{j}.neighbors)
                    if obj.neighbors{j}.neighbors{k}.K==obj
                        obj.idx(j)=k;
                    end
                end
            end
            obj=obj.costBuild();
            obj=obj.sysBuild();            
            Options = sdpsettings('solver',obj.solver,'verbose',obj.verbose,'savesolveroutput',1);
            % Creazione del controller Yalmip
            obj=obj.setController(Options);
            
            
        end

        function obj = optimize(obj, x0,ref)
            
            % if nargin<2 || isempty(ref)
            %     obj.x{end}=zeros(obj.n,1);
            %     obj.u{end}=zeros(obj.m,1);
            %     obj.xN{end}=zeros(obj.n*obj.Ni,1);
            %     obj.uN{end}=zeros(obj.m*obj.Ni,1);
            % end
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

            %sol_r=[obj.A(:,:,1)-eye(obj.n) obj.B(:,:,1) obj.An(:,:,1) obj.Bn(:,:,1);obj.C(:,:,1) zeros(obj.p,obj.m+obj.nN+obj.mN)]\[zeros(obj.n,1);ref];
            %In{end+1} = sol_r(1:obj.n,1);
            
            In{end+1} = ref;
            In{end+1} = obj.Rho;
             %In{end+1} = [1;obj.Th.chebyCenter().x];
             %In{end+1} = [1;obj.ThC.chebyCenter().x];
            In{end+1} = [1;zeros(obj.Th.Dim,1)];
            In{end+1} = [1;zeros(obj.ThC.Dim,1)];
            % Solve the optimization problem
            % Optimization problem depends also on th and thC estimates
            if ~all([obj.sConvergence,obj.nConvergence],'all')
                [Sol,~,~,~,~,Data] = obj.controller(In); % Solve the optimization problem
                obj.u_star = Sol{5}(:,1);
                if ~isempty(obj.An)
                    obj.gammaX = [Sol{4};   %x_sol
                                  Sol{6}];  %xN_sol 
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

    end

    methods (Access=public)

        function obj = sysBuild(obj)
            
            % Constraints initialization
            Ak=obj.paramEval(obj.Ak,obj.th');
            A=obj.paramEval(obj.A,obj.th');
            An=obj.paramEval(obj.An,obj.th');
            B=obj.paramEval(obj.B,obj.th');
            H=permute(pagemtimes(permute(obj.Hbar,[1 3 2]),squeeze(obj.th_v)'),[1 3 2]);
            HC=permute(pagemtimes(permute(obj.HCbar,[1 3 2]),squeeze(obj.thC_v)'),[1 3 2]);
            Gb=vertcat(obj.G,zeros(size(obj.Fc,1),size(obj.G,2)));
            obj.constraints = (obj.V.A*obj.x{1}<=obj.alpha{1}); % Initial Tube
            Akn=obj.paramEval(obj.Akn,obj.th');
            %Fnb=vertcat(obj.Fn,zeros(size(obj.Fc,1),size(obj.Fn,2)));
            %% An must be defined evaluated at the estimated th (the paramenter related to the neighbors)
            %% 
            if iscell(obj.xa) && obj.track

                for k=1:obj.N
                    obj.constraints = [obj.constraints, ...
                                    obj.x{k+1} == Ak*(obj.x{k}-obj.xa{k}) + B*obj.c{k}+ Akn*(obj.xN{k}-obj.xaN{k}) + obj.xa{k+1},...
                                    ismember([obj.xa{k};obj.ua{k};obj.xaN{k}],obj.Zs),...
                                    ...[A-eye(obj.n),B,An]*[obj.xa{k};obj.ua{k};obj.xaN{k}] ==zeros(obj.n,1),... %Artificial equilibria
                                    obj.xa{k+1}==[A,B,An]*[obj.xa{k};obj.ua{k};obj.xaN{k}],... %Artificial Trajectories
                                    ...[A B An;obj.C zeros(obj.p,obj.m) zeros(obj.p,obj.nN)]*[obj.xr{k};obj.ur{k};obj.xrN{k}] == [obj.xr{k+1};obj.r(:,k)],...
                                    obj.u{k} == obj.K*(obj.x{k}-obj.xa{k})+obj.Kn*(obj.xN{k}-obj.xaN{k})+obj.c{k}+obj.ua{k}];  % Control Policy
                
                                    %% Tube evolution must depend on the neighbors states
                    for j=1:obj.v
                    obj.constraints = [ obj.constraints, ...
                                        H(:,:,j)*obj.alpha{k}+obj.V.A*(obj.xa{k+1}-obj.Ak_v(:,:,j)*obj.xa{k})+obj.V.A*obj.Bv(:,:,j)*obj.c{k}+...
                                        obj.V.A*obj.Akn_v(:,:,j)*(obj.xN{k}-obj.xaN{k})+obj.V.A*obj.Anv(:,:,j)*obj.xaN{k}+obj.w_bar<=obj.alpha{k+1}]; % Tube evoluition
                    end
                                        %% Also tightening depends on the neighbors states
                    for j=1:obj.vC
                    obj.constraints = [ obj.constraints, ...
                                        HC(:,:,j)*obj.alpha{k}-Gb*obj.K*obj.xa{k}+Gb*obj.c{k}+Gb*obj.ua{k}+...
                                        (Gb*obj.Kn)*obj.xN{k}-(Gb*obj.Kn)*obj.xaN{k}<=ones(size(Gb,1),1)];         % Constraint tightening                                                                                        % Tracking constraints
                    end    
                end

                % Terminal tube must depend on the neighbors states
                for j=1:obj.v
                    obj.constraints = [ obj.constraints, ...
                                        H(:,:,j)*obj.alpha{end}+obj.V.A*(eye(obj.n)-obj.Ak_v(:,1:obj.n,j))*obj.xa{end}...
                                        -obj.V.A*obj.Bv(:,:,j)*obj.Kn*obj.xaN{end}+obj.w_bar<=obj.alpha{end}];         % Terminal Tube
                end
                for j=1:obj.vC
                    obj.constraints = [ obj.constraints, ...
                                        HC(:,:,j)*obj.alpha{end}-Gb*obj.K*obj.xa{end}-Gb*obj.Kn*obj.xaN{end}+Gb*obj.ua{end}<=ones(size(Gb,1),1)]; % Terminal constraints
                end
                obj.constraints = [obj.constraints, ...   
                                ismember([obj.xa{end};obj.ua{end};obj.xaN{end}],obj.Zs),...
                                [A-eye(obj.n),B,An]*[obj.xa{end};obj.ua{end};obj.xaN{end}] ==zeros(obj.n,1),...
                                ...[A-eye(obj.n) B An;obj.C zeros(obj.p,obj.m) zeros(obj.p,obj.nN)]*[obj.xr{end};obj.ur{end};obj.xrN{end}] == [zeros(obj.n,1);obj.r(:,end)]
                                ...obj.C*obj.xr{end}==obj.r(:,end)
                                ];            % Artificial equilibrium point stationarity
            else
                for k=1:obj.N
                    obj.constraints = [obj.constraints, ...
                                    obj.x{k+1} == Ak*(obj.x{k}-obj.x{end}) + B*obj.c{k}+B*obj.u{end}+ Akn*(obj.xN{k}-obj.xN{end}) + A*obj.x{end} + An*obj.xN{end},...
                                    obj.u{k} == obj.K*(obj.x{k}-obj.x{end})+obj.Kn*(obj.xN{k}-obj.xN{end})+obj.c{k}+obj.u{end}];  % System dynamics
                
                                    %% Tube evolution must depend on the neighbors states
                    for j=1:obj.v
                    obj.constraints = [ obj.constraints, ...
                                        H(:,:,j)*obj.alpha{k}+obj.V.A*(eye(obj.n)-obj.Ak_v(:,:,j))*obj.x{end}+obj.V.A*obj.Bv(:,:,j)*obj.c{k}+...
                                        obj.V.A*obj.Akn_v(:,:,j)*(obj.xN{k}-obj.xN{end})+obj.V.A*obj.Anv(:,:,j)*obj.xN{end}+obj.w_bar<=obj.alpha{k+1}]; % Tube evoluition
                    end
                                        %% Also tightening depends on the neighbors states
                    for j=1:obj.vC
                    obj.constraints = [ obj.constraints, ...
                                        HC(:,:,j)*obj.alpha{k}-Gb*obj.K*obj.x{end}+Gb*obj.c{k}+Gb*obj.u{end}+...
                                        (Gb*obj.Kn)*obj.xN{k}-(Gb*obj.Kn)*obj.xN{end}<=ones(size(Gb,1),1)];         % Constraint tightening                                                                                        % Tracking constraints
                    end    
                end
                
                if obj.track 
                    % Terminal tube must depend on the neighbors states
                    for j=1:obj.v
                        obj.constraints = [ obj.constraints, ...
                                            H(:,:,j)*obj.alpha{end}+obj.V.A*(eye(obj.n)-obj.Ak_v(:,1:obj.n,j))*obj.x{end}...
                                            -obj.V.A*obj.Bv(:,:,j)*obj.Kn*obj.xN{end}+obj.w_bar<=obj.alpha{end}];         % Terminal Tube
                    end
                    for j=1:obj.vC
                        obj.constraints = [ obj.constraints, ...
                                            HC(:,:,j)*obj.alpha{end}-Gb*obj.K*obj.x{end}-Gb*obj.Kn*obj.xN{end}+Gb*obj.u{end}<=ones(size(Gb,1),1)]; % Terminal constraints
                    end
                    obj.constraints = [obj.constraints, ...   
                                    ...obj.xN{obj.N+1} == obj.xN{end}, ...
                                    ismember([obj.xa;obj.ua;obj.xaN],obj.Zs),...
                                    [A-eye(obj.n),B,An]*[obj.xa;obj.ua;obj.xaN] ==zeros(obj.n,1),...
                                    ...[obj.xr;obj.ur,obj.xrN]==zeros(obj.n+obj.m+obj.nN,1)];
                                    ...[A-eye(obj.n) B An;obj.C zeros(obj.p,obj.m) zeros(obj.p,obj.nN)]*[obj.xr;obj.ur;obj.xrN] == [zeros(obj.n,1);obj.r]
                                    ];            % Artificial equilibrium point stationarity
                end
            end
        end
        function [Inputs,Outputs] = setArguments(obj)
        
            [Inputs,Outputs] = setArguments@yDMPC(obj);
            Inputs{end+1} = obj.th;                                  % Position end-1
            Inputs{end+1} = obj.thC;                                % Position end
        end
%         function obj = setController(obj, Options)
%             Inputs = {obj.x{1}};                                    % Position 1
%             Outputs = {obj.J,obj.K_cost,obj.ADMM_cost};             % Positions 1  2  3
%             if iscell(obj.xa)&&iscell(obj.xr)
%                 Outputs{end+1} = [obj.x{:}, obj.xa{:}, obj.xr{:}];  % Position 4
%                 Outputs{end+1} = [obj.u{:}, obj.ua{:}];             % Position 5
%             elseif isa(obj.xr,'sdpvar')&&isa(obj.xa,'sdpvar')
%                 Outputs{end+1} = [obj.x{:}, obj.xa, obj.xr];        % Position 4
%                 Outputs{end+1} = [obj.u{:}, obj.ua];                % Position 5
%             end
%             if ~isempty(obj.An)
%                 Inputs{end+1} = obj.xN{1};                          % Position 2
%                 Inputs{end+1} = obj.xN_bar;                         % Position 3
%                 Inputs{end+1} = obj.lambdax;                        % Position 4
%                 if iscell(obj.xaN)&&iscell(obj.xrN)
%                     Outputs{end+1} = [obj.xN{:}, obj.xaN{:}, obj.xrN{:}];                        % Position 6
%                 elseif isa(obj.xrN,'sdpvar')&&isa(obj.xaN,'sdpvar')
%                     Outputs{end+1} = [obj.xN{:}, obj.xaN, obj.xrN];                        % Position 6
%                 end
%             end
%             if ~isempty(obj.Bn)
%                 Inputs{end+1} = obj.uN_bar;                         % Position 2/5
%                 Inputs{end+1} = obj.lambdau;                        % Position 3/6
%                 Outputs{end+1} = [obj.uN{:}, obj.uaN];                       % Position 6/7
%             end
%             
%             Inputs{end+1} = obj.r;                                  % Position end-3
%             Inputs{end+1} = obj.rho;                                % Position end-2
%             Inputs{end+1} = obj.th;                                  % Position end-1
%             Inputs{end+1} = obj.thC;                                % Position end
%             
%             obj.controller = optimizer(obj.constraints, obj.J, Options, Inputs, Outputs);        
%         end

        function obj = tightConstraints(obj)

            Fb=cat(1,repmat(obj.F,1,1,obj.vC),obj.Fc*obj.Cv);
            Gb=vertcat(obj.G,zeros(size(obj.Fc,1),size(obj.G,2)));
            DF=Fb(:,:,1)+Gb*obj.K;
            if obj.qC>0
                DF=[DF;Fb(:,:,2:end)-Fb(:,:,1)];
            end
            % Parametric uncertainty [3]-[Section 5.1]
            Aineq=[]; bineq=zeros(obj.na*obj.vC,1);
            for k=1:obj.vC
                Aineq=[Aineq;-kron(eye(obj.na),obj.thC_v(k,:))];
            end
            Aeq=kron(obj.V.A',eye(obj.qC+1));
            beq=reshape(permute(DF,[3 2 1]),[(obj.qC+1)*obj.n 1 size(DF,1)]);
            %options=mskoptimset; options.Display='off';
            options = mskoptimset('MSK_IPAR_NUM_THREADS', 6,'MaxIter',1e4);
            for i=1:size(Fb,1)
                f=zeros(1,1,obj.vC);
                for j=1:obj.vC 
                    cost=kron(ones(1,obj.na),obj.thC_v(j,:))';
                    [hbar, f(j)]=linprog(cost, Aineq, bineq, Aeq, beq(:,:,i),[],[],options);
                     %[hbar, f(j)]=linprog(cost, Aineq, bineq, Aeq, beq(:,:,i),[],[],optimoptions('linprog','Display','off'));
                    if f(j)>=max(f)
                        obj.HCbar(:,:,i)=reshape(hbar,[obj.qC+1 obj.na]);
                    end
                end
            end
            obj.HCbar=permute(obj.HCbar,[3 2 1]);
            
            
        end
        
        function obj = tubeInclusion(obj)

            % With the following we want to derive the code for 
            % S[i] = {x| V*x[i] <= alpha[i]} and S[i+1] = {x| V*x[i+1] <= alpha[i+1]} (Multiplicative uncertainty) [4] [5]
            % or S[i] = {x| V*x[i] <= alpha[i]} and S[i+1](theta) = {x| V*x[i+1](theta) <= alpha[i+1]} (Parametric uncertainty) [3] [6]
            % such that S[i] \subseteq S[i+1] for all i=1,...,N-1
            
            % Simple Multiplicative uncertainty has been parametrized setting parameter as one in each vertex of the set i.e. Th.V=eye(q)
            % in this way it can be manages as a parametric uncertainty
            %options=mskoptimset; options.Display='off';
            options = mskoptimset('MSK_IPAR_NUM_THREADS', 6,'MaxIter',1e4);
            for i = 1 : obj.na
                [~, obj.w_bar(i,:)]=linprog(-obj.V.A(i,:),obj.W.A,obj.W.b,[],[],[],[],options);
                 %[~, obj.w_bar(i,:)]=linprog(-obj.V.A(i,:),obj.W.A,obj.W.b,[],[],[],[],optimoptions('linprog','Display','off'));
            end
        
            % Parametric uncertainty [3]-[Section 5.1]
            Aineq=[]; bineq=zeros(obj.na*obj.v,1);
            for k=1:obj.v
                Aineq=[Aineq;-kron(eye(obj.na),obj.th_v(k,:))];
            end
            Aeq=kron(obj.V.A',eye(obj.q+1));
            beq=reshape(permute(pagemtimes(obj.V.A,obj.Ak(:,1:obj.n,:)),[3 2 1]),[(obj.q+1)*obj.n 1 obj.na]);
            %options=mskoptimset; options.Display='off';
            for i=1:obj.na
                f=zeros(1,1,obj.v);
                for j=1:obj.v 
                    cost=kron(ones(1,obj.na),obj.th_v(j,:))';
                    [hbar, f(j)]=linprog(cost, Aineq, bineq, Aeq, beq(:,:,i),[],[],options);
                     %[hbar, f(j)]=linprog(cost, Aineq, bineq, Aeq, beq(:,:,i),[],[],optimoptions('linprog','Display','off'));
                    if f(j)>=max(f)
                        obj.Hbar(:,:,i)=reshape(hbar,[obj.q+1 obj.na]);
                    end
                end
            end
            obj.Hbar=permute(obj.Hbar,[3 2 1]);


        end
        
        function H = paramEval(obj,A, th)
            if isa(th,"single")||isa(th,"double")
                H=permute(pagemtimes(th,permute(A,[3 1 2])),[2 3 1]);
            elseif isa(th,'sdpvar')
                H=0;
                for i=1:length(th)
                    H=H+th(i)*A(:,:,i);
                end
            end
        end

        function [th_v, v] = setParamSet(obj,Th)
            th_v = [ones(size(Th.V,1),1) Th.V];        
            v=size(th_v,1);
        end
    end
end
