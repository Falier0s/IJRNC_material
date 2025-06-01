
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
% FILE: yRMPC.m
%==========================================================================
% AUTHOR: Fabio Faliero
% 
% DESCRIPTION: Robust Linear Quadratic MPC using YALMIP
% This function is intended to build a simple MPC controller using YALMIP
% for both regulation and tracking problems using Polytopic Tubes for system 
% subject to additive, multiplicative or parametric uncertainty or a combination of them.
% Literature:
% [1] J. Kholer, E. Andina, R. Soloperto, M. A. Muller, F. Allgower, "Linear robust adaptive model predictive control: Computational complexity and conservatism - extended version", 2020
% [2] M. V. Kothare, V. Balakrishnan, M. Morari, "Robust constrained model predictive control using linear matrix inequalities", Automatica, 1996
% [3] T. Peschke, D. Gorges, "Robust adaptive tube tracking model predictive control for piece-wise constant reference signals", Int J Robust Nonlinear Control, 2023
% [4] J. Fleming, B. Kouvaritakis, M. Cannon, "Robust Tube MPC for Linear SystemsWith Multiplicative Uncertainty", IEEE Transactions on Automatic Control VOL.60 NO.4, April 2015
% [5] T. Peschke, D. Gorges, "Robust Tube-Based Tracking MPC for Linear Systems with Multiplicative Uncertainty", IEEE 58th Conference on Decision and Control (CDC), 2019
% [6] X. Lu, M. Cannon, "Robust Adaptive Tube Model Predictive Control", American Control Conference (ACC), 2019
% 
%
% COPYRIGHT (c) 2024 Politecnico di Torino
% All rights reserved.


classdef yRMPC < yMPC
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
        w_bar
        Hbar
        HCbar   
    end

    properties (Access=public)
        th_v                % vertices of the uncertain state parameter set (transposed)
        v                   % number of vertices of the uncertain state parameter set
        thC_v               % vertices of the uncertain output parameter set (transposed)
        vC                  % number of vertices of the uncertain output parameter set
        na                  % number of inequalities constraining the tube evolution
    end

    
    methods
        
        function obj = yRMPC(sys, Q, R, N, options)

            obj@yMPC(sys, Q, R, N, options); % To be checked probably must be initialized with the nominal system

            C_mask = squeeze(any(sys.C ~= 0, [1, 2]));
            AB_mask= squeeze(any([any(sys.A ~= 0, [1, 2]),any(sys.B ~= 0, [1, 2])],2));

            obj.A=obj.A(:,:,AB_mask);
            obj.B=obj.B(:,:,AB_mask);
            obj.C=obj.C(:,:,C_mask);
            
            obj.q=sum(AB_mask)-1;
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

            obj=obj.addHConstraints();
            if ~isempty(ThC)
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

        function obj = initialize(obj,mode)%
        
            % In the initialization we do all the OFFLINE computations
            % and we build the OCP problem to be solved online using YALMIP

            % Offline computations:
            % 1. Compute the tube gain K_tube and the terminal weight P    **OFFLINE:STEP1**
            % 2. Compute the tube shape V as a contractive set
            % 3. Compute the bounding worst case disturbace (i.e. w_bar),
            %    tube evolution condition (i.e. matrix H) 
            %    and the constraint tightening (i.e. matrix Hc)

            if all(obj.K==0,'all')
                obj=obj.stabGain(obj.Av,obj.Bv,obj.W.V,mode);  
            end
                                                             % 1. Compute the tube gain K_tube and the terminal weight P

            obj.Ak=obj.A+pagemtimes(obj.B,obj.K);

            obj.Ak_v=permute(pagemtimes(permute(obj.Ak,[1 3 2]),squeeze(obj.th_v)'),[1 3 2]);                                    % Define closed loop system evaluated at the vertices of the uncertain set
            %obj=obj.contractSet(obj.Ak_v,...
            %        Polyhedron([obj.F;obj.G*obj.K],ones(size([obj.F;obj.G*obj.K],1),1)));                                              % 2. Compute the tube shape V as a contractive set
            obj=obj.contractSet(obj.Ak_v,...[], ...
                                Polyhedron(obj.Z.A*[eye(obj.n);obj.K],obj.Z.b),obj.lam);%,Polyhedron());                                                % 2. Compute the tube shape V as a contractive set
            obj.V=obj.Xa;
            obj.na=length(obj.V.A);
            obj.alpha=sdpvar(repmat(obj.na,1,obj.N+1),ones(1,obj.N+1));
            obj=obj.tubeInclusion();                                                    % 3. Compute the bounding worst case disturbace (i.e. w_bar), tube evolution condition (i.e. matrix Hbar) 
            obj=obj.tightConstraints();                                                 % and the constraint tightening (i.e. matrix HCbar)
            % Same as yMPC
            % OCP construction and solver definition

            obj=obj.costBuild();
            obj=obj.sysBuild();            
            if strcmp(obj.solver,'gurobi')
                Options = sdpsettings('solver',obj.solver,'verbose',obj.verbose,'savesolveroutput',1,'warmstart',1);
            else
                Options = sdpsettings('solver',obj.solver,'verbose',obj.verbose,'savesolveroutput',1);
            end
            % Creazione del controller Yalmip
            obj=obj.setController(Options);
        end
        
        function obj = solve(obj, x0,ref)
            % Risoluzione del problema MPC
            obj.x0 = x0;
            if nargin < 3 || isempty(ref)
                ref = 0;
            end
            % [Sol,~,~,~,~,Data] = obj.controller(obj.x0,ref,[1;obj.Th.chebyCenter().x],[1;obj.ThC.chebyCenter().x]);
            [Sol,~,~,~,~,Data] = obj.controller(obj.x0,ref,[1;zeros(obj.Th.Dim,1)],[1;zeros(obj.ThC.Dim,1)]);%,AT
            obj=obj.warmstart();
            obj.time=Data.solvertime;
            obj.u_star = Sol{1};
            obj.J_star = Sol{2};
            obj.verboseMode();
            
        end
    end
    
    methods (Access=public)

        function obj = sysBuild(obj)

            % Constraints Parameters

            Ak=obj.paramEval(obj.Ak,obj.th');
            A=obj.paramEval(obj.A,obj.th');
            B=obj.paramEval(obj.B,obj.th');
            
            H=permute(pagemtimes(permute(obj.Hbar,[1 3 2]),squeeze(obj.th_v)'),[1 3 2]);
            HC=permute(pagemtimes(permute(obj.HCbar,[1 3 2]),squeeze(obj.thC_v)'),[1 3 2]);
            Gb=vertcat(obj.G,zeros(size(obj.Fc,1),size(obj.G,2)));

            obj.constraints = (obj.V.A*obj.x{1}<=obj.alpha{1}); % Initial Tube
            if iscell(obj.xa) && obj.track
                for k=1:obj.N
                    obj.constraints = [ obj.constraints, ...
                                        obj.x{k+1} == Ak*(obj.x{k}-obj.xa{k}) + B*obj.c{k}+obj.xa{k+1}, ...                      % System dynamics
                                        ...obj.x{k+1} == Ak*(obj.x{k}-obj.xa) + B*obj.c{k}+A*obj.xa+B*obj.ua, ...                      % System dynamics
                                        obj.Zs.A*[obj.xa{k};obj.ua{k}]<=obj.Zs.b,...
                                        ...[A-eye(obj.n),B]*[obj.xa{k};obj.ua{k}]==zeros(obj.n,1),... %Artificial equilibria
                                        obj.xa{k+1}==[A B]*[obj.xa{k};obj.ua{k}] ,...   % Artificial Trajectories
                                        ...[A B;obj.C zeros(obj.p,obj.m)]*[obj.xr{k};obj.ur{k}] == [obj.xr{k+1};obj.r(:,k)],...
                                        obj.u{k} == obj.K*(obj.x{k}-obj.xa{k})+obj.c{k}+obj.ua{k}];
                    for j=1:obj.v
                    obj.constraints = [ obj.constraints, ...
                                        H(:,:,j)*obj.alpha{k}+obj.V.A*(obj.xa{k+1}-obj.Ak_v(:,:,j)*obj.xa{k})+obj.V.A*obj.Bv(:,:,j)*obj.c{k}+obj.w_bar<=obj.alpha{k+1}]; % Tube evoluition
                    end
                    for j=1:obj.vC
                    obj.constraints = [ obj.constraints, ...
                                        HC(:,:,j)*obj.alpha{k}-Gb*obj.K*obj.xa{k}+Gb*obj.c{k}+Gb*obj.ua{k}<=ones(size(Gb,1),1)];         % Constraint tightening                                                                                          % Tracking constraints
                    end
                end

                for j=1:obj.v
                    obj.constraints = [ obj.constraints, ...
                                        H(:,:,j)*obj.alpha{end}+obj.V.A*(eye(obj.n)-obj.Ak_v(:,:,j))*obj.xa{end}+obj.w_bar<=obj.alpha{end}];         % Terminal Tube
                end
                for j=1:obj.vC
                    obj.constraints = [ obj.constraints, ...
                                        HC(:,:,j)*obj.alpha{end}-Gb*obj.K*obj.xa{end}+Gb*obj.ua{end}<=ones(size(Gb,1),1)];                       % Terminal constraints
                end
                obj.constraints = [ obj.constraints, ...
                                    ...obj.Zs.A*[obj.xa;obj.ua]<=obj.Zs.b,...
                                    ...[A-eye(obj.n),B]*[obj.xa;obj.ua]==zeros(obj.n,1),...
                                    obj.Zs.A*[obj.xa{end};obj.ua{end}]<=obj.Zs.b,...
                                    [A-eye(obj.n),B]*[obj.xa{end};obj.ua{end}]==zeros(obj.n,1),...
                                    ...[A-eye(obj.n),B;obj.C zeros(obj.p,obj.m)]*[obj.xr{end};obj.ur{end}]==[zeros(obj.n,1);obj.r(:,end)]
                                    ...obj.C*obj.xr{end}==obj.r(:,end)
                                    ];                                                  % Artificial equilibrium point
                
            else
                for k=1:obj.N
                    obj.constraints = [ obj.constraints, ...
                                        obj.x{k+1} == Ak*(obj.x{k}-obj.xa) + B*obj.c{k}+A*obj.xa+B*obj.ua, ...                      % System dynamics
                                        obj.u{k} == obj.K*(obj.x{k}-obj.xa)+obj.c{k}+obj.ua];
                    for j=1:obj.v
                    obj.constraints = [ obj.constraints, ...
                                        H(:,:,j)*obj.alpha{k}+obj.V.A*(eye(obj.n)-obj.Ak_v(:,:,j))*obj.xa+obj.V.A*obj.Bv(:,:,j)*obj.c{k}+obj.w_bar<=obj.alpha{k+1}]; % Tube evoluition
                    end
                    for j=1:obj.vC
                    obj.constraints = [ obj.constraints, ...
                                        HC(:,:,j)*obj.alpha{k}-Gb*obj.K*obj.xa+Gb*obj.c{k}+Gb*obj.ua<=ones(size(Gb,1),1)];         % Constraint tightening                                                                                          % Tracking constraints
                    end
                end
                
                if obj.track
                    for j=1:obj.v
                        obj.constraints = [ obj.constraints, ...
                                            H(:,:,j)*obj.alpha{end}+obj.V.A*(eye(obj.n)-obj.Ak_v(:,:,j))*obj.xa+obj.w_bar<=obj.alpha{end}];         % Terminal Tube
                    end
                    for j=1:obj.vC
                        obj.constraints = [ obj.constraints, ...
                                            HC(:,:,j)*obj.alpha{end}-Gb*obj.K*obj.xa+Gb*obj.ua<=ones(size(Gb,1),1)];                       % Terminal constraints
                    end
                    obj.constraints = [ obj.constraints, ...
                                        obj.Zs.A*[obj.xa;obj.ua]<=obj.Zs.b,...
                                        [A-eye(obj.n),B]*[obj.xa;obj.ua]==zeros(obj.n,1),...
                                        ...[A-eye(obj.n),B;obj.C zeros(obj.p,obj.m)]*[obj.xr;obj.ur]==[zeros(obj.n,1);obj.r]
                                        ];                                                  % Artificial equilibrium point
                end
            end

        end
        
        function obj = setController(obj, Options)
            Inputs={obj.x{1},obj.r,obj.th,obj.thC};%,obj.INV
            obj.controller = optimizer(obj.constraints, obj.J, Options, Inputs, {obj.u{1},obj.J});
        end

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
            %options=mskoptimset(); options.Display='off';
            options = mskoptimset('MSK_IPAR_NUM_THREADS', 6,'MaxIter',1e6);
            for i=1:size(Fb,1)
                f=zeros(1,1,obj.vC);
                for j=1:obj.vC 
                    cost=kron(ones(1,obj.na),obj.thC_v(j,:))';
                    %[hbar, f(j)]=linprog(cost, Aineq, bineq, Aeq, beq(:,:,i),[],[],options);
                     [hbar, f(j)]=linprog(cost, Aineq, bineq, Aeq, beq(:,:,i),[],[],optimoptions('linprog','Display','off'));
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
            %options=mskoptimset(); options.Display='off';
            options = mskoptimset('MSK_IPAR_NUM_THREADS', 6,'MaxIter',1e5);
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
            beq=reshape(permute(pagemtimes(obj.V.A,obj.Ak),[3 2 1]),[(obj.q+1)*obj.n 1 obj.na]);
            %options=mskoptimset();% options.Display='off';
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

        function obj = warmstart(obj)

            obj=warmstart@yMPC(obj);
            cellfun(@(x,y) warmstart(x,y),{obj.alpha{1:obj.N}},cellfun(@value,{obj.alpha{2:obj.N+1}},'UniformOutput',false),'UniformOutput',false);
            
        end

    end

end