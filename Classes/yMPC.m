
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
% FILE: yMPC.m
%==========================================================================
% AUTHOR: Fabio Faliero
% 
% DESCRIPTION: Linear Quadratic MPC using YALMIP
% This function is intended to build a simple MPC controller using YALMIP
% for both regulation and tracking problems of linear systems. 
% % Literature:
% [1] C. Dorea, J.C. Hennet, "(a,b)-invariant polyhedral sets for the linear discrete-time systems", 1999   
% [2] M. Hovd, S. Olaru, G. Bitsoris, "Low Complexity Constraint Control Using Contractive Sets", IFAC 19th World Congress, 2014
% [3] F. Bianchini, S. Miani, " Set-theoretic methods in control", System & Control: Foundations & Applications, 2008
%
% COPYRIGHT (c) 2024 Politecnico di Torino
% All rights reserved They belong to GD
% The crazy person who wrote this thing belongs to GD 4ever



classdef yMPC < handle
    properties
        A                   % State matrix
        B                   % Input matrix
        C                   % Output matrix
        Ak                  % Closed loop state matrix
        Q                   % State cost matrix
        R                   % Input cost matrix
        P                   % Terminal state cost matrix
        S
        T
        N                   % Prediction horizon
        Nc                  % Control horizon
        x0                  % Initial state
        Ts                  % Sampling time
        u                   % Input variable 
        x                   % State variable
        xr
        ur
        xa
        ua
        z                   % concatenated state and input variable
        r                   % Reference
        K                   % Feedback gain
        n                   % Number of states
        m                   % Number of inputs
        p                   % Number of outputs
        Mn                  % SVDecomposed matrix for [x_s;u_s]=Mn*nu
        Nn                  % SVDecomposed matrix     y_r=Nn*nu
        Mx
        X
        Xa
        Zs
        s
        U
        Y
        relax
        Z
        Zc
        nu
        svd_flag
        constraints         % Optimization Problem Constraints
        J                   % Optimization Problem Objective
        J_star              % Optimal Objective
        u_star              % Optimal input
        solver              % Optimization Solver
        verbose             % Verbose level
        controller          % YAMLIP optimizer object
        partSol             % Particularized solution   
    end

    properties (Access=public)
        
        timer               % Measeures time for the optimization process
        name
        time
        track
        sigma               
        F                   % Constraints matrix Fx+Gu<=1
        G                   % Constraints matrix Fx+Gu<=1
        Fc                  % Constraints matrix Fx+Gu<=1 for the output constraints
        lam                 % Contractivity parameter
        prjctr
    end

    
    methods
        
        function obj = yMPC(sys, Q, R, N, options)
            
            % initialize the MPC object with the system matrices and the cost
            % System matrices initialization
            if ~isfield(options,'Ts')||isempty(options.Ts)
                obj.A = sys.A;
                obj.B = sys.B;
                obj.C = sys.C;
                obj.Ts = 0;
            else
                I=zeros(size(sys.A,1),size(sys.A,1),size(sys.A,3));
                I(:,:,1)=eye(size(sys.A,1));
                obj.Ts = options.Ts;
                obj.A = obj.Ts*sys.A+I;
                obj.B = obj.Ts*sys.B;
                obj.C = sys.C;
            end
            % Prediction horizon initialization

            obj.N = N;

            % System dimensions

            obj.n=size(obj.A,1);
            obj.m=size(obj.B,2);
            obj.p=size(obj.C,1);

            % Decision Variables

            obj.x = sdpvar(repmat(obj.n,1,obj.N+1),ones(1,obj.N+1));
            obj.u = sdpvar(repmat(obj.m,1,obj.N),ones(1,obj.N));
             %obj.u{obj.N+1} = zeros(obj.m, 1);                           % Last input is not a decision variable

            obj.z = sdpvar(obj.n+obj.m,1);                               % Concatenated state and input sequence
            % Cost function initialization
            
            obj.J = 0;


            % Cost matrices initialization

            obj.Q = Q;
            obj.R = R;
            obj.S = Q.*1e4;%50*(obj.C(:,:,1)*Q*obj.C(:,:,1)');
            obj.T = R.*1e4;


            % Check for options

            % If no solver defined use quadprog

            if ~isfield(options, 'solver')||isempty(options.solver)
                obj.solver = 'quadprog';
            else
                obj.solver=options.solver;
            end
            % If no verbose defined use 0
            if ~isfield(options,'relax')||isempty(options.relax)
                obj.relax=false;
            else
                obj.relax=options.relax;
            end

            if ~isfield(options, 'verbose')||isempty(options.verbose)
                obj.verbose = 0;
            else 
                obj.verbose = options.verbose;
            end


            % Define reference signal and the artficial equilibrium variable as the last in the state and input sequence
            % If no reference defined use 0
            %obj.r=sdpvar(obj.p,1);
            %obj.xr=sdpvar(obj.n,1);
            %obj.ur=sdpvar(obj.m,1);

            if isfield(options, 'ref') && all(logical(options.ref))
                obj.track=true;
                %obj.xa = sdpvar(obj.n,1);
                %obj.ua = sdpvar(obj.m,1);
                if strcmp(lower(options.ref), 'traj') || strcmp(lower(options.ref), 'trajectory')
                    obj.xa = sdpvar(repmat(obj.n,1,obj.N+1),ones(1,obj.N+1));
                    obj.ua = sdpvar(repmat(obj.m,1,obj.N),ones(1,obj.N));
                    obj.xr = sdpvar(repmat(obj.n,1,obj.N+1),ones(1,obj.N+1));
                    obj.ur = sdpvar(repmat(obj.m,1,obj.N),ones(1,obj.N));
                    obj.r  = sdpvar(obj.p,obj.N+1);
                     %obj.r  = sdpvar(repmat(obj.p,1,obj.N+1),ones(1,obj.N+1));
                else   
                    obj.xa = sdpvar(obj.n,1);
                    obj.ua = sdpvar(obj.m,1);
                    obj.xr = sdpvar(obj.n,1);
                    obj.ur = sdpvar(obj.m,1);
                    obj.r=sdpvar(obj.p,1);
                end
            else
                obj.track=false;
                obj.xa = zeros(obj.n,1);
                obj.ua = zeros(obj.m,1);
                obj.xr = zeros(obj.n,1);
                obj.ur = zeros(obj.m,1);
                obj.r=sdpvar(obj.p,1);
            end

            if isfield(options, 'name')&&~isempty(options.name)
                obj.name = options.name;
            else
                obj.name = [class(obj), '_', num2str(randi(1000))];
            end
            if isfield(options, 'svd')&&~isempty(options.svd)
                obj.svd_flag = options.svd;
            else
                obj.svd_flag = false;
            end

            if isfield(options, 'Nc')&&~isempty(options.Nc)
                obj.Nc = options.Nc;
            else
                obj.Nc = N;
            end

            if ~isfield(options,'sigma')||isempty(options.sigma)
                obj.sigma = 0.95;
            else
                obj.sigma = options.sigma;
            end

            if ~isfield(options, 'K')||isempty(options.K)
                obj.K = zeros(size(sys.B, 2), size(sys.A, 1));
            end

            if isfield(options, 'xBounds')&&~isempty(options.xBounds)
                lxb=options.xBounds(:,1);
                uxb=options.xBounds(:,2);
                if length(lxb)==1&&length(uxb)~=obj.n
                    lxb=lxb*ones(obj.n,1);
                    uxb=uxb*ones(obj.n,1);
                end
            % Substitute inf with large scalar, i.e. 1e12    
                lxb(isinf(lxb))=-1e12;
                uxb(isinf(uxb))=1e12;
            else
               lxb=-1e12*ones(obj.n,1);
               uxb=1e12*ones(obj.n,1);
            end
            
            if isfield(options, 'uBounds')&&~isempty(options.uBounds)
                lub=options.uBounds(:,1);
                uub=options.uBounds(:,2);
                if length(lub)==1&&length(uub)~=obj.m
                    lub=lub*ones(obj.m,1);
                    uub=uub*ones(obj.m,1);
                end
            % Substitute inf with large scalar, i.e. 1e12    
               lub(isinf(lub))=-1e12;
               uub(isinf(uub))=1e12;
            else
               lub=-1e12*ones(obj.m,1);
               uub=1e12*ones(obj.m,1);
            end

            if isfield(options, 'yBounds')&&~isempty(options.yBounds)
                lyb=options.yBounds(:,1);
                uyb=options.yBounds(:,2);
                if length(lyb)==1&&length(uyb)~=obj.p
                    lyb=lyb*ones(obj.p,1);
                    uyb=uyb*ones(obj.p,1);
                end
            % Substitute inf with large scalar, i.e. 1e12    
               lyb(isinf(lyb))=-1e12;
               uyb(isinf(uyb))=1e12;
            else
                lyb=-1e12*ones(obj.p,1);
                uyb=1e12*ones(obj.p,1);
            end
            
            obj.X=Polyhedron('lb',lxb,'ub',uxb);
            obj.U=Polyhedron('lb',lub,'ub',uub);
            obj.Y=Polyhedron('lb',lyb,'ub',uyb);

            % Define the Constraint Space

            obj.Z = obj.X*obj.U;
            obj.Zs= (obj.sigma*obj.Z);
            obj.Zc=Polyhedron([ismember(obj.C(:,:,1)*obj.z(1:obj.n,:),obj.Y)]);
            
            %obj.svd_flag=false;
            obj.time=0;
        end
        
        function obj = addHConstraints(obj, varargin)
            % Add hard constraints to the optimization problem
            % obj.constraints = [obj.constraints, Hconstraints];
            if nargin < 2 || isempty(varargin)
               warning('No constraints defined, the function will build the maticies F, Fc and G.');
               obj.Fc=obj.Zc.A(:,1:obj.p)./obj.Zc.b;       obj.Fc=obj.Fc(abs(obj.Zc.b)~=1e12,:)./obj.Zc.b(abs(obj.Zc.b)~=1e12);
               obj.F=obj.Z.A(:,1:obj.n);                   obj.F=obj.F(abs(obj.Z.b)~=1e12,:)./obj.Z.b(abs(obj.Z.b)~=1e12);      
               obj.G=obj.Z.A(:,obj.n+1:obj.n+obj.m);       obj.G=obj.G(abs(obj.Z.b)~=1e12,:)./obj.Z.b(abs(obj.Z.b)~=1e12);    
            else
               if nargin>2 && islogical(varargin{end}) % Output constraint flag
                   flag=varargin{end}; % copy the variable
                   varargin=varargin{1:end-1}; % remove the last element
               end

               if nargin == 2
                   Hconstraints = varargin{1};
                   if isa(Hconstraints,'constraint')||isa(Hconstraints,'lmi')
                       if size(depends(Hconstraints),2)<obj.n+obj.m
                           error(['Invalid constraints definition. Constraints must be defined as one of the following:\n' ...
                               '1) A function of the concatenated state and input variable z, i.e., Hconstraints = [F G]*z <= b\n' ...
                               '2) A matrix A = [F G] and a vector b\n' ...
                               '3) A Polyhedron object of the z variable\n\n' ...
                               'Note: Check the dimensions of the Polyhedron object. If Z.A contains zeros in an entire column, ' ...
                               'a specific constraint for that variable must be added, bounded by 1e12, or using mode 2).']);    
                       end%
                       Hconstraints = Polyhedron([Hconstraints, obj.Z.A*obj.z<=obj.Z.b]); 
                   elseif isa(Hconstraints,'Polyhedron')
                       if size(Hconstraints.A,2)<obj.n+obj.m
                           error(['Invalid constraints definition. Constraints must be defined as one of the following:\n' ...
                               '1) A function of the concatenated state and input variable z, i.e., Hconstraints = [F G]*z <= b\n' ...
                               '2) A matrix A = [F G] and a vector b\n' ...
                               '3) A Polyhedron object of the z variable\n\n' ...
                               'Note: Check the dimensions of the Polyhedron object. If Z.A contains zeros in an entire column, ' ...
                               'a specific constraint for that variable must be added, bounded by 1e12, or using mode 2).']);    
                       end
                   end
               elseif nargin==3
                   A=varargin{1}; %#ok<*PROPLC>
                   b=varargin{2};
                   if size(A,2)<obj.n+obj.m
                       error(['Invalid constraints definition. Constraints must be defined as one of the following:\n' ...
                           '1) A function of the concatenated state and input variable z, i.e., Hconstraints = [F G]*z <= b\n' ...
                           '2) A matrix A = [F G] and a vector b\n' ...
                           '3) A Polyhedron object of the z variable\n\n' ...
                           'Note: Check the dimensions of the Polyhedron object. If Z.A contains zeros in an entire column, ' ...
                           'a specific constraint for that variable must be added, bounded by 1e12, or using mode 2).']);     
                   end
                   Hconstraints = [A, b];
                   Hconstraints = [Hconstraints, obj.Z.A*obj.z<=obj.Z.b];
                   Hconstraints=Polyhedron(Hconstraints);   
               end  
               if flag
                   obj.Zc=obj.Zc.intersect(Hconstraints);
                   % Process Fc matrix excluding infinite constraints
               else
                   obj.Z=obj.Z.intersect(Hconstraints);
                   % Process F and G matrices excluding infinite constraints
               end
            end
  
        end

        function obj = addSConstraints(obj, Sconstraints) %#ok<*INUSD>
            % TO BE DONE - SOFT CONSTRAINT MANAGING with slack variables?
            % Aggiunta di vincoli soft al problema di ottimizzazione
            % obj.constraints = [obj.constraints, constraints];
        end  

        function obj = initialize(obj,mode) 
            obj=obj.addHConstraints();
            if obj.svd_flag
                obj=obj.decompose();
            end
            if all(obj.K==0,'all')
                obj=obj.stabGain(obj.A,obj.B,[],mode);
            %    [obj.K,obj.P]=dlqr(obj.A,obj.B,obj.Q,obj.R);
            end
            %obj.lam=.99;
            obj.Ak = obj.A+obj.B*obj.K;
            if ~obj.svd_flag
                obj.Mn=eye(obj.n+obj.m);
            else
                obj.Mx=[eye(obj.n),zeros(obj.n,obj.m)]*obj.Mn;
    
            end
            L = [-obj.K eye(obj.m)]*obj.Mn;
            Aa=[obj.Ak obj.B*L;zeros(size(obj.Mn,2),obj.n) eye(size(obj.Mn,2))];
            
            obj=obj.contractSet(Aa, ...[],.99
                                Polyhedron([obj.Z.A*[eye(obj.n), zeros(obj.n,size(obj.Mn,2));obj.K, L];[zeros(size(obj.Z.A,1),obj.n) obj.Zs.A*obj.Mn]],[obj.Z.b;obj.Zs.b]),1);%,Polyhedron());
            
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
        end
        
        function obj = solve(obj, x0,ref)
            % Risoluzione del problema MPC
            obj.x0 = x0;
            
            if nargin < 3 || isempty(ref)
                ref = 0;
            end

            [Sol,~,~,~,~,Data] = obj.controller(obj.x0,ref);
            obj=obj.warmstart();
            obj.time=Data.solvertime;
            obj.u_star = Sol{1};
            obj.J_star = Sol{2};
            obj.verboseMode();
            
        end

    end
    
    methods (Access=public)

        function obj = costBuild(obj)

            StageC=0;
            TerminalC=0;
            OffsetC=0;

            if ~obj.track
                for k=1:obj.n
                    StageC=StageC + (obj.x{k})'*obj.Q*(obj.x{k}) + (obj.u{k})'*obj.R*(obj.u{k});
                end
                TerminalC = TerminalC+(obj.x{obj.N+1})'*obj.P*(obj.x{obj.N+1});
            else
                if iscell(obj.xa)
                    for k=1:obj.N
                        StageC = StageC + (obj.x{k}-obj.xa{k})'*obj.Q*(obj.x{k}-obj.xa{k}) + (obj.u{k}-obj.ua{k})'*obj.R*(obj.u{k}-obj.ua{k});
                        OffsetC = OffsetC + (obj.C*obj.xa{k}-obj.r(:,k))'*obj.C*(1e0.*obj.Q)*obj.C'*(obj.C*obj.xa{k}-obj.r(:,k));
                        %OffsetC = OffsetC + (obj.xa{k}-obj.xr{k})'*(1e0.*obj.Q)*(obj.xa{k}-obj.xr{k});
                        %OffsetC = OffsetC + (obj.xa-obj.xr{k})'*(1e2.*obj.Q)*(obj.xa-obj.xr{k});
                    end
                    TerminalC = TerminalC + (obj.x{obj.N+1}-obj.xa{end})'*obj.P*(obj.x{obj.N+1}-obj.xa{end});
                    %OffsetC = OffsetC + (obj.xa-obj.xr{end})'*(1e1.*obj.P)*(obj.xa-obj.xr{end});
                    %OffsetC = OffsetC + (obj.xa{end}-obj.xr{end})'*(1e0.*obj.P)*(obj.xa{end}-obj.xr{end});
                    OffsetC = OffsetC + (obj.C*obj.xa{end}-obj.r(:,end))'*obj.C*(5e0.*obj.Q)*obj.C'*(obj.C*obj.xa{end}-obj.r(:,end));
                else
                    for k=1:obj.N
                        StageC = StageC + (obj.x{k}-obj.xa)'*obj.Q*(obj.x{k}-obj.xa) + (obj.u{k}-obj.ua)'*obj.R*(obj.u{k}-obj.ua);
                    end
                    if obj.svd_flag
                        OffsetC = OffsetC + (obj.nu{1}-obj.nu{2})'*obj.Mx'*(5e1.*obj.Q)*obj.Mx*(obj.nu{1}-obj.nu{2});
                    else
                        OffsetC = OffsetC + (obj.C*obj.xa-obj.r)'*obj.C*(5e1.*obj.Q)*obj.C'*(obj.C*obj.xa-obj.r);
                       % OffsetC = OffsetC + (obj.xa-obj.xr)'*(1e4.*obj.P)*(obj.xa-obj.xr);
                    end
                    TerminalC = TerminalC + (obj.x{obj.N+1}-obj.xa)'*(obj.P)*(obj.x{obj.N+1}-obj.xa);
                end

            end

            obj.J = obj.J + StageC + TerminalC + OffsetC;
            
        end

        function obj = stabGain(obj,A,B,w_v,cost)
        
            %Compute the feedback gain and the terminal cost matrix

            X=sdpvar(obj.n,obj.n);
            Y=sdpvar(obj.m,obj.n);
            gamma=sdpvar;
            l=sdpvar;
            t=sdpvar;
            tol=1e-8;

            %Define the LMIs

            v=size(A,3);
            % Define the LMIs
            L1=[];            
            
            T=[X*obj.Q^.5         Y'*obj.R^.5
               zeros(obj.n)    zeros(obj.n,obj.m)];

            for j=1:v
                L1=[L1, [[X (A(:,:,j)*X+B(:,:,j)*Y)'; A(:,:,j)*X+B(:,:,j)*Y X] T;T' gamma*eye(obj.n+obj.m)]>=tol];
            end
            
            % Contractiveness

            L2=[];
            for j=1:v
                L2=[L2,[l*X A(:,:,j)*X+B(:,:,j)*Y; (A(:,:,j)*X+B(:,:,j)*Y)' l*X]>=tol];

            end

        
            %Constraint Satisfaction


            L3=[];

            Fbar=[obj.F;obj.Fc];

            Gbar=[obj.G;zeros(size(obj.Fc,1),size(obj.G,2))];

            for i = 1:size(Fbar,1)
                L3=[L3, [1 Fbar(i,:)*X+Gbar(i,:)*Y;(Fbar(i,:)*X+Gbar(i,:)*Y)' X]>=tol];
            end


            % Noise Attenuation

            L4=[];
            if ~isempty(w_v)

                for j=1:v
                    for k=1:size(w_v,1)
                        L4=[L4 [t*X zeros(obj.n,1) (A(:,:,j)*X+B(:,:,j)*Y)';zeros(1,obj.n) 1-t w_v(k,:);A(:,:,j)*X+B(:,:,j)*Y w_v(k,:)' X]>=tol];
                    end
                end
            end
            
            %Solve the LMI

           
            if strcmp(cost,'volume') % Minimum Volume Cost
                K_LMI=optimizer([X>=tol.*eye(size(X,1)),L1,L2,L3,L4,gamma>=tol],-log(det(X)),sdpsettings('solver','mosek','verbose',0),{l,t},{X,Y});
            elseif strcmp(cost,'performance') % Stabilizing Cost
                K_LMI=optimizer([X>=tol.*eye(size(X,1)),L1,L2,L3,L4,gamma>=tol],gamma,sdpsettings('solver','mosek','verbose',0),{l,t},{X,Y});
            elseif strcmp(cost,'LQR')
                K_LMI=optimizer([X>=tol.*eye(size(X,1)),L1,L2,L3,L4,gamma>=tol],trace(X),sdpsettings('solver','mosek','verbose',0),{l,t},{X,Y});
            end
            obj.lam=.95;
            tau=0.95;
            flag=true;
            while flag
            
                Sol=K_LMI(obj.lam,tau);
                obj.P=Sol{1}\eye(obj.n);
                obj.K=Sol{2}*obj.P;

                for i=1:v
                    rho(i)=max(abs(eig(A(:,:,v)+B(:,:,v)*obj.K)));
                end
                    
                    %rho
                
                if any(rho>=1)||any(isnan(rho))
                    obj.lam=obj.lam.*.9;
                    tau=tau.*.9;
                elseif obj.lam<=1.01*(3+max(rho))/4
                    flag=false;
                    obj.lam=max(rho);
                else
                    obj.lam=(3+max(rho))/4;
                    tau=tau.*.99;
                end
            end

        end

        function obj = decompose(obj)

            % In here we want to find the minimum parameter decomposition of the steady states tuple (x_s,u_s) 
            % such that those states can be associated to a given setpoint y_r. i.e. E*[x_s;u_s]=F*y_r


            E = [obj.A-eye(obj.n) obj.B;obj.C zeros(obj.p,obj.m)];
            F = [zeros(obj.n,obj.p);eye(obj.p)];

            [U,S,V]=svd(E);

            % Now we want to define the matricies M and N such that:
            % [x_s;u_s]=M*\nu and y_r=N*\nu  where \nu is the new set of minimal required parameters
            Vp=V(:,rank(S)+1:end); 
            Up=U(:,rank(S)+1:end);
            if size(U,2)==(obj.p+obj.n)
                G=eye(obj.p);
            elseif size(U,2)<(obj.p+obj.n)
                G=(F'*Up)%_p; % add the orthonormal complement of the matrix to be defined
                G=G(:,rank(S)+1:end);
            end

            if size(U,2)<(obj.m+obj.n)
                obj.Mn=[V*pinv(S)*U'*F*G Vp]; % Requires the orthonormal complement Vp to be defined
                obj.Nn=[G zeros(obj.p,obj.m+obj.n-size(U,2))];
            elseif size(U,2)==(obj.m+obj.n)
                obj.Mn=V*pinv(S)*U'*F*G;
                obj.Nn=G;
            end

            obj.nu=sdpvar(repmat(size(obj.Mn,2),1,2),ones(1,2));
            obj.svd_flag=true;
            
        end

        function obj = sysBuild(obj)

            % Constraints initialization
            if obj.relax
                obj.s=sdpvar(size([obj.F; obj.Fc * obj.C], 1), 1);
                obj.J=obj.J+0e0*sum(obj.s);
                obj.constraints = [obj.s>=0];
            else
                obj.s=zeros(size([obj.F; obj.Fc * obj.C], 1), 1);
                obj.constraints=[];
            end
            for k=1:obj.N
               obj.constraints = [ obj.constraints, ...
                                   [[obj.F;obj.Fc*obj.C],[obj.G;zeros(size(obj.Fc,1),size(obj.G,2))]]*[obj.x{k};obj.u{k}]<=ones(size([obj.F;obj.Fc*obj.C],1),1)+obj.s, ... % inequality constraints Fx+Gu<=1
                                   obj.x{k+1} == obj.A*obj.x{k} + obj.B*obj.u{k}];   % System dynamics
                                     
            end
            if obj.track

                if iscell(obj.xa)
                    for k=1:obj.N
                        obj.constraints = [obj.constraints, ...
                                            ...ismember([obj.x{obj.N+1};obj.xa{k};obj.ua{k}],obj.Xa),...
                                            obj.xa{k+1}==[obj.A obj.B]*[obj.xa{k};obj.ua{k}] ,... % Artificial trajectories
                                            ...[obj.A-eye(obj.n) obj.B]*[obj.xa{k};obj.ua{k}] == zeros(obj.n,1),... % Artificial equilibria
                                            ...[obj.A obj.B;obj.C zeros(obj.p,obj.m)]*[obj.xr{k};obj.ur{k}] == [obj.xr{k+1};obj.r(:,k)]
                                            ];
                    end
                    obj.constraints = [obj.constraints, ...
                                        ...ismember([obj.x{obj.N+1};obj.xa;obj.ua],obj.Xa),...
                                        ismember([obj.x{obj.N+1};obj.xa{end};obj.ua{end}],obj.Xa),...
                                        [obj.A-eye(obj.n) obj.B]*[obj.xa{end};obj.ua{end}] == zeros(obj.n,1),...
                                        ...[obj.A-eye(obj.n) obj.B]*[obj.xa;obj.ua] == zeros(obj.n,1),...
                                        ...[obj.A-eye(obj.n) obj.B;obj.C zeros(obj.p,obj.m)]*[obj.xr{end};obj.ur{end}] == [zeros(obj.n,1);obj.r(:,end)]];
                                        ...obj.C*obj.xr{end}==obj.r(:,end)
                                        ];
                elseif obj.svd_flag
                    obj.constraints = [ obj.constraints, ismember([obj.x{obj.N+1};obj.nu{1}],obj.Xa),...
                                                        [obj.xa;obj.ua]==obj.Mn*obj.nu{1},...
                                                        obj.r==obj.Nn*obj.nu{2}];
                else
                    obj.constraints = [obj.constraints, ismember([obj.x{obj.N+1};obj.xa;obj.ua],obj.Xa),...
                                                        [obj.A-eye(obj.n) obj.B]*[obj.xa;obj.ua] == zeros(obj.n,1),...
                                                        ...[obj.A-eye(obj.n) obj.B;obj.C zeros(obj.p,obj.m)]*[obj.xr;obj.ur] == [zeros(obj.n,1);obj.r]
                                                        ];
                end
            end

        end

        function obj= contractSet(obj,Aa,X0,lam)

            v=size(Aa,3);

            % Initialize the admissible Set 

            t=0;
            Contractive=false;
            na=size(Aa,2);
            
            obj.Xa=X0;
            A=eye(size(Aa,1));
            while ~Contractive
                t=t+1;
                A=pagemtimes(A,Aa./lam);
                Xa_1=Polyhedron([obj.Xa.A; ...
                                 reshape(pagemtimes(X0.A,permute(A,[1 3 2])),[],na)],...
                                 [obj.Xa.b;repmat(X0.b,v,1)]).minHRep();
                fprintf("Iteration: %d, Constraints: %d\n", t, size(Xa_1.A, 1));    
                if  Xa_1.eq(obj.Xa)
                    Contractive=true;
                end
                obj.Xa=Xa_1;
                 %if Contractive
                 %    obj.Xa.outerApprox();
                 %end
            end
        end

        function obj = setController(obj, Options)
            Inputs={obj.x{1},obj.r};
            Outputs={obj.u{1},obj.J};%,[obj.x{:}],[obj.u{:}]
            obj.controller = optimizer(obj.constraints, obj.J, Options, Inputs, Outputs);
        end

        function verboseMode(obj)
            % Save data in a file
            obj= obj.updateTime();
            obj.saveToFile();
            
            if obj.verbose == 1 % If verbose mode is on, print the data
                obj.printVerbose();
            end
        end
        
        function obj = updateTime(obj)
            %Empty
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
                fprintf(fid, 'OpTime, OpInput, OpCost\n');
            end
            
            fprintf(fid, '%f, %s, %f\n', obj.time, mat2str(obj.u_star), obj.J_star); % Write the data
            fclose(fid); % Close the file
        end
        
        function printVerbose(obj)
            fprintf(' [%s] |\t Solver: %s  |\t Optimization Time: %f |\t Optimal Input: %s |\t Optimal Cost: %f\n', obj.name, obj.solver, obj.time, mat2str(obj.u_star), obj.J_star);
        end

        function obj = warmstart(obj)

            cellfun(@(x,y) warmstart(x,y),{obj.x{1:obj.N}},cellfun(@value,{obj.x{2:obj.N+1}},'UniformOutput',false),'UniformOutput',false);
            cellfun(@(x,y) warmstart(x,y),{obj.u{1:obj.N-1}},cellfun(@value,{obj.u{2:obj.N}},'UniformOutput',false),'UniformOutput',false);

        end

    end
    
end