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
% FILE: Agent.m
%==========================================================================
% AUTHOR: Fabio Faliero
% 
% DESCRIPTION: This class defines the agents in a distributed system 
% 
% COPYRIGHT (c) 2024 Politecnico di Torino
% All rights reserved.


% The Agent class is used to define a dynamic system in state-space form defined as:
%               x+ = Ax + Bu + sum(Anj*xj)




classdef Agent < handle
    properties
        dx
        u
        y
        x0
        x1
        K
        Ts      % Integration time
        Tc      % Control time
        xn
        neighbors
        zeta
        An % Lista delle matrici di interazione con i vicini
        sys % Sistema dinamico
        x_s
        z_1        
        Th
        ThC
    end
    properties (SetAccess = private)
        storage
        name
        old_dx
        prev_u
        q
        qC
        th
        thC
        y0
        W
        E
        simulMode
        AB_mask
    end
    
    methods
        function obj = Agent(sys, x0, Options)
            
            if isfield(Options, 'propagation')&&~isempty(Options.simulMode)
                obj.simulMode = Options.simulMode;
            else
                obj.simulMode = 'Discrete';
            end

            obj.Ts = Options.Ts;

            if strcmp(obj.simulMode,'Discrete') 
                sys.A = pagemtimes(sys.A,obj.Ts);
                sys.B = pagemtimes(sys.B,obj.Ts);
                sys.A(:,:,1)=sys.A(:,:,1)+eye(size(sys.A,1));
            end

                obj.sys = sys;      % The object contains the parametric model of the system
            

            if isfield(Options, 'W')&&~isempty(Options.W)
                obj.W = Options.W;
            else
                obj.W = Polyhedron();
            end

            if isfield(Options, 'E')&&~isempty(Options.E)
                obj.E = Options.E;
            else
                obj.E = Polyhedron();
            end

            if isfield(Options, 'name')&&~isempty(Options.name)
                obj.name = Options.name;
            else
                obj.name = [class(obj), '_', num2str(randi(1000))];
            end

            if isfield(Options, 'Th')&&~isempty(Options.Th)
                obj.Th = Options.Th;
            else
                obj.Th = Polyhedron();
            end
            if isfield(Options, 'ThC')&&~isempty(Options.thC)
                obj.ThC = Options.Th;
            else
                obj.ThC = Polyhedron();
            end
            obj.q=obj.Th.Dim;   
            obj.qC=obj.ThC.Dim;
            C_mask = squeeze(any(sys.C ~= 0, [1, 2]));
            obj.AB_mask= squeeze(any([any(sys.A ~= 0, [1, 2]),any(sys.B ~= 0, [1, 2])],2));
            if isfield(Options,'th')
                obj.th = [1;Options.th];
                sys.A= squeeze(pagemtimes(permute(sys.A(:,:,obj.AB_mask),[1 3 2]),obj.th));
                sys.B= squeeze(pagemtimes(permute(sys.B(:,:,obj.AB_mask),[1 3 2]),obj.th));
            else
                obj.th = [1;obj.Th.randomPoint()];
                sys.A= squeeze(pagemtimes(permute(sys.A(:,:,obj.AB_mask),[1 3 2]),obj.th));
                sys.B= squeeze(pagemtimes(permute(sys.B(:,:,obj.AB_mask),[1 3 2]),obj.th));
            end
            if isfield(Options,'thC')
                obj.thC = [1;Options.thC];
                sys.C= squeeze(pagemtimes(permute(sys.C(:,:,C_mask),[1 3 2]),obj.thC));
            else
                obj.thC = [1;obj.ThC.randomPoint()];
                sys.C= permute(pagemtimes(permute(sys.C(:,:,C_mask),[1 3 2]),obj.thC),[1 3 2]);
            end

            obj.x_s= [sys.A(:,:,1)-eye(size(sys.A,1)) sys.B(:,:,1);sys.C(:,:,1) sys.D(:,:,1)];
            obj.dx = @(t, x, u, d) sys.A(:,:,1)*x + sys.B(:,:,1)*u + d; % The system dynamics evaluated at real parameters
            obj.y = @(x,e) sys.C(:,:,1) * x+e;
            obj.x0 = x0;
            obj.neighbors = {};
            obj.An = [];
            obj.z_1=[obj.x0;zeros(size(sys.B,2),1)];
            obj.y0=sys.C(:,:,1)*x0+obj.E.randomPoint;
            obj.x1=x0;
        end

        function obj = setNeighbors(obj, neighbor, Anj)
            %Take the pointer of neighbor's object and the interaction matrix Anj
            obj.neighbors{end+1} = neighbor;
            obj.An = [obj.An Anj(:,:,obj.AB_mask)];
            if length(obj.neighbors) == 1
                obj.old_dx = obj.dx;
            end
            obj = obj.addNeighbourIngfluence();
            obj.z_1=vertcat(obj.z_1,neighbor.x0);
        end

        function obj = setController(obj, controller)
            obj.K = controller;
%             if isa(obj.K, 'yMPC')
            if isa(obj.K, 'yMPC')   
                obj.Tc=obj.K.Ts; % Set the control time
                if obj.Tc < obj.Ts
                    error('Control time must be a multiple of the integration time Ts');
                elseif mod(obj.Tc, obj.Ts) ~= 0
                    old_Ts = obj.Ts;
                    obj.Ts = obj.Tc / round(obj.Tc / obj.Ts);
                    warning('The integration time %f (s) Ts will be set as the nearest divisor of the control time Tc at %f (s)', old_Ts, obj.Ts);              
                end
            end
            if isa(obj.K, 'yDRMPC')
                obj.K=obj.K.getNeighbors(obj.neighbors, obj.An,[]);%,[-10 10].*ones(size(obj.An,2),1));%
%             elseif isa(obj.K, 'yDMPC')
%                 obj.K=obj.K.getNeighbors(obj.neighbors, squeeze(pagemtimes(permute(obj.An,[1 3 2]),[1;obj.Th.chebyCenter().x])));%,[-10 10].*ones(size(obj.An,2),1));%
            elseif isa(obj.K, 'yDMPC')
                obj.K=obj.K.getNeighbors(obj.neighbors, squeeze(pagemtimes(permute(obj.An,[1 3 2]),[1;obj.Th.chebyCenter().x])),[]);%,[-10 10].*ones(size(obj.An,2),1));%
            end
            
        end

        function u=control(obj,ref)
        
        % Controlla se il controller è una classe di MPC o una sottoclasse di essa
            
            if nargin < 2 || isempty(ref)
                ref = 0;
            end

            if isa(obj.K, 'yMPC')
%             if isa(obj.K, 'yMPC2')
                Input = {obj.x0,ref};
                if isa (obj.K, 'yRAMPC')||isa (obj.K, 'yDRAMPC')
                    Input = {obj.x0,ref,obj.Th,obj.ThC,obj.y0};
                end
                obj.K = obj.K.solve(Input{:});
                if isa (obj.K, 'yRAMPC')||isa (obj.K, 'yDRAMPC')
                    obj.z_1 = obj.K.z_1;
                end
                u = obj.K.u_star;
            else
                % setpoint
                if size(obj.x_s,2)>1
                    obj.x_s=obj.x_s\[zeros(size(obj.x0,1),1);ref];
                    obj.x_s=obj.x_s(1:size(obj.x0,1),1);
                end
                u = obj.K * (obj.x0-obj.x_s);
            end
        
        end

        function [obj, x, y, u] = simulate(obj,u,t)
            
           
            if ~isempty(obj.An)&&~isempty(obj.neighbors)
                %append the neighbors states to the input vector
                for i = 1:length(obj.neighbors)
                    u = [u; repmat(obj.neighbors{i}.x0,1,size(u,2))]; %#ok<AGROW>
                end
            end
            if strcmp(obj.simulMode,'Discrete') 
                x=obj.dx([],obj.x0,u,obj.W.randomPoint());
            else
                [~,x] = RungeKutta(@(t, x, u) obj.dx(t, x, u, obj.W.randomPoint()), [0 t], obj.x0, obj.Ts, u);
                x = x(:,2:end);
            end       
            u = u(1:end-size(obj.An,2),:);
            
            y = obj.y(x,obj.E.randomPoint);
            obj.y0=y;
            % Aggiorna lo stato
            obj.x1 = x(:,end);
        end
        function obj= update(obj)
            
%             obj.z_1=[obj.x0;u]; % Update the z_{k-1} for the next iteration
            obj.x0=obj.x1;
%             obj.y0=obj.y(obj.x0,obj.E.randomPoint()); % Add Output Error
        end

        function obj = identifyThetaSets(obj)
        %% TO BE DONE: Distributed SetMembership

            % Set Membership Identification of the agent PUI

            C_mask = squeeze(any(obj.sys.C ~= 0, [1, 2]));
            obj.AB_mask= squeeze(any([any(obj.sys.A ~= 0, [1, 2]),any(obj.sys.B ~= 0, [1, 2])],2));
            A=obj.sys.A(:,:,obj.AB_mask);
            B=obj.sys.B(:,:,obj.AB_mask);
            C=obj.sys.C(:,:,C_mask);
            if ~isempty(obj.An)
                An=obj.Tc.*obj.An;
            else
                An=[];
            end
            if ~strcmp(obj.simulMode,'Discrete')
                I=zeros(size(A,1),size(A,1),size(A,3));
                I(:,:,1)=eye(size(A,1));
                A = obj.Tc*A+I;
                B = obj.Tc* B;
                 %obj.sys.A
                C = C;
                 %obj.z_1=[obj.x0;obj.K.u_star];
            end
            
            D=squeeze(pagemtimes([A,B,An],obj.z_1));
            Dc=squeeze(pagemtimes(C,obj.x0));
        
            d=D(:,1)-obj.x0; %TO CHECK
            dC=Dc(:,1)-obj.y0;
        
            D=D(:,2:end);
            Dc=Dc(:,2:end);

            % Non Falsified Parameters Sets
            Phi=Polyhedron();
            PhiC=Polyhedron();
            if (obj.q > 0)
                Phi=Polyhedron(-obj.W.A*D,obj.W.b+obj.W.A*d);
            end
        
            if (obj.qC > 0)
                PhiC=Polyhedron(-obj.E.A*Dc,obj.E.b+obj.E.A*dC);
            end

            % Intersect Previous PUI with NFP sets
%             if ~isempty(obj.neighbors)
%                 paramIdxA = find(obj.K.mask(2:end,:));
%                 for j=1:length(obj.neighbors)
%                     paramIdxN=find(obj.neighbors{j}.K.mask(2:end,:));
%                     % compare parameters in common with the other agents
%                     [commonParam, commonIdxA, commonIdxN] = intersect(paramIdxA, paramIdxN);
%                     if ~isempty(commonParam)
% 
%                         H_exp=zeros(size(obj.Th.A,1),size(obj.neighbors{j}.Th.A,2));
%                         H_exp(:,commonIdxN)=obj.Th.A(:,commonIdxA);
%                         obj.neighbors{j}.K.Th=obj.neighbors{j}.Th&Polyhedron(H_exp,obj.Th.b);
%                         obj.neighbors{j}.K.Th=obj.neighbors{j}.Th.minHRep;
%                         H_exp=zeros(size(obj.neighbors{j}.Th.A,1),size(obj.Th.A,2));
%                         H_exp(:,commonIdxA)=obj.neighbors{j}.Th.A(:,commonIdxN);
%                         Th_til=obj.Th&Polyhedron(H_exp,obj.neighbors{j}.Th.b);
% %                         obj.K.Th=obj.Th.minHRep;
%                         
%                     end
%                 end
%             else
            Th_til=obj.Th;
%             end
            if ~Phi.isEmptySet
                Th_til=obj.Th&Phi;
            end
%             Th_til=Th_til&Phi;
            ThC_til=obj.ThC&PhiC;

            % Cartesian Product of the two sets To build a single LP problem

            S=obj.Th*obj.ThC;        h=size(S.A,1); %h=size([S.A;2*(dec2bin(0:2^S.Dim-1)-'0')-1],1);
            S_til=Th_til*ThC_til;    h_til=size(S_til.A,1);

            % Define the LP problem 
            % - c: objective function coefficients

            c=[ones(h,1);zeros(h_til*h,1)];

            % - A_in, b_in: inequality constraints

            A_in=[-eye(h),  kron(eye(h),S_til.b')
                  zeros(h_til*h,h), -eye(h_til*h)];
            b_in=[zeros(h,1);zeros(h_til*h,1)];
        
            % - A_eq, b_eq: equality constraints

            A_eq=[zeros(h*(obj.q+obj.qC),h),kron(eye(h),S_til.A')];
        
            b_eq=[vec(S.A')];
%             b_eq=[vec([S.A;2*(dec2bin(0:2^S.Dim-1)-'0')-1]')];

            % Solve the LP problem
            options = mskoptimset('MSK_IPAR_NUM_THREADS', 6,'MaxIter',1e8);
            %options=optimoptions('linprog','Display','off');
            Sol=linprog(c',A_in,b_in,A_eq,b_eq,[],[],options);
            v=size(obj.Th.A,1);
            vC=size(obj.ThC.A,1);
%             Th_up=Polyhedron();
%             ThC_up=Polyhedron();
            
            
            if ~isempty(Sol)
                Th_up=Polyhedron(obj.Th.A,Sol(1:v,:));
                ThC_up=Polyhedron(obj.ThC.A,Sol(v+1:v+vC,:));  
                if ~isempty(obj.neighbors)
                    A_tmp=[];
                    b_tmp=[];
                    paramIdxA = find(obj.K.mask(2:end,:));
                    for j=1:length(obj.neighbors)
                        paramIdxN=find(obj.neighbors{j}.K.mask(2:end,:));
                        % compare parameters in common with the other agents
                        [commonParam, commonIdxA, commonIdxN] = intersect(paramIdxA, paramIdxN);
                        if ~isempty(commonParam)
                            % Project the neighbour set on the agent parameter
                            for k=1:length(commonIdxN)
                                Th_tmp=obj.neighbors{j}.Th.projection(commonIdxN(k));
                                A_tmp(end+1:end+length(Th_tmp.A),commonIdxA(k))=Th_tmp.A;
                                b_tmp(end+1:end+length(Th_tmp.A),1)=Th_tmp.b;
                            end
                         end
                    end
                    ThN=Polyhedron(A_tmp,b_tmp);
                    if obj.Th.contains(ThN) && obj.Th.contains(Th_up) 
                        Th_int = Th_up&ThN;
                        if Th_up.isEmptySet || ThN.isEmptySet
                            if Th_up.isEmptySet
                                obj.Th=ThN;
                            elseif ThN.isEmptySet
                                obj.Th=Th_up;
                            end
                        elseif ~Th_int.isEmptySet
                           obj.Th=Th_up&ThN;
                        end
                    end
                elseif obj.Th.contains(Th_up)  
                 %if obj.Th.contains(Th_up)
                    if ~Th_up.isEmptySet
                        obj.Th=Th_up;
    
                    end
                end
                if obj.ThC.contains(ThC_up)
                    obj.ThC=ThC_up;
                end
            end
            if obj.Th.isEmptySet
                pause(0.00001)
            end

            %% Distributed update of Theta sets
            % This section is intended for Distributed Robust Adaptive MPC
            % where learned set information can be intersected with the ones
            % learned by neighbours and learn the set fastely.
            % In details, one want to build two sets Th_up which is the set learned by the agent
            % and Th_N which is the set joining the information from the neighbours.
            % Th_N must be build knowing that the neighbours parameters and their respective 
            % sets do not necessarily refears to the same parameters of the agent and the relative 
            % labeling of those parameter may be different in the parameter vector.     
            % For that reason each neighbour set must be projected on the agent parameter and than 
            % Th_N obtained with a coerent carthesian product of each parameter in the agent set.
            % Than one can define the resulting set as the intersection of Th_N and Th_up.
            
        end


        function obj = addNeighbourIngfluence(obj)
            % Aggiorna la funzione dx per includere l'influenza del nuovo vicino
            if strcmp(obj.simulMode,'Discrete')
                An=obj.Ts.*obj.An;
            else
                An=obj.An;
            end

            obj.dx = @(t, x, u, w) obj.old_dx(t, x, u(1:end-size(An,2)),w) + squeeze(pagemtimes(permute(An,[1 3 2]),obj.th)) * u(end-size(An,2)+1:end);
        end
    end
end


