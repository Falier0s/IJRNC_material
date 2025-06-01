
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
% FILE: yDRAMPC.m
%==========================================================================
% AUTHOR: Fabio Faliero
% 
% DESCRIPTION: Distributed Linear Quadratic MPC using YALMIP
% This function is intended to build a simple MPC controller using YALMIP
% for both regulation and tracking problems of a distributed system using ADMM. 
% 
% COPYRIGHT (c) 2024 Politecnico di Torino
% All rights reserved.




classdef yDRAMPC < yDRMPC
    properties
        D
        Dc
        th_hat
        thC_hat
        z_1
    end
    
    methods
        function obj = yDRAMPC(sys, Q, R, N, options)
            % Chiamata al costruttore della superclasse
            obj@yDRMPC(sys, Q, R, N, options);
            obj.th_hat = []; %obj.Th.chebyCenter().x;%zeros(size(obj.q+1,1),1);
            obj.thC_hat = []; % obj.ThC.chebyCenter().x;%zeros(size(obj.qC+1,1),1);
            obj.z_1 = [];%zeros(obj.n+obj.m,1);
        end
        function obj = solve(obj, x0, ref,Th, ThC,y0)

            switch obj.status 
                case 0
                    obj = obj.optimize(x0,ref,Th, ThC,y0);
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
        function obj = optimize(obj, x0,ref,Th, ThC,y0)
            
            if nargin<2 || isempty(ref)
                obj.x{end}=zeros(obj.n,1);
                obj.u{end}=zeros(obj.m,1);
                obj.xN{end}=zeros(obj.n*obj.Ni,1);
                obj.uN{end}=zeros(obj.m*obj.Ni,1);
            end
            obj.Th = Polyhedron(Th.A(any(Th.A~=0,2),:),Th.b(any(Th.A~=0,2)));%obj.Th&Th;
            obj.ThC = Polyhedron(ThC.A(any(ThC.A~=0,2),:),ThC.b(any(ThC.A~=0,2)));%obj.ThC&ThC;
            if isempty(obj.th_hat)||isempty(obj.thC_hat)
                obj.th_hat = obj.Th.chebyCenter().x;%zeros(size(obj.q+1,1),1);
                obj.thC_hat = obj.ThC.chebyCenter().x;%zeros(size(obj.qC+1,1),1);
            end
            obj = obj.reset(x0);
            obj = obj.filter(x0,y0);
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
            %sol_r=[permute(pagemtimes(permute(obj.A-eye(obj.n),[1 3 2]),[1;obj.th_hat]),[1 3 2]),...
            %       permute(pagemtimes(permute(obj.B,[1 3 2]),[1;obj.th_hat]),[1 3 2]),...
            %       permute(pagemtimes(permute(obj.An,[1 3 2]),[1;obj.th_hat]),[1 3 2]);...  permute(pagemtimes(permute(obj.Bn,[1 3 2]),[1;obj.th_hat]),[1 3 2]);
            %       permute(pagemtimes(permute(obj.C,[1 3 2]),[1;obj.thC_hat]),[1 3 2]),...
            %       zeros(obj.p,obj.m+obj.nN+obj.mN)]\[zeros(obj.n,1);ref];
            %In{end+1} = sol_r(1:obj.n,1);
            In{end+1} = ref;
            In{end+1} = obj.Rho;
            In{end+1} = [1;obj.th_hat];
            In{end+1} = [1;obj.thC_hat];
            % Solve the optimization problem
            % Optimization problem depends also on th and thC estimatess
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
            obj.z_1=[x0;obj.u_star;obj.x0N];
        end

        % Receive the the self variables from the neighbours and update the self consensus variables

        function obj = filter(obj,x0,y0)

            mu=.05;         
            D=[];
            Dc=[];
            if isempty(obj.z_1)
                obj.z_1=[x0;zeros(obj.m,1);obj.x0N];
            end

            for i=1:obj.q
                D=[D,[obj.A(:,:,i+1),obj.B(:,:,i+1),obj.An(:,:,i+1)]*obj.z_1];
            end
        
            for i=1:obj.qC
                Dc = [Dc,obj.C(:,:,i+1)*obj.x0]; 
            end
        
            H=blkdiag(obj.Th.A,obj.ThC.A);
            h=[obj.Th.b;obj.ThC.b];
            if obj.q > 0
                th_til(1:obj.q) = obj.th_hat(1:obj.q) + mu * D' * (obj.x0 - obj.paramEval([obj.A(:,:,2:end), obj.B(:,:,2:end), obj.An(:,:,2:end)], obj.th_hat(1:obj.q,:)') * obj.z_1);
            end
            
            if obj.qC > 0
                th_til(obj.q+1:end) = obj.thC_hat(obj.q+1:end) + mu * Dc' * (y0 - obj.paramEval(obj.C(:,:,2:end), obj.thC_hat(q+1:end,:)') * obj.x0);
            end
            %options=mskoptimset();
            options = mskoptimset('MaxIter',1e8);
            %th_var=sdpvar(obj.q+obj.qC,1);
            
             th_new=quadprog(eye(obj.q+obj.qC),-th_til,H,h,[],[],[],[],[obj.Th.chebyCenter().x;obj.ThC.chebyCenter().x],options);
             %th_new=quadprog(eye(obj.q+obj.qC),-th_til,H,h,[],[],[],[],[obj.Th.chebyCenter().x;obj.ThC.chebyCenter().x],optimoptions("quadprog","Display","off"));
            if ~isempty(th_new)
                obj.th_hat = th_new(1:obj.q);
                obj.thC_hat = th_new(obj.q+1:end);
            end
            
        end
    
    end
end
