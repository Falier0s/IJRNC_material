
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
% DESCRIPTION: Robust Adaptive Linear Quadratic MPC using YALMIP
% This function is intended to build a simple MPC controller using YALMIP
% for both regulation and tracking problems using Polytopic Tubes for system 
% subject to additive, multiplicative or parametric uncertainty or a combination of them.
% Literature:
% [1] T. Peschke, D. Gorges, "Robust adaptive tube tracking model predictive control for piece-wise constant reference signals", Int J Robust Nonlinear Control, 2023
% [2] X. Lu, M. Cannon, "Robust Adaptive Tube Model Predictive Control", American Control Conference (ACC), 2019
% [3] Parsi
% [4] Kholer
% COPYRIGHT (c) 2024 Politecnico di Torino
% All rights reserved.


classdef yRAMPC < yRMPC
    properties
        D
        Dc
        th_hat
        thC_hat
        z_1
    end

    
    methods
        
        function obj = yRAMPC(sys, Q, R, N, options)
            
            obj@yRMPC(sys, Q, R, N, options);
            obj.th_hat = []; %obj.Th.chebyCenter().x;%zeros(size(obj.q+1,1),1);
            obj.thC_hat = []; % obj.ThC.chebyCenter().x;%zeros(size(obj.qC+1,1),1);
            obj.z_1 = [];%zeros(obj.n+obj.m,1);

        end

        function obj = solve(obj, x0,ref,Th, ThC,y0)

            % Update the Polyhedron of uncertainty
            obj.Th = Polyhedron(Th.A,Th.b);%obj.Th&Th;
            obj.ThC = Polyhedron(ThC.A,ThC.b);%obj.ThC&ThC;
            if isempty(obj.th_hat)||isempty(obj.thC_hat)
                obj.th_hat = obj.Th.chebyCenter().x;%zeros(size(obj.q+1,1),1);
                obj.thC_hat = obj.ThC.chebyCenter().x;%zeros(size(obj.qC+1,1),1);
            end
                % Filter the new estimate from the set
            obj.x0 = x0;
            obj = obj.filter(x0,y0);
             %obj.th_hat=obj.Th.chebyCenter().x;
             %obj.thC_hat=obj.ThC.chebyCenter().x;
            % Risoluzione del problema MPC
            
            if nargin < 3 || isempty(ref)
                ref = 0;
            end
             %AT=[eye(obj.n),zeros(obj.n,obj.m)]*inv([paramEval(obj.A,[1;obj.th_hat]')-eye(obj.n) paramEval(obj.B,[1;obj.th_hat]'); paramEval(obj.C,[1;obj.thC_hat]') zeros(obj.p,obj.m)]);
             %sol_r=[squeeze(pagemtimes(permute(obj.A,[1 3 2]),[1;obj.th_hat]))-eye(obj.n) squeeze(pagemtimes(permute(obj.B,[1 3 2]),[1;obj.th_hat]));squeeze(pagemtimes(permute(obj.C,[1 3 2]),[1;obj.thC_hat])) zeros(obj.p,obj.m)]\[zeros(obj.n,1);ref];
             %[Sol,~,~,~,~,Data] = obj.controller(obj.x0,sol_r(1:obj.n,1),[1;obj.th_hat],[1;obj.thC_hat]);
            [Sol,~,~,~,~,Data] = obj.controller(obj.x0,ref,[1;obj.th_hat],[1;obj.thC_hat]);
            obj=obj.warmstart();
            obj.time=Data.solvertime;
            obj.u_star = Sol{1};
            obj.J_star = Sol{2};
            obj.verboseMode();
            obj.z_1=[x0;obj.u_star]; % Update the z_{k-1} for the next iteration
            
        end

        function obj = filter(obj,x0,y0)

            mu=.05;         
            D=[];
            Dc=[];
            if isempty(obj.z_1)
                obj.z_1=[x0;zeros(obj.m,1)];
            end

            for i=1:obj.q
                D=[D,[obj.A(:,:,i+1),obj.B(:,:,i+1)]*obj.z_1];
            end
        
            for i=1:obj.qC
                Dc = [Dc,obj.C(:,:,i+1)*obj.x0]; 
            end
        
            H=blkdiag(obj.Th.A,obj.ThC.A);
            h=[obj.Th.b;obj.ThC.b];
            if obj.q > 0
                th_til(1:obj.q) = obj.th_hat(1:obj.q) + mu * D' * (obj.x0 - obj.paramEval([obj.A(:,:,2:end), obj.B(:,:,2:end)], obj.th_hat(1:obj.q,:)') * obj.z_1);
            end
            
            if obj.qC > 0
                th_til(obj.q+1:end) = obj.thC_hat(obj.q+1:end) + mu * Dc' * (y0 - obj.paramEval(obj.C(:,:,2:end), obj.thC_hat(q+1:end,:)') * obj.x0);
            end
            %options=mskoptimset();
            options = mskoptimset('MaxIter',1e6);
            th_new=quadprog(eye(obj.q+obj.qC),-th_til,H,h,[],[],[],[],[obj.Th.chebyCenter().x;obj.ThC.chebyCenter().x],options);%optimoptions('quadprog','Display','off')
            obj.th_hat = th_new(1:obj.q);
            obj.thC_hat = th_new(obj.q+1:end);
            
        end

    end
end