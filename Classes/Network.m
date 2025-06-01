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
% FILE: Network.m
%==========================================================================
% AUTHOR: Fabio Faliero
% 
% DESCRIPTION: This class defines the agents in a distributed system 
% 
% COPYRIGHT (c) 2024 Politecnico di Torino
% All rights reserved.


% The Network class is used to define a network of agents in a distributed system





classdef Network < handle
    properties
        agents
        type
        name
    end
      
    methods
        function obj = Network(type)
            if nargin == 0
                obj.type = 'centralized';
            else
                obj.type = type;
            end
            obj.agents = {};
            if strcmp(obj.type,'centralized')
                obj.name="C_"+num2str(randi(1000));
            else
                obj.name="D_"+num2str(randi(1000));
            end
        end
        
        function obj=addAgent(obj, Agent)
            obj.agents{end+1} =  Agent;
        end

        function [X,Y,U] = run(obj,T,r,simType)
            
            %preallocate memory
            X = zeros(size(obj.agents{1}.sys.A,1), T/obj.agents{1}.Ts, length(obj.agents));
            Y = zeros(size(obj.agents{1}.sys.C,1), T/obj.agents{1}.Ts, length(obj.agents));
            U = zeros(size(obj.agents{1}.sys.B,2), T/obj.agents{1}.Ts, length(obj.agents));
%             R = zeros(size(obj.agents{1}.sys.C,1), T/obj.agents{1}.Ts, length(obj.agents));
            % Check agents integration step size Ts and if they are different override them
            % using the minimum one and print a warning
            Ts = obj.agents{1}.Ts;
            for i = 2:length(obj.agents)
                if obj.agents{i}.Ts ~= Ts
                    warning('Agents %d integration step size is different. Overriding with the minimum one', i);
                    Ts = min(Ts, obj.agents{i}.Ts);
                end
            end

            % Check agents control step size Tc and if they are different override them
            % using the maximum one and print a warning
            Tc = obj.agents{1}.Tc;
            for i = 2:length(obj.agents)
                if obj.agents{i}.Tc ~= Tc
                    warning('Agents %d control step size is different. Overriding with the maximum one', i);
                    Tc = max(Tc, obj.agents{i}.Tc);
                end
            end    
            
            % reference signal must be a vector of function handles of times (One function for each agent)

            if nargin<3
                r =cell(1,length(obj.agents));
            else
                u=cell(1,length(obj.agents));

                % Initialize every input u as a zero of input dimension coming from the agent B matrix
                for i = 1:length(obj.agents)
                    u{i} = zeros(size(obj.agents{i}.sys.B,2),1);
                end
                % compute for the ith agent from time 0 to T
                for i = 1:length(obj.agents)
                    r{i} = r{i}(0:Ts:(T+10)); % overwrites the reference function handle with the reference signal at each time step
                end
            end
            t=0;

            for i = 1:length(obj.agents)
                X(:,1,i) = obj.agents{i}.x0;
            end
            switch simType
                case 'realTime'
                    % Real time simulation
                    u_span=cell(1,length(obj.agents));
                    u_old = u;

                    while t<=T/Ts

                        %save initial time
                        if mod(t, Tc/Ts) == 0
                            t0 = t;
                        end
                        
                        % Compute control signal for each agent

                            for i = 1:length(obj.agents)
                                % TO BE DONE: Distributed SetMembership  
                                if isa (obj.agents{i}.K, 'yRAMPC')||isa (obj.agents{i}.K, 'yDRAMPC')
                                    
                                    % Identify new parameter sets
                                    obj.agents{i}.identifyThetaSets();
                                   
                                    %Update the controller sets
                                    obj.agents{i}.K.Th = obj.agents{i}.Th;
                                    obj.agents{i}.K.ThC = obj.agents{i}.ThC;
                                end
                            end

                            for i = 1:length(obj.agents)
                                % Distributed control  
                                if ~isempty(obj.agents{i}.K.N)
                                    N=obj.agents{i}.K.N;
                                else
                                    N=0;
                                end
                                if size(obj.agents{i}.K.r,2)>1
                                    t_slice=t+1:t+N+1;
                                else
                                    t_slice=t+1;
                                end
                                u{i}=obj.agents{i}.control(r{i}(:,t_slice));%r{i}(:,t+1:t+N+1
                            end
                            cT=obj.agents{1}.K.time;
                            for i = 2:length(obj.agents)
                                cT=max(cT,obj.agents{i}.K.time);
                            end
                            cStep=ceil(cT/Ts);
                            
                        % update the system using the old input signal until the cStep   
                        
                        if t+cStep-t0 >= Tc/Ts
                            % reduce cStep leaving only one step free
                            cStep = Tc/Ts-(t-t0+1);
                            warning('Solution is not in real time');
                            rStep=1;
                            flag=false;
                        else
                            rStep = Tc/Ts-cStep-(t-t0); %remaining steps
                            flag=true;
                        end
                        cT=cStep*Ts;
                        for i = 1:length(obj.agents)
                            u_span{i}=repmat(u_old{i},1,cStep);
                        end


                        % if converged append from cStep+1 to the end the new control step (when  mod(t,(Tc/Ts))==0) the new control signal
                        if obj.checkConvergence()
                            flag=true;
                            u_old = u;
                            for i = 1:length(obj.agents)
                                u_span{i} = [u_span{i} repmat(u{i},1,rStep)];   % append the new control signal
                            end
                            cT = cT+rStep*Ts;
                        end
                        if flag
                            for i = 1:length(obj.agents)
                                [obj.agents{i}, X(:,t+2:t+1+cT/Ts,i), Y(:,t+1:t+cT/Ts,i), U(:,t+1:t+cT/Ts,i)] = obj.agents{i}.simulate(u_span{i},cT);
                            end
                            
                            for i = 1:length(obj.agents)
                                obj.agents{i}=obj.agents{i}.update();
                            end
                        t=t+cT/Ts;
                        end
                        

                    end
                otherwise 

            % Non Real time simulation
                    called=false;
                    while t<=T/Ts
                        % Compute control signal

                        if mod(t, Tc/Ts) == 0   % if time zero or multiple of Agents.Tc
                            if (isa (obj.agents{i}.K, 'yRAMPC')||isa (obj.agents{i}.K, 'yDRAMPC'))&&(~called)
                                for i = 1:length(obj.agents)
                                % TO BE CHECKED: Distributed SetMembership  
                                
                                    called=true;
                                    % Identify new parameter sets
                                    if t~=0
                                        obj.agents{i}.identifyThetaSets();
                                    end
                                    %Update the controller sets
                                    obj.agents{i}.K.Th = obj.agents{i}.Th;
                                    obj.agents{i}.K.ThC = obj.agents{i}.ThC;
                                end
                            end

                            for i = 1:length(obj.agents)
                                if ~isempty(obj.agents{i}.K.N)
                                    N=obj.agents{i}.K.N;
                                else
                                    N=0;
                                end
                                if size(obj.agents{i}.K.r,2)>1
                                    t_slice=t+1:t+N+1;
                                else
                                    t_slice=t+1;
                                end

                                u{i}=obj.agents{i}.control(r{i}(:,t_slice));%r{i}(:,t+1:t+N+1
                            end
                            
                            
                        end
                        
                        
                        if obj.checkConvergence()
                            for i = 1:length(obj.agents)
                                [obj.agents{i}, X(:,t+2,i), Y(:,t+1,i), U(:,t+1,i)] = obj.agents{i}.simulate(u{i},Ts);
                            end
                            if t==100
                                pause(0.0001)
                            end
%                              subplot(1,2,1)
%                              plot(t,X(1,t+1,i),'r.')
%                              subplot(1,2,1)
%                              plot(t,X(2,t+1,i),'r.')
%                              %plot(t,r{i}(:,t+1),'g*')
%                              hold on
%                              
                            for i = 1:length(obj.agents)
                                obj.agents{i}=obj.agents{i}.update();
                            end
                            fprintf('%d \n',t)
                            t=t+1;
                            called=false;
                            for j=1:length(obj.agents)
                
                                if isa (obj.agents{i}.K, 'yRAMPC') ||  isa(obj.agents{i}.K, 'yDRAMPC')
                                    Th_Vrep{j}(:,:,t)=obj.agents{j}.Th.V;
                                    th_hat{j}(:,t)=obj.agents{j}.K.th_hat; 
                                else
                                    th_hat{j}=[];
                                    Th_Vrep{j}=[];
                                end
                            end
                        end
                        
                    end
                    % pause(0.001);
                    if isa (obj.agents{1}.K, 'yRAMPC') ||  isa(obj.agents{1}.K, 'yDRAMPC')
                        save(obj.name+"_th_log.mat","Th_Vrep","th_hat")
                    end
            end


        end

        function flag = checkConvergence(obj)
            if strcmp(obj.type, 'centralized')
                flag=true;
            elseif strcmp(obj.type, 'distributed')
                %initialize flag as a vectpr of false
                flag = zeros(1,length(obj.agents));
                for i = 1:length(obj.agents)
                    flag(i) = all([obj.agents{i}.K.gConvergence,obj.agents{i}.K.resetted]);  %% == obj.agents{i}.K.Ni;
                end
                flag= all(flag,"all");
                if flag
                    for i=1:length(obj.agents)
                        obj.agents{i}.K.resetted=false;
                    end
                end
            end
        end
    end
end




