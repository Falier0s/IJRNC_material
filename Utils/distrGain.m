
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
% FILE: distrGain.m
%==========================================================================
% AUTHOR: Fabio Faliero
% 
% DESCRIPTION: Distributed terminal cost and stabilizing gain for interconnected systems.
% Literature: 
% [1] C. Conte, C.N. Jones, M. Morari, M.N. Zeilinger, "Distributed syntesis and stability of Cooperative Distributed model predictive control for linear systems", Automatica, 2016
% COPYRIGHT (c) 2024 Politecnico di Torino
% All rights reserved.


function [] = distrGain(mode,centr,varargin)
    % The function takes as input the DMPC DRMPC or DRAMPC objects and
    % computes the stabilizing gain for the interconnected system
    M=length(varargin);
    I=eye(M);

    tol=1e-5;
    
    % The centralized system has centr.n states and centr.m inputs
    % Define the lifting matrices T{j} V{j} W{j}  to select the subsystem states and inputs
    %initialize the cell
    T=cell(M,1);
    V=cell(M,1);
    W=cell(M,1);
    n= centr.n/M;
    for j=1:M
        % The matrix T{j} selects the states of the j-th subsystem i.e. x{j}=T{j}x
        % The matrix W{j} selects the states of the j-th subsystem and its neighbours i.e. [x{j-1};x{j};x{j+1}]=W{j}x but the first and last subsystems are treated differently
        T{j}=kron(I(j,:),eye(n));
        if j == 1
        I_slice = I(j:j+1, :); % usa solo le righe j e j+1
        elseif j == M
            I_slice = I(j-1:j, :); % usa solo le righe j-1 e j
        else
            I_slice = I(j-1:j+1, :); % usa le righe j-1, j e j+1
        end
        % Calcolo la matrice W{j} usando il prodotto di Kronecker
        W{j} = kron(I_slice, eye(n));
        % The matrix V{j} selects the inputs of the j-th subsystem i.e. u{j}=V{j}u
        V{j}=I(j,:);
        Q{j}=varargin{j}.Q;
        E{j}=sdpvar(n);
        
        Y{j}=sdpvar(varargin{j}.m,varargin{j}.n+varargin{j}.nN);
        
        En{j}=E{j};
        Eb{j}=E{j};
    end
   
    % Define the matrices E Y G F En Eb
    % E{j} correspond to the inv(P{j}) local terminal cost (i.e. Vf{j} = x{j}'*P{j}*x{j}) of the j-th subsystem and is a decision variable - sdpvar
    % En{j} correspond to W{j}*E{j}*W{j}' and depend on the decision variables E{j} 
    % Y{j} correspond to the [K{j-1},K{j},K{j+1}]*En{j} = Kn{j}*En{j} gain matrix of the j-th subsystem and is a decision variable - sdpvar
    % G{j} correspond to the local relaxation cost (i.e. Vr{j} = [x{j-1};x{j};x{j+1}]'*G{j}*[x{j-1};x{j};x{j+1}] depending on neighbouring subsystems) and is a decision variable - sdpvar
    % F{j} correspond to En{j}*G{j}*En{j} and depend on the decision variables E{j} and G{j}
    % Eb{j} = W{j}*T{j}'*E{j}*T{j}*W{j}' and depend on the decision variables E{j}  
    % Define the matrices S as decision variables that as the same size of F
    for i=1:M
        
        if i~=1
            En{i}=blkdiag(E{i-1},En{i});
            En{i-1}=blkdiag(En{i-1},E{i});
            Eb{i}=blkdiag(zeros(n),Eb{i});
            Eb{i-1}=blkdiag(Eb{i-1},zeros(n));
            Q{i}=blkdiag(zeros(varargin{i-1}.n),Q{i});
            Q{i-1}=blkdiag(Q{i-1},zeros(varargin{i}.n));
        end
        F{i}=sdpvar(varargin{i}.n+varargin{i}.nN,varargin{i}.n+varargin{i}.nN);
        R{i}=varargin{i}.R;
        S{i}=sdpvar(varargin{i}.n);
        if size(varargin{i}.A,3)>1
            An{i}=varargin{i}.Av;
            B{i}=varargin{i}.Bv;
            v{i}=varargin{i}.v;
        else
            An{i}=varargin{i}.A;
            B{i}=varargin{i}.B;
            v{i}=1;
        end
        

    end
    
    % Defining An{i} = [A{i-1} A{i} A{i+1}] (depending on the neighbouring subsystems that is different for the first and the last agents) 
    % As defined in [1] at Eq. 19a the LMIs can be defined such that
    % [Eb{i}+F{i}  En{i}*An{i}'+Y{i}'*B{i}' En{i}*Q{i}^(1/2) Y{i}'*R{i}^(1/2);
    %  An{i}*En{i}+B{i}*Y{i} E{i} zeros(size(E{i},1),size(Q{i},2)) zeros(size(E{i},1),size(R{i},2));
    %  (Q{i}^(1/2))*En{i} zeros(size(Q{i},1),size(E{i},2)) eye(size(Q{i},1)) zeros(size(Q{i},1),size(R{i},2));  
    %  (R{i}^(1/2))*Y{i} zeros(size(R{i},1),size(E{i},2)) zeros(size(R{i},1),size(Q{i},2)) eye(size(R{i},1))]>=0 
    % Then it follow the comditions of Eq. 21a and 21b of [1] respectively
    % F{i}<=S{i}
    % sum_{i=1}^{M} W{i}'*S{i}*W{i} <=0 <=> sum_{j\in N_i} T{j}*W{j}'*S{i}*W{j}*T{j}' <=0 \forall i=1,...,M where N_i is the set of neighbours of the i-th subsystem including the subsystem itself
 
    % Define S as a blkdiag matrix of sdpvar
    % in detail is blkdiag(sdpvar(varargin{j-1).n),sdpvar(varargin{j}.n),sdpvar(varargin{j+1}.n)) for j=1,...,M first and last has to be treated differently
    % and define An{j} as horizontal concatenation of A{j-1} A{j} A{j+1} for j=1,...,M except for the first and last subsystems

    C=0;
    J=0;
    for i=1:M
    
        if i~=1
            S{i}=blkdiag(sdpvar(varargin{i-1}.n),S{i});
            S{i-1}=blkdiag(S{i-1},sdpvar(varargin{i}.n));
            if size(varargin{i}.A,3)>1
                An{i}=[varargin{i}.Anv(:,1:varargin{i-1}.n,:),An{i}];
                An{i-1}=[An{i-1},varargin{i-1}.Anv(:,end-varargin{i}.n+1:end,:)];
            else
                An{i}=[varargin{i}.An(:,1:varargin{i-1}.n),An{i}];
                An{i-1}=[An{i-1},varargin{i-1}.An(:,end-varargin{i}.n+1:end)];
            end
        end
        
    end

    % Define the LMI constraints
    for i=1:M
        C=C+W{i}'*S{i}*W{i};
    end

   C=[C<=0];


    for i=1:M
        
        C=[C, F{i}<=S{i}];    
        C=[C, E{i}>=tol*eye(size(E{i}))];      % positive definite
        for j=1:v{i}
            C=[C, [Eb{i}+F{i}  [En{i}*An{i}(:,:,j)'+Y{i}'*B{i}(:,:,j)' En{i}*Q{i}^(1/2) Y{i}'*R{i}^(1/2)];
                [En{i}*An{i}(:,:,j)'+Y{i}'*B{i}(:,:,j)' En{i}*Q{i}^(1/2) Y{i}'*R{i}^(1/2)]' blkdiag(E{i},eye(varargin{i}.n+varargin{i}.nN+varargin{i}.m))]>=0];
        end   
        if strcmp(mode,'volume')
            J=J-log(det(E{i}));
        elseif strcmp(mode,'LQR')
            J=[];%J+(trace(E{i}));
        end
        
    end 

    options=sdpsettings('solver','mosek','verbose',0);

    K_LMI=optimize(C,J,options);
    check=0;
    centr.K=zeros(centr.m,centr.n);
    for i=1:M
        centr.P=blkdiag(centr.P,inv(value(E{i})));
        varargin{i}.P=inv(value(E{i}));
        KN{i}=value(Y{i}*inv(value(En{i})));
        K=value(Y{i})*inv(value(En{i}));
        G{i}= inv(value(En{i}))*value(F{i})*inv(value(En{i}));
        %check=check+W{i}'*value(F{i})*W{i};
        
        if i~=1
            varargin{i}.Kn= K(:,1:varargin{i-1}.n);
            centr.K(pastm+1:pastm+varargin{i}.m,pastn-varargin{i-1}.n+1:pastn-varargin{i-1}.n+size(K,2))=K;
            K(:,1:varargin{i-1}.n)=[];
            pastm=pastm+varargin{i}.m;
            pastn=pastn+varargin{i}.n;
        elseif i==1
            centr.K(1:varargin{i}.m,1:size(K,2))=K;
            pastn=varargin{i}.n;
            pastm=varargin{i}.m;
        end
        varargin{i}.K=K(:,1:varargin{i}.n);
        K(:,1:varargin{i}.n)=[];
        varargin{i}.Kn=[varargin{i}.Kn,K];
    end

    v=size(centr.A,3);
    for i=1:v
        rho(i)=max(abs(eig(centr.A(:,:,i)+centr.B(:,:,i)*centr.K)));
    end
    centr.lam=max(rho)+(0.999-max(rho)).*9/10;
    for i=1:M
        varargin{i}.lam=max(rho)+(0.999-max(rho)).*10/10;
    end

    
end