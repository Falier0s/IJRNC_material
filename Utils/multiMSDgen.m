
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
% FILE: multiMSDgen.m
%==========================================================================
% AUTHOR: Fabio Faliero
% 
% DESCRIPTION: This utility function generates a multiple Mass-Spring-Damper system
%              given the masses spring and damper vectors or only partial information 
%              to generate a random system with the given fixed properties. 
%              The function outputs the ss model of the system.
%              If uncertain system requested the function will output q+1 systems 
%              according to the parameter uncertainty.
% 
% COPYRIGHT (c) 2024 Politecnico di Torino
% All rights reserved.

function [ssModels, Th, ThC] = multiMSDgen(n, m, k, c, uncertainty)
    % FILE: multiMSDgen.m
    %==========================================================================
    % DESCRIPTION: This utility function generates a multiple Mass-Spring-Damper system
    %              given the masses, spring, and damper vectors or only partial information 
    %              to generate a random system with the given fixed properties. 
    %              The function outputs the state-space model of the system.
    %              If uncertain system requested, the function will output q+1 systems 
    %              according to the parameter uncertainty.
    
    % INPUT:
    %   n - Number of subsystems (masses)
    %   m - Vector or scalar of masses (1xN or scalar)
    %   k - Vector or scalar of spring constants (1x(N-1) or scalar)
    %   c - Vector or scalar of damper coefficients (1x(N-1) or scalar)
    %   uncertainty - Struct with fields 'mass', 'spring', 'damper', 'output' 
    %                 representing the percentage of uncertainty (e.g., 0.1 for 10%)
    %
    % OUTPUT:
    %   ssModels - Struct containing nominal and uncertain state-space models
    %   warningMsg - Any warning messages related to the deduction of system parameters
    
        % Initialize warning message
        warningMsg = '';
        
        % Deduce number of subsystems if n is not provided
        if isempty(n)
            if ~isempty(m)
                n = length(m);
            elseif ~isempty(k)
                n = length(k) + 1;
            elseif ~isempty(c)
                n = length(c) + 1;
            else
                error('Insufficient information to deduce the number of subsystems (n).');
            end
            warningMsg = sprintf('Warning: Number of subsystems (n) deduced as %d.', n);
        end
        
        % Handle scalar inputs for m, k, c by replicating the values
        if isscalar(m)
            m = repmat(m, 1, n);
        elseif length(m) ~= n
            error('Length of mass vector (m) must match the number of subsystems (n).');
        end
        
        if isscalar(k)&&(n>2)
            k = repmat(k, 1, n-1);
        elseif length(k) ~= (n-1)
            error('Length of spring constant vector (k) must be n-1.');
        end
        
        if isscalar(c)&&(n>2)
            c = repmat(c, 1, n-1);
        elseif length(c) ~= (n-1)
            error('Length of damper coefficient vector (c) must be n-1.');
        end
    
        % Create symbolic variables
        syms m_s [n 1] real 
        syms k_s [n-1 1] real
        syms c_s [n-1 1] real
        
        % Initialize symbolic matrices
        for i=1:n-1
            A_Aux=kron([-1 1;1 -1],[0 0;k_s(i) c_s(i)]);
            if i==1
                A=A_Aux;
            else
                A=blkdiag(A,zeros(2))+blkdiag(zeros(2*(i-1)),A_Aux);
            end
        end
        A=A./kron(m_s,ones(2,1))+kron(eye(n),[0 1;0 0]);
        B=kron(eye(n)./m_s,[0 1]');
        C=kron(eye(n),[1 0]);
        u_out=false;

        % Number of models (nominal + uncertainties)
        p = 1;
        if nargin == 5 && ~isempty(uncertainty)
            
            if isfield(uncertainty, 'mass') 
                if isscalar(uncertainty.mass)
                    uncertainty.mass=uncertainty.mass.*ones(n,1);
                elseif length(uncertainty.mass) ~= n
                    error("Length of mass uncertainty must be consistent with the number of agents.")
                end
                p=p+sum(uncertainty.mass~=0);
            else
                uncertainty.mass=zeros(n,1);
            end
            if isfield(uncertainty, 'spring') 
                if isscalar(uncertainty.spring)
                    uncertainty.spring=uncertainty.spring.*ones(n-1,1);
                elseif length(uncertainty.spring) ~= n-1
                    error("Length of spring uncertainty must be consistent with the number of agents.")
                end
                p=p+sum(uncertainty.spring~=0);
            else
                uncertainty.spring=zeros(n-1,1);
            end            
            if isfield(uncertainty, 'damper') 
                if isscalar(uncertainty.damper)
                    uncertainty.damper=uncertainty.damper.*ones(n-1,1);
                elseif length(uncertainty.damper) ~= n-1
                    error("Length of damper uncertainty must be consistent with the number of agents.")
                end
                p=p+sum(uncertainty.damper~=0);
            else
                uncertainty.damper=zeros(n-1,1);
            end    
               
            if isfield(uncertainty, 'output') 
                p = p + 1;
                u_out=true;
            end
        end
        
        th_m=1./m_s(uncertainty.mass~=0);
        th_k=[];
        th_c=[];
        for i=1:n-1
            if uncertainty.spring(i)~=0
                th_k=[th_k;k_s(i)/m_s(i);k_s(i)/m_s(i+1)];
            end
            if uncertainty.damper(i)~=0
                th_c(i)=[th_c;k_s(i)/m_s(i);k_s(i)/m_s(i+1)];
            end
        end
        th_k=unique(simplify(subs(th_k,m_s(uncertainty.mass==0),ones(size(m(uncertainty.mass==0),1),1))));
%         th_k=unique(simplify(subs(th_k,m_s(uncertainty.mass==0),m(uncertainty.mass==0))));
        th_c=unique(simplify(subs(th_c,m_s(uncertainty.mass==0),ones(size(m(uncertainty.mass==0),1),1))));
%         th_c=unique(simplify(subs(th_c,m_s(uncertainty.mass==0),m(uncertainty.mass==0))));
        th=[th_m;th_k;th_c];
        
        th_nom=double(subs(th,[m_s;k_s;c_s],[m;k;c]));
        th_u=[];
        th_l=[];
        th_u=double(subs(th,[m_s;k_s;c_s],[m.*(1-uncertainty.mass);k.*(1+uncertainty.spring);c.*(1+uncertainty.damper)]))-th_nom;
        th_l=double(subs(th,[m_s;k_s;c_s],[m.*(1+uncertainty.mass);k.*(1-uncertainty.spring);c.*(1-uncertainty.damper)]))-th_nom;


        A=subs(A,[m_s(uncertainty.mass==0);k_s(uncertainty.spring==0);c_s(uncertainty.damper==0)],[m(uncertainty.mass==0);k(uncertainty.spring==0);c(uncertainty.damper==0)]);
        for i=1:length(th)
            A(:,:,i+1)=subs(A(:,:,1)-subs(A(:,:,1),th(i),0),th(i),1);
            B(:,:,i+1)=subs(B(:,:,1)-subs(B(:,:,1),th(i),0),th(i),1);
        end
        A(:,:,1)=subs(A(:,:,1),[m_s;k_s;c_s],[m;k;c]);
        B(:,:,1)=subs(B(:,:,1),[m_s;k_s;c_s],[m;k;c]);
        A=double(A);
        B=double(B);
        
        C(:,:,2:size(A,3))=zeros(n,2*n,p-1);
        thC_l=[];
        thC_u=[];
        if u_out
            C(:,:,end+1)=C(:,:,1);
            thC_u=ones(size(C,1),1).*uncertainty.output;
            thC_l=-ones(size(C,1),1).*uncertainty.output;
        end


        ssModels=ss(A,B,C,zeros(n,n,size(A,3)));


        if ~isempty(th_u)&&~isempty(th_l)
            Th=Polyhedron('lb',th_l,'ub',th_u);
        else
            Th=Polyhedron();
        end
        if ~isempty(thC_u)&&~isempty(thC_l)
            ThC=Polyhedron('lb',thC_l,'ub',thC_u);
        else
            ThC=Polyhedron();
        end

    end

                    

                    