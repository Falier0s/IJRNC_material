%% MPC
clear 
clc
close all
yalmip('clear')

addpath('../Utils')
rng(754)

%% - Initialization
masses = ones(5, 1);
spring = [1.4 2.6 1.8 2.2]';
damper = [2 2 1.65 1.8]';
uncertainty.spring = ([.18 .23 .22 .15]' ./ spring);
n = 5;

[sys, Th, ~] = multiMSDgen(n, masses, spring, damper, uncertainty);

subSys1 = ss(sys.A(1:2, 1:2, :), sys.B(1:2, 1, :), sys.C(1, 1:2, :), sys.D(1, 1, :));
subSys2 = ss(sys.A(3:4, 3:4, :), sys.B(3:4, 2, :), sys.C(2, 3:4, :), sys.D(2, 2, :));
subSys3 = ss(sys.A(5:6, 5:6, :), sys.B(5:6, 3, :), sys.C(3, 5:6, :), sys.D(3, 3, :));
subSys4 = ss(sys.A(7:8, 7:8, :), sys.B(7:8, 4, :), sys.C(4, 7:8, :), sys.D(4, 4, :));
subSys5 = ss(sys.A(9:10, 9:10, :), sys.B(9:10, 5, :), sys.C(5, 9:10, :), sys.D(5, 5, :));

w_bnd=0.05;
e_bnd = 0;

Ts=0.25;
Tc=0.25;

W = Polyhedron('lb', -w_bnd * ones(2 * n, 1), 'ub', w_bnd * ones(2 * n, 1));
E = Polyhedron('lb', -e_bnd * ones(n, 1), 'ub', e_bnd * ones(n, 1));

th_real = [1.26 2.39 2.0 2.34]';
th_real = th_real - spring;

% th_real = Th.randomPoint;

sysNom = ss(sys.A(:, :, 1), sys.B(:, :, 1), sys.C(:, :, 1), sys.D(:, :, 1));
subNomSys1 = ss(sysNom.A(1:2, 1:2, :), sysNom.B(1:2, 1, :), sysNom.C(1, 1:2, :), sysNom.D(1, 1, :));
subNomSys2 = ss(sysNom.A(3:4, 3:4, :), sysNom.B(3:4, 2, :), sysNom.C(2, 3:4, :), sysNom.D(2, 2, :));
subNomSys3 = ss(sysNom.A(5:6, 5:6, :), sysNom.B(5:6, 3, :), sysNom.C(3, 5:6, :), sysNom.D(3, 3, :));
subNomSys4 = ss(sysNom.A(7:8, 7:8, :), sysNom.B(7:8, 4, :), sysNom.C(4, 7:8, :), sysNom.D(4, 4, :));
subNomSys5 = ss(sysNom.A(9:10, 9:10, :), sysNom.B(9:10, 5, :), sysNom.C(5, 9:10, :), sysNom.D(5, 5, :));
      
% Initialize 
opt1=agentOpt();    % Agent 1 initialize options
opt2=agentOpt();    % Agent 2 initialize options
opt3=agentOpt();    % Agent 3 initialize options
opt4=agentOpt();    % Agent 4 initialize options
opt5=agentOpt();    % Agent 5 initialize options
optC=agentOpt();    % Centralized system initialize options
opt1.name='subSystem1';
opt2.name='subSystem2';
opt3.name='subSystem3';
opt4.name='subSystem4';
opt5.name='subSystem5';
optC.name='centrSystem';
opt1.Ts=Ts;
opt2.Ts=Ts;
opt3.Ts=Ts;
opt4.Ts=Ts;
opt5.Ts=Ts;
optC.Ts=Ts;
opt1.W=W.projection(1:2);
opt2.W=W.projection(3:4);
opt3.W=W.projection(5:6);
opt4.W=W.projection(7:8);
opt5.W=W.projection(9:10);
optC.W=W;
opt1.E=E.projection(1);
opt2.E=E.projection(2);
opt3.E=E.projection(3);
opt4.E=E.projection(4);
opt5.E=E.projection(5);
optC.E=E;
opt1.Th=Th.projection(1);
opt2.Th=Th.projection(1:2);
opt3.Th=Th.projection(2:3);
opt4.Th=Th.projection(3:4);
opt5.Th=Th.projection(4);
optC.Th=Th;
opt1.th=th_real(1);
opt2.th=th_real(1:2);
opt3.th=th_real(2:3);
opt4.th=th_real(3:4);
opt5.th=th_real(4);
optC.th=th_real;

%x0=[2 1 2 -1 0 0 1 -2 0 1]';
umax=25;
xmax=[2.22 6]';
x0=repmat([0.6;0],n,1).*(rand(2*n,1).*repmat(xmax,n,1)-repmat(xmax,n,1)./2);

cAg=Agent(sys,x0,optC);

dAg1=Agent(subSys1,x0(1:2),opt1);
dAg2=Agent(subSys2,x0(3:4),opt2);
dAg3=Agent(subSys3,x0(5:6),opt3);
dAg4=Agent(subSys4,x0(7:8),opt4);
dAg5=Agent(subSys5,x0(9:10),opt5);
dAg1=dAg1.setNeighbors(dAg2,sys.A(1:2,3:4,:));
dAg2=dAg2.setNeighbors(dAg1,sys.A(3:4,1:2,:));
dAg2=dAg2.setNeighbors(dAg3,sys.A(3:4,5:6,:));
dAg3=dAg3.setNeighbors(dAg2,sys.A(5:6,3:4,:));
dAg3=dAg3.setNeighbors(dAg4,sys.A(5:6,7:8,:));
dAg4=dAg4.setNeighbors(dAg3,sys.A(7:8,5:6,:));
dAg4=dAg4.setNeighbors(dAg5,sys.A(7:8,9:10,:));
dAg5=dAg5.setNeighbors(dAg4,sys.A(9:10,7:8,:));
% Nominal Model for simple MPC design

tend=50;

r0=pw_ref(repmat([-2.2 2.2],n,1),tend,1/Ts,4,'custom');
cRef = @(t) r0(:, min(floor(t / Ts) + 1, size(r0, 2)));
dRef1 = @(t) r0(1, min(floor(t / Ts) + 1, size(r0(1,:), 2)));
dRef2 = @(t) r0(2, min(floor(t / Ts) + 1, size(r0(2,:), 2)));
dRef3 = @(t) r0(3, min(floor(t / Ts) + 1, size(r0(3,:), 2)));
dRef4 = @(t) r0(4, min(floor(t / Ts) + 1, size(r0(4,:), 2)));
dRef5 = @(t) r0(5, min(floor(t / Ts) + 1, size(r0(5,:), 2)));

%% - Controller Design

cKopt=mpcOpt();
cKopt.ref='traj';
cKopt.Ts=Tc;
cKopt.W=W;
cKopt.verbose=0;
cKopt.solver='OSQP';
cKopt.uBounds=repmat([-umax umax],n,1);
cKopt.xBounds=repmat([-xmax xmax],n,1);
cKopt.relax=true;
dKopt=dmpcOpt();
dKopt.ref='traj';
dKopt.Ts=Tc;
dKopt.W=W.projection(1:2);

dKopt.verbose=0;
dKopt.solver='OSQP';
dKopt.uBounds=[-umax umax];
dKopt.xBounds=[-xmax xmax];
dKopt.updtMode = 'fixed';
dKopt.rho=10;
dKopt.maxIter =500;
dKopt.relax=true;

% Centralized MPC

K = yMPC(sysNom,eye(2*n)*10,eye(n)*5,10,cKopt);
cAg=cAg.setController(K); 
 
% Distributed MPC

K1= yDMPC(subNomSys1,eye(2)*10,5,10,dKopt);
dAg1=dAg1.setController(K1);

K2= yDMPC(subNomSys2,eye(2)*10,5,10,dKopt);
dAg2=dAg2.setController(K2);

K3= yDMPC(subNomSys3,eye(2)*10,5,10,dKopt);
dAg3=dAg3.setController(K3);

K4= yDMPC(subNomSys4,eye(2)*10,5,10,dKopt);
dAg4=dAg4.setController(K4);

K5= yDMPC(subNomSys5,eye(2)*10,5,10,dKopt);
dAg5=dAg5.setController(K5);

distrGain('volume',K,K1,K2,K3,K4,K5);
% distrGain('performance',K,K1,K2,K3,K4,K5);
K = K.initialize;
K1 = K1.initialize;
K2 = K2.initialize;
K3 = K3.initialize;
K4 = K4.initialize;
K5 = K5.initialize;
%% -- Run&Plot Centralized MPC

cN=Network('centralized');
cN=cN.addAgent(cAg);
[x,y,u]=cN.run(tend,{cRef},'');

figure

for i=1:n*2
    if mod(i, 2) == 1
        j = 1;
    else
        j = 2;
    end
    l = ceil(i / 2);
    subplot(2,n,i)
    hold on
    plot(0:Ts:tend+Ts,x(i,1:tend/Ts+2),'r')
    plot(0:Ts:tend+Ts,xD(j,1:tend/Ts+2,l),'b')
    grid on
end


figure

for i=1:n
    
    subplot(1,n,i)
    hold on
    plot(0:Ts:tend,u(i,1:tend/Ts+1),'r')
    plot(0:Ts:tend,uD(:,1:tend/Ts+1,i),'b')
    grid on
end

I=eye(n);
figure
for i=1:n
    subplot(1,n,i)
    hold on
    plot(0:Ts:tend,y(i,1:tend/Ts+1),'r*')
    plot(0:Ts:tend,yD(:,1:tend/Ts+1,i),'b*')
    plot(0:Ts:tend,I(i,:)*cRef(0:Ts:tend),'m')
    grid on
end

%% -- Run&Plot Distributed MPC

dN=Network('distributed');  
dN=dN.addAgent(dAg1);
dN=dN.addAgent(dAg2);
dN=dN.addAgent(dAg3);
dN=dN.addAgent(dAg4);
dN=dN.addAgent(dAg5);
[xD,yD,uD]=dN.run(tend,{dRef1,dRef2,dRef3,dRef4,dRef5},'');
 
close all

figure

for i=1:n*2
    if mod(i, 2) == 1
        j = 1;
    else
        j = 2;
    end
    l = ceil(i / 2);
    subplot(2,n,i)
    hold on
    plot(0:Ts:tend+Ts,x(i,1:tend/Ts+2),'r')
    plot(0:Ts:tend+Ts,xD(j,1:tend/Ts+2,l),'b')
    grid on
end

figure

for i=1:n

    subplot(1,n,i)
    hold on
    plot(0:Ts:tend,u(i,1:tend/Ts+1),'r')
    plot(0:Ts:tend,uD(:,1:tend/Ts+1,i),'b')
    grid on
end

I=eye(n);
figure
for i=1:n
    subplot(1,n,i)
    hold on
    plot(0:Ts:tend,y(i,1:tend/Ts+1),'r*')
    plot(0:Ts:tend,yD(:,1:tend/Ts+1,i),'b*')
    plot(0:Ts:tend,I(i,:)*cRef(0:Ts:tend),'m')
    grid on
end

%save("MPC.mat")

%% Robust Adaptive MPC

clear 
clc
close all
yalmip('clear')

addpath('../Utils')
rng(754)

%% - Initialization

masses = ones(5, 1);
spring = [1.4 2.6 1.8 2.2]';
damper = [2 2 1.65 1.8]';
uncertainty.spring = ([.18 .23 .3 .15]' ./ spring);
n = 5;

[sys, Th, ~] = multiMSDgen(n, masses, spring, damper, uncertainty);

subSys1 = ss(sys.A(1:2, 1:2, :), sys.B(1:2, 1, :), sys.C(1, 1:2, :), sys.D(1, 1, :));
subSys2 = ss(sys.A(3:4, 3:4, :), sys.B(3:4, 2, :), sys.C(2, 3:4, :), sys.D(2, 2, :));
subSys3 = ss(sys.A(5:6, 5:6, :), sys.B(5:6, 3, :), sys.C(3, 5:6, :), sys.D(3, 3, :));
subSys4 = ss(sys.A(7:8, 7:8, :), sys.B(7:8, 4, :), sys.C(4, 7:8, :), sys.D(4, 4, :));
subSys5 = ss(sys.A(9:10, 9:10, :), sys.B(9:10, 5, :), sys.C(5, 9:10, :), sys.D(5, 5, :));

w_bnd=0.05;
e_bnd = 0;

Ts=0.25;
Tc=0.25;

W = Polyhedron('lb', -w_bnd * ones(2 * n, 1), 'ub', w_bnd * ones(2 * n, 1));
E = Polyhedron('lb', -e_bnd * ones(n, 1), 'ub', e_bnd * ones(n, 1));

th_real = [1.26 2.39 2.0 2.34]';
th_real = th_real - spring;

% th_real = Th.randomPoint;

sysNom = ss(sys.A(:, :, 1), sys.B(:, :, 1), sys.C(:, :, 1), sys.D(:, :, 1));
subNomSys1 = ss(sysNom.A(1:2, 1:2, :), sysNom.B(1:2, 1, :), sysNom.C(1, 1:2, :), sysNom.D(1, 1, :));
subNomSys2 = ss(sysNom.A(3:4, 3:4, :), sysNom.B(3:4, 2, :), sysNom.C(2, 3:4, :), sysNom.D(2, 2, :));
subNomSys3 = ss(sysNom.A(5:6, 5:6, :), sysNom.B(5:6, 3, :), sysNom.C(3, 5:6, :), sysNom.D(3, 3, :));
subNomSys4 = ss(sysNom.A(7:8, 7:8, :), sysNom.B(7:8, 4, :), sysNom.C(4, 7:8, :), sysNom.D(4, 4, :));
subNomSys5 = ss(sysNom.A(9:10, 9:10, :), sysNom.B(9:10, 5, :), sysNom.C(5, 9:10, :), sysNom.D(5, 5, :));
      
% Initialize 
opt1=agentOpt();    % Agent 1 initialize options
opt2=agentOpt();    % Agent 2 initialize options
opt3=agentOpt();    % Agent 3 initialize options
opt4=agentOpt();    % Agent 4 initialize options
opt5=agentOpt();    % Agent 5 initialize options
optC=agentOpt();    % Centralized system initialize options
opt1.name='subSystem1';
opt2.name='subSystem2';
opt3.name='subSystem3';
opt4.name='subSystem4';
opt5.name='subSystem5';
optC.name='centrSystem';
opt1.Ts=Ts;
opt2.Ts=Ts;
opt3.Ts=Ts;
opt4.Ts=Ts;
opt5.Ts=Ts;
optC.Ts=Ts;
opt1.W=W.projection(1:2);
opt2.W=W.projection(3:4);
opt3.W=W.projection(5:6);
opt4.W=W.projection(7:8);
opt5.W=W.projection(9:10);
optC.W=W;
opt1.E=E.projection(1);
opt2.E=E.projection(2);
opt3.E=E.projection(3);
opt4.E=E.projection(4);
opt5.E=E.projection(5);
optC.E=E;
opt1.Th=Th.projection(1);
opt2.Th=Th.projection(1:2);
opt3.Th=Th.projection(2:3);
opt4.Th=Th.projection(3:4);
opt5.Th=Th.projection(4);
optC.Th=Th;
opt1.th=th_real(1);
opt2.th=th_real(1:2);
opt3.th=th_real(2:3);
opt4.th=th_real(3:4);
opt5.th=th_real(4);
optC.th=th_real;

umax=25;
xmax=[2.22 6]';
x0=repmat([0.6;0],n,1).*(rand(2*n,1).*repmat(xmax,n,1)-repmat(xmax,n,1)./2);

cAg=Agent(sys,x0,optC);

dAg1=Agent(subSys1,x0(1:2),opt1);
dAg2=Agent(subSys2,x0(3:4),opt2);
dAg3=Agent(subSys3,x0(5:6),opt3);
dAg4=Agent(subSys4,x0(7:8),opt4);
dAg5=Agent(subSys5,x0(9:10),opt5);
dAg1=dAg1.setNeighbors(dAg2,sys.A(1:2,3:4,:));
dAg2=dAg2.setNeighbors(dAg1,sys.A(3:4,1:2,:));
dAg2=dAg2.setNeighbors(dAg3,sys.A(3:4,5:6,:));
dAg3=dAg3.setNeighbors(dAg2,sys.A(5:6,3:4,:));
dAg3=dAg3.setNeighbors(dAg4,sys.A(5:6,7:8,:));
dAg4=dAg4.setNeighbors(dAg3,sys.A(7:8,5:6,:));
dAg4=dAg4.setNeighbors(dAg5,sys.A(7:8,9:10,:));
dAg5=dAg5.setNeighbors(dAg4,sys.A(9:10,7:8,:));
% Nominal Model for simple MPC design

tend=50;

r0=pw_ref(repmat([-2.2 2.2],n,1),tend,1/Ts,4,'custom');
cRef = @(t) r0(:, min(floor(t / Ts) + 1, size(r0, 2)));
dRef1 = @(t) r0(1, min(floor(t / Ts) + 1, size(r0(1,:), 2)));
dRef2 = @(t) r0(2, min(floor(t / Ts) + 1, size(r0(2,:), 2)));
dRef3 = @(t) r0(3, min(floor(t / Ts) + 1, size(r0(3,:), 2)));
dRef4 = @(t) r0(4, min(floor(t / Ts) + 1, size(r0(4,:), 2)));
dRef5 = @(t) r0(5, min(floor(t / Ts) + 1, size(r0(5,:), 2)));

%% - Controller Design

cKopt=mpcOpt();
cKopt.ref='traj';
cKopt.Ts=Tc;
cKopt.W=W;
cKopt.verbose=0;
cKopt.solver='OSQP';
cKopt.uBounds=repmat([-umax umax],n,1);
cKopt.xBounds=repmat([-xmax xmax],n,1);

dKopt=dmpcOpt();
dKopt.ref='traj';
dKopt.Ts=Tc;
dKopt.W=W.projection(1:2);

dKopt.verbose=0;
dKopt.solver='OSQP';
dKopt.uBounds=[-umax umax];
dKopt.xBounds=[-xmax xmax];
dKopt.updtMode = 'fixed';
dKopt.rho=10;
dKopt.maxIter =500;

K = yRAMPC(sys,eye(2*n)*10,eye(n)*5,10,cKopt);
cAg=cAg.setController(K);


% Distributed MPC

K1= yDRAMPC(subSys1,eye(2)*10,5,10,dKopt);
dAg1=dAg1.setController(K1);

dKopt.W=W.projection(3:4);
K2= yDRAMPC(subSys2,eye(2)*10,5,10,dKopt);
dAg2=dAg2.setController(K2);

dKopt.W=W.projection(5:6);
K3= yDRAMPC(subSys3,eye(2)*10,5,10,dKopt);
dAg3=dAg3.setController(K3);

dKopt.W=W.projection(7:8);
K4= yDRAMPC(subSys4,eye(2)*10,5,10,dKopt);
dAg4=dAg4.setController(K4);

dKopt.W=W.projection(9:10);
K5= yDRAMPC(subSys5,eye(2)*10,5,10,dKopt);
dAg5=dAg5.setController(K5);


K = K.initUnc(Th,[]);
K1 = K1.initUnc(opt1.Th,[]);
K2 = K2.initUnc(opt2.Th,[]);
K3 = K3.initUnc(opt3.Th,[]);
K4 = K4.initUnc(opt4.Th,[]);
K5 = K5.initUnc(opt5.Th,[]);

distrGain('volume',K,K1,K2,K3,K4,K5);

K = K.initialize();
K1 = K1.initialize();
K2 = K2.initialize();
K3 = K3.initialize();
K4 = K4.initialize();
K5 = K5.initialize();

%% -- Run&Plot Centralized RAMPC

cN=Network('centralized');
cN=cN.addAgent(cAg);
[x,y,u]=cN.run(tend,{cRef},'');

figure

for i=1:n*2
    if mod(i, 2) == 1
        j = 1;
    else
        j = 2;
    end
    l = ceil(i / 2);
    subplot(2,n,i)
    hold on
    plot(0:Ts:tend+Ts,x(i,1:tend/Ts+2),'r')
%     plot(0:Ts:tend+Ts,xD(j,1:tend/Ts+2,l),'b')
    grid on
end


figure

for i=1:n
    
    subplot(1,n,i)
    hold on
    plot(0:Ts:tend,u(i,1:tend/Ts+1),'r')
%     plot(0:Ts:tend,uD(:,1:tend/Ts+1,i),'b')
    grid on
end

I=eye(n);
figure
for i=1:n
    subplot(1,n,i)
    hold on
    plot(0:Ts:tend,y(i,1:tend/Ts+1),'r*')
%     plot(0:Ts:tend,yD(:,1:tend/Ts+1,i),'b*')
    plot(0:Ts:tend,I(i,:)*cRef(0:Ts:tend),'m')
    grid on
end

%% -- Run&Plot Distributed RAMPC

dN=Network('distributed');  
dN=dN.addAgent(dAg1);
dN=dN.addAgent(dAg2);
dN=dN.addAgent(dAg3);
dN=dN.addAgent(dAg4);
dN=dN.addAgent(dAg5);
[xD,yD,uD]=dN.run(tend,{dRef1,dRef2,dRef3,dRef4,dRef5},'');
 
close all

figure

for i=1:n*2
    if mod(i, 2) == 1
        j = 1;
    else
        j = 2;
    end
    l = ceil(i / 2);
    subplot(2,n,i)
    hold on
    plot(0:Ts:tend+Ts,x(i,1:tend/Ts+2),'r')
    plot(0:Ts:tend+Ts,xD(j,1:tend/Ts+2,l),'b')
    grid on
end

figure

for i=1:n

    subplot(1,n,i)
    hold on
    plot(0:Ts:tend,u(i,1:tend/Ts+1),'r')
    plot(0:Ts:tend,uD(:,1:tend/Ts+1,i),'b')
    grid on
end

I=eye(n);
figure
for i=1:n
    subplot(1,n,i)
    hold on
    plot(0:Ts:tend,y(i,1:tend/Ts+1),'r*')
    plot(0:Ts:tend,yD(:,1:tend/Ts+1,i),'b*')
    plot(0:Ts:tend,I(i,:)*cRef(0:Ts:tend),'m')
    grid on
end
%save("RAMPC.mat")
% LogFun({'yRAMPC_318.csv'},{'yDRAMPC_313.csv','yDRAMPC_359.csv','yDRAMPC_623.csv','yDRAMPC_630.csv','yDRAMPC_655.csv'},Ts,tend)
%movefile *.mat Logs\Complete_5ag\PiecewiseT
%movefile *.csv Logs\Complete_5ag\PiecewiseT

