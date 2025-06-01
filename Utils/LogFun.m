function LogFun(mpcFiles, dmpcFiles,Ts,Tend)
    % LogReadUtility reads MPC and DMPC log files and plots the required information.
    % Inputs:
    %   mpcFiles  - Cell array of strings containing the filenames of MPC CSV files.
    %   dmpcFiles - Cell array of strings containing the filenames of DMPC CSV files.

    % Read MPC files
    t=0:Ts:Tend;
    t=t';
    mpcData = cell(size(mpcFiles));
    for i = 1:length(mpcFiles)
        mpcData{i} = readtable(mpcFiles{i},"Delimiter",',');
    end

    % Read DMPC files
    dmpcData = cell(size(dmpcFiles));
    for i = 1:length(dmpcFiles)
        dmpcData{i} = readtable(dmpcFiles{i},"Delimiter",',');
    end

    % Find the index of the first iteration in the first DMPC file
    idx = find(dmpcData{1}.Iteration == 1);
%     idx = idx(1:min(length(idx),length(t)));
    %% Cost compare
    ADMMcost = zeros(height(dmpcData{1}), 1);
    Consensus = zeros(height(dmpcData{1}), 1);
    for i = 1:length(dmpcData)
        ADMMcost = ADMMcost + dmpcData{i}.kCost;
        Consensus= Consensus+dmpcData{i}.consCost;
    end
    ADMMcost = [ADMMcost(idx(2:end)-1,:); ADMMcost(end)];
    Consensus = [Consensus(idx(2:end)-1,:); Consensus(end)];
    Centcost = 0.5.*mpcData{1}.OpCost;
    e = zeros(length(ADMMcost), 1);
    for i = 1:length(ADMMcost)
        e(i) = rms(Centcost(i) - ADMMcost(i));  
    end

    figure
    subplot(2,2,1)
    plot(t,Centcost(1:length(ADMMcost)),t(1:length(ADMMcost)),ADMMcost)
    grid on
    legend('Centralized','Distributed')
    title('Main Problem Cost')
    subplot(2,2,3)
    plot(t(1:length(Consensus)),sqrt(Consensus.^2))
    grid on
    title('ADMM summed Costs')
    subplot(2,2,[2 4])
    plot(e)
    grid on
    title('RMS Error')

    %% Time compare
    time = zeros(height(dmpcData{1}), length(dmpcData));
    for i = 1:length(dmpcData)
        time(:, i) = dmpcData{i}.ElapsTime_s_;
    end
    dTime = mean([time(idx(2:end)-1,:); time(end,:)], 2);
    cTime = mpcData{1}.OpTime;
    Iter=[];
    for i=1:length(dmpcData)
        Iter=[Iter,[dmpcData{i}.Iteration(idx(2:end)-1,:);dmpcData{i}.Iteration(end,1)]];
    end
    Iter = mean(Iter,2);
    figure
    subplot(2,1,1)
    plot(t(1:length(dTime)),dTime,t,cTime)
    set(gca, 'YScale', 'log');
    grid on
    legend('Distributed','Centralized')
    title('Optimisation Time')
    subplot(2,1,2)
    plot(t(1:length(Iter)),Iter)
    grid on
    title('ADMM Iterations')

    %% Residuals
%     figure
%     for i = 1:length(dmpcData)
%         subplot(2,1,1)
%         plot(t,[dmpcData{i}.PrimalRes(idx(2:end)-1,:);dmpcData{i}.PrimalRes(end,:)]);
%         hold on
%         grid on
%            set(gca, 'YScale', 'log');
%         subplot(2,1,2)
%         plot(t,[dmpcData{i}.DualRes(idx(2:end)-1,:);dmpcData{i}.DualRes(end,:)]);
%         hold on
%         grid on
%            set(gca, 'YScale', 'log');
%     end
end