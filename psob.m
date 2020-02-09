function [bic,cost,score,history] = psob(nPop,data,lamda,miu,omega)
%   data    n*m
%   lamda   the weight of gene Volume
%   miu     the weight of condition Volume
    %% Problem Definiton
    data = choseData(data);
    n = size(data,1);
    m = size(data,2);
    nVar = n+m;        % Number of Unknown (Decision) Variables

    VarSize = [1 nVar];         % Matrix Size of Decision Variables

    VarMin = 0;	        % Lower Bound of Decision Variables
    VarMax = 1;         % Upper Bound of Decision Variables
    c2bT = 0.5;         % decide 0 or 1 threshold
    costFun = @calc_fit4;
    %% Parameters of PSO

    MaxIt = 300;   % Maximum Number of Iterations
    early_stopping_cnt = 0;
    early_stopping_maxcnt = 40;

    % nPop = 20;     % Population Size (Swarm Size)

    w = 0.6;           % Intertia Coefficient
    wdamp = 1;         % Damping Ratio of Inertia Coefficient
    c1 = 1.80;         % Personal Acceleration Coefficient
    c2 = 1.00;         % Social Acceleration Coefficient

    % The Flag for Showing Iteration Information
    ShowIterInfo = 1;    

    MaxVelocity = 0.2*(VarMax-VarMin);
    MinVelocity = -MaxVelocity;
    
    %% Initialization

    % The Particle Template
    empty_particle.Position = [];
    empty_particle.Velocity = [];
    empty_particle.Cost = [];
    empty_particle.Best.Position = [];
    empty_particle.Best.Cost = [];

    % Create Population Array
    particle = repmat(empty_particle, nPop, 1);

    % Initialize Global Best
    GlobalBest.Cost = inf;

    % Initialize Population Members
    for i=1:nPop

        % Generate Random Solution
        particle(i).Position = unifrnd(VarMin, VarMax, VarSize);

        % Initialize Velocity
        particle(i).Velocity = zeros(VarSize);

        % Evaluation
        bic = conti2bit(particle(i).Position,c2bT);
        particle(i).Cost = costFun(bic,data,lamda,miu,omega);

        % Update the Personal Best
        particle(i).Best.Position = particle(i).Position;
        particle(i).Best.Cost = particle(i).Cost;
        % Update Global Best
        if particle(i).Best.Cost < GlobalBest.Cost
            GlobalBest = particle(i).Best;
        end
    end

    % Array to Hold Best Cost Value on Each Iteration
    history = zeros(MaxIt, 1);


    GBestCostOld = GlobalBest.Cost;
%% Main Loop of PSO

    for it=1:MaxIt

        for i=1:nPop

            % Update Velocity
            particle(i).Velocity = w*particle(i).Velocity ...
                + c1*rand(VarSize).*(particle(i).Best.Position - particle(i).Position) ...
                + c2*rand(VarSize).*(GlobalBest.Position - particle(i).Position);

            % Apply Velocity Limits
            particle(i).Velocity = max(particle(i).Velocity, MinVelocity);
            particle(i).Velocity = min(particle(i).Velocity, MaxVelocity);
            
            % Update Position
            particle(i).Position = particle(i).Position + particle(i).Velocity;
            
            % Apply Lower and Upper Bound Limits
              particle(i).Position = max(particle(i).Position, VarMin);
              particle(i).Position = min(particle(i).Position, VarMax);

            % Evaluation
            bic = conti2bit(particle(i).Position,c2bT);
            particle(i).Cost = costFun(bic,data,lamda,miu,omega);

            % Update Personal Best
            if particle(i).Cost < particle(i).Best.Cost

                particle(i).Best.Position = particle(i).Position;
                particle(i).Best.Cost = particle(i).Cost;

                % Update Global Best
                if particle(i).Best.Cost < GlobalBest.Cost
                    GlobalBest = particle(i).Best;
                end            

            end
        end

        % early stopping
        change =  GBestCostOld - GlobalBest.Cost;
        if change < GBestCostOld/10000
            early_stopping_cnt = early_stopping_cnt + 1;
        else
            early_stopping_cnt = 0;
        end
        % Store the Best Cost Value
        history(it) = GlobalBest.Cost;
        GBestCostOld = GlobalBest.Cost;

        % Display Iteration Information
        if ShowIterInfo
            disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(history(it))]);
        end

        % Damping Inertia Coefficient
        w = w * wdamp;
        if early_stopping_cnt > early_stopping_maxcnt
            disp('psob early stoping');
            break
        end
    end % end of iteration
    
    if it < MaxIt
        for i =it:MaxIt
            history(i) = GlobalBest.Cost;
        end
    end
    bic = conti2bit(GlobalBest.Position,c2bT);
    cost = GlobalBest.Cost;
    score = calc_bench(bic,data);
end