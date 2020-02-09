function [bic,cost,score,history] = qpsob(nPop,data,lamda,miu,omega)

% INPUT:
%   data    n*m
%   lamda   the weight of gene Volume
%   miu     the weight of condition Volume
% OUTPUT:
%   xmin    : best solution found
%   fmin    : function value at the best solution, f(xmin)
%   histout : record of function evaluations and fitness value by iteration

%% Problem Definiton
data = choseData(data);
n = size(data,1);
m = size(data,2);
D = n+m;        % Number of Unknown (Decision) Variables


lb = 0;	        % Lower Bound of Decision Variables
ub = 1;         % Upper Bound of Decision Variables

ShowIterInfo =1;

%% Parameters of QPSO

maxit = 300;          % Maximum Number of Iterations
early_stopping_cnt = 0;
early_stopping_maxcnt = 40;
% nPop = 20;            % Population Size (Swarm Size)
w1 = 0.5;
w2 = 1.0;

c1 = 1.5;
c2 = 1.5;

c2bT = 0.5;         % decide 0 or 1, threshold
% Initializing solution
x = unifrnd(lb,ub,[nPop,D]);

% Evaluate initial population
pbest = x;

costFun = @calc_fit4;


f_x = 500*ones(nPop,1);

for i = 1:nPop
    bic = conti2bit(x(i,:),c2bT);
    f_x(i) = costFun(bic,data,lamda,miu,omega);
end

f_pbest = f_x;

[~,g] = min(f_pbest);
gbest = pbest(g,:);
f_gbest = f_pbest(g);


history = zeros(maxit,1);

f_gbest_old = f_gbest;
for  it = 1:maxit 

    alpha = (w2 - w1) * (maxit - it)/maxit + w1;
    mbest = sum(pbest)/nPop;

    for i = 1:nPop

        fi = rand(1,D);

        p = (c1*fi.*pbest(i,:) + c2*(1-fi).*gbest)/(c1 + c2);

        u = rand(1,D);
        
        b = alpha*abs(x(i,:) - mbest);
        v = log(1./u);

        if rand < 0.5
            x(i,:) = p + b .* v;
        else
            x(i,:) = p - b .* v;
        end
        
        % Keeping bounds
        x(i,:) = max(x(i,:),lb);
        x(i,:) = min(x(i,:),ub);

        bic = conti2bit(x(i,:),c2bT);
        f_x(i) = costFun(bic,data,lamda,miu,omega);
        

        if f_x(i) < f_pbest(i)
            pbest(i,:) = x(i,:);
            f_pbest(i) = f_x(i);
        end

        if f_pbest(i) < f_gbest
            gbest = pbest(i,:);
            f_gbest = f_pbest(i);
        end

    end
    % early stopping
    change = f_gbest_old - f_gbest;
    if change < f_gbest_old/10000
        early_stopping_cnt = early_stopping_cnt + 1;
    else
        early_stopping_cnt = 0;
    end
    % Display Iteration Information
    if ShowIterInfo
        disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(f_gbest)]);
    end
    history(it) = f_gbest;
    f_gbest_old = f_gbest;

    if early_stopping_cnt > early_stopping_maxcnt
        disp('qpsob early stoping');
        break
    end
end % end of iteration
if it < maxit
    for i =it:maxit
        history(i) = f_gbest;
    end
end
bic = conti2bit(gbest,c2bT);
cost = f_gbest;

% costs = zeros(nPop,1);
% bics = zeros(nPop,D);
% for i = 1:nPop
%     costs(i) = f_x(i);
%     bics(i,:) = conti2bit(x(i,:),c2bT);
% end
score = calc_bench(bic,data);


end