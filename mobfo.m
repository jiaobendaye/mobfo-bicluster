%*********************细菌觅食算法**********************
%%%%%%%%%%%%%%%%%%%-----BFA算法-----%%%%%%%%%%%
clear;
clc;

global data omega lambda miu Pset_size c2bTh Lb Ub
data = choseData(1);
lambda = 0.7;
miu = 0.1;
omega =0.2;
Pset_size = 20;
n = size(data,1);
m = size(data,2);
c2bTh = 0.5;
%-----(1)初始化参数-----
p = n+m;    % 搜索范围的维度
Ub = ones(p,1);   %上界
Lb = zeros(p,1);  %下界
s = 100;   % 细菌的个数
Nc = 50;  % 趋化的次数
Ns = 4;   % 趋化操作中单向运动的最大步数
C(:,1) = 0.001*ones(s,1);  % 翻转选定方向后，单个细菌前进的步长
Nre = 4;    % 复制操作步骤数
Ned = 2;    % 驱散(迁移)操作数
Sr = s/2;   % 每代复制（分裂）数
Ped = 0.25; % 细菌驱散(迁移)概率
beta=3/2;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
for i = 1:s     % 产生初始细菌个体的位置
    P(:,i,1) = rand(p,1);
end

[~, best, pareto_set] = pickMinMSR(P(:,:,1));
%------------------细菌趋药性算法循环开始---------------------
%-----(2)驱散(迁移)操作开始-----
for l = 1:Ned        
    %-----(3)复制操作开始-----
    for k = 1:Nre   
        %-----(4)趋化操作(翻转或游动)开始-----
        for j = 1:Nc  
            disp(['j: ',num2str(j),' k: ',num2str(k),' l: ',num2str(l)])
            %-----(4.1)对每一个细菌分别进行以下操作-----

            for i = 1:s
                %-----(4.2)计算函数J(i,j,k,l)，表示第i个细菌在第l次驱散第k次
                %----------复制第j次趋化时的适应度值-----
                J(i,j,:) = myCost(P(:,i,j));
                %-----(4.4)保存细菌目前的适应度值，直到找到更好的适应度值取代之-----
                Jlast = J(i,j,:);
                %-----(4.5)翻转，产生一个随机向量PHI(i),代表翻转后细菌的方向-----
                Delta(:,i) = (2*round(rand(p,1))-1).*rand(p,1);
                % PHI表示翻转后选择的一个随机方向上前进
                PHI = Delta(:,i)/sqrt(Delta(:,i)'*Delta(:,i));
                % %-----(4.6)移动，向着翻转后细菌的方向移动一个步长，并且改变细菌的位置-----
                % P(:,i,j+1) = P(:,i,j) + C(i)*PHI;
                %levy 飞行
                % s = P(:,i,j);
                u=randn(size(P(:,i,j)))*sigma;
                v=randn(size(P(:,i,j)));
                step=u./abs(v).^(1/beta);
                stepsize=0.01*step.*(P(:,i,j)-best);
                % P(:,i,j+1) = s + stepsize.*randn(size(P(:,i,j)));
                P(:,i,j+1) = simplebounds(P(:,i,j) + stepsize.*randn(size(P(:,i,j))));

                %-----(4.7)计算细菌当前位置的适应度值-----
                J(i,j+1,:) = myCost(P(:,i,j+1));
                %-----(4.8)游动-----
                m = 0; % 给游动长度计数器赋初始值
                while(m < Ns) % 未达到游动的最大长度，则循环
                    m = m + 1;
                    % 新位置的适应度值是否更好？如果更好，将新位置的适应度值
                    % 存储为细菌i目前最好的适应度值
                    % if(J(i,j+1,k,l) < Jlast)
                    if(is_better(Jlast, J(i,j+1,:), lambda, miu, omega))
                        Jlast = J(i,j+1,:);  %保存更好的适应度值
                        % 在该随机方向上继续游动步长单位,修改细菌位置
                        P(:,i,j+1) =simplebounds(P(:,i,j+1) + stepsize.*randn(size(P(:,i,j+1))));
                        % 重新计算新位置上的适应度值
                        J(i,j+1,:) = myCost(P(:,i,j+1));
                    else
                        % 否则，结束此次游动
                        m = Ns;
                    end
                end
                J(i,j,:) = Jlast; % 更新趋化操作后的适应度值
            end  % 如果i<N，进入下一个细菌的趋化，i=i+1
            min_J = squeeze(min(J(:,j,:)))'
            [~,best,~] = pickMinMSR(P(:,:,j+1));
            %-----(5)如果j<Nc，此时细菌还处于活跃状态，进行下一次趋化，j=j+1-----
        end
        %----------------下面进行复制操作----------------
        %-----(6)复制-----
        %-----(6.1)根据所给的k和l的值，将每个细菌的适应度值按升序排序-----
        [Jhealth, ~, news] = pickMinMSR(P(:,:,Nc+1));
        pareto_set = updatePareto(pareto_set, news);
        % Jhealth = sum(J(:,:,k,l),2);  % 给每个细菌设置健康函数值
        [Jhealth,sortind] = sort(Jhealth); % 按健康函数值升序排列函数
        P(:,:,1) = P(:,sortind,Nc+1);
        %-----(6.2)将代价小的一半细菌分裂成两个，代价大的一半细菌死亡-----
        for i = 1:Sr
            % 健康值较差的Sr个细菌死去，Sr个细菌分裂成两个子细菌，保持个体总数的s一致性
            P(:,i+Sr,1) = P(:,i,1);
        end
    %-----(7)如果k<Nre，转到(3)，进行下一代细菌的趋化-----
    end
    %-----(8)趋散，对于每个细菌都以Ped的概率进行驱散，但是驱散的细菌群体的总数
    %--------保持不变，一个细菌被驱散后，将被随机重新放置到一个新的位置-----
    for m = 1:s
        % 产生随机数，如果既定概率大于该随机数，细菌i灭亡，随机产生新的细菌i
        if(Ped > rand)
            P(:,m,1) = rand(p,1);
        end
    end
end  % 如果l<Ned，转到(2)，否则结束
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%报告
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pareto_scores = calc_scores(pareto_set, data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cost = myCost(x)
    %input:
    %   x           (n+m)*nPop
    %output:
    %   cost        nPop*[msr, -gv, -cv]
    global data c2bTh
    bic = conti2bit(x', c2bTh);
    cost = calc_scores(bic, data);
end
function s=simplebounds(s)
    global Lb Ub
    % Apply the lower bound
    ns_tmp=s;
    I=ns_tmp<Lb;
    ns_tmp(I)=Lb(I);
    
    % Apply the upper bounds 
    J=ns_tmp>Ub;
    ns_tmp(J)=Ub(J);
    % Update this new move 
    s=ns_tmp;
end
function [domained_cnt, best, paretoSet] = pickMinMSR(Popu)
    %input
    % Popu        (n+m)*nPop continue
    %output
    %   domained_cnt    nPop*1
    %   best  min msr solution in pso-pareto set, (n＋m)*1
    %   paretoSet       size*(n+m)  0/1bits 
    global data c2bTh
    bics = conti2bit(Popu, c2bTh)';
    scores = calc_scores(bics, data);
    domained_cnt  = domain_fit(scores);

    min_cnt =  min(domained_cnt);
    index = domained_cnt == min_cnt;

    psoPareto = Popu(:,index);
    psoScores = scores(index,:);
    % sort by msr
    [~, sortind] = sortrows(psoScores);
    best = psoPareto(:,sortind(1));
    if min_cnt == 0
        paretoSet = conti2bit(psoPareto, c2bTh)'; 
    else
        paretoSet = [];
    end
end
