function [domained_cnt, best, paretoSet] = pickMinMSR(Popu)
    %input
    % Popu        (n+m)*nPop continue
    %output
    %   domained_cnt    nPop*1
    %   best  min msr solution in pso-pareto set, (n＋m)*1
    %   paretoSet       size*(n+m)  0/1bits 
    bics = conti2bit(Popu)';
    scores = calc_scores(bics);
    domained_cnt  = domain_fit(scores);

    min_cnt =  min(domained_cnt);
    index = domained_cnt == min_cnt;

    psoPareto = Popu(:,index);
    psoScores = scores(index,:);
    % sort by msr
    [~, sortind] = sortrows(psoScores);
    best = psoPareto(:,sortind(1));
    if min_cnt == 0
        paretoSet = conti2bit(psoPareto)'; 
    else
        paretoSet = [];
    end
end