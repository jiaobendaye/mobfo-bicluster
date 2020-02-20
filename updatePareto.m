function Pset_new  = updatePareto(Pset, news)
    %input:
    %     Pset         old_size*(n+m) 0/1 bits
    %     news            N*(n+m)
    % output:
    %     Pset_new     new_size*(n+m)
    global Pset_size weight
    old_size = size(Pset,1);
    %remove duplications
    tmp = unique([Pset; news], 'rows');
    if size(tmp,1) == old_size && size(tmp,1) <= Pset_size
        Pset_new = Pset;
        return
    end
    scores = calc_scores(tmp);
    domained_cnt  = domain_fit(scores);
    pareto_set = tmp(domained_cnt == 0,:);
    if size(pareto_set,1)<=Pset_size
        Pset_new = pareto_set;
    else
        scores = abs(scores(domained_cnt == 0,:));
        total = sum(scores);
        ave = 1 / size(pareto_set,1);
        percent = scores ./ total;
        difference = percent - ave;
        final = difference *  weight;
        %final 越大越好, 取最大的Pset_size个
        [~, I] = sort(final, 'desc');
        Pset_new = pareto_set(I(1:Pset_size),:);
    end
end
    