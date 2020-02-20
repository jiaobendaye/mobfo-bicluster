function Bool = is_better(old, new)
    % input:
    %     old       [msr, -gv, -cv]
    %     new       [msr, -gv, -cv]
    %     lambda, miu, meage    weights
    % output:
    %     if new is better than old, return 1, else 0;
    global lambda miu omega
    wei = [-lambda, miu, omega]';
    if old == new
        Bool = 0;
        return
    end
    if is_domain(old,new)
        Bool = 1;
    elseif is_domain(new,old)
        Bool = 0;
    else
        % normalization 
        both = abs([old; new]);
        total = sum(both);
        percentage = both ./total;
        difference = percentage(2,:) - percentage(1,:);
        final = difference * wei;
        Bool = final > 0;
    end
    
end